"""

 Load or calculate morphometrics

 To calculate morphometrics, we do PCA on some lobula landmark neurons (like LTs) and then
 fit polynomial surfaces to get the idea of where the layers are

"""
## Packages
from neuprint import Client
import pandas as pd
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA
## My own modules
import modules.getbodyids as getbodyids
import modules.getsynapses as getsynapses
import modules.visualize as visualize
import modules.utility as utility

# Return morphology (depth + spread) dataframes -- either saved or new
def getmorphology(**kwargs):
    # You can specify a substring of filename to automatically select saved
    # morphology data -- provide nonexistent name to calculate something anew
    if 'filename' in kwargs:
        filename = kwargs.get('filename')
    else:
        filename = ''

    if 'landmarkname' in kwargs:
        landmarkname = kwargs.get('landmarkname')
    else:
        landmarkname = ''

    # sometimes you might want to try different binning parameters for morphology
    # with the same bodyId list. Use this flag for such cases
    if 'recaluculateMorphology' in kwargs:
        recalcFlag = kwargs.get('recalculateMorphology')
    else:
        recalcFlag = 0

    # just making explicit what is being called...
    print('Running getmorphology...')

    # First, list the existing morphology matrices
    # assumption is that depth/spread matrices are generated as pairs
    depthlist = glob.glob('.\\data\\depth\\depth_*.csv')

    # take the ones whose name contains the specified filename
    newlist = [file for file in depthlist if filename in file and landmarkname in file]

    if newlist and not recalcFlag:
        if len(newlist)==1:
            ind = 0
        else:
            print('Found multiple dataset that matches the filename provided')
            utility.print_indexed(newlist)
            ind = int(input('Select which one you want to use : '))

        depth = pd.read_csv(newlist[ind])
        _, depth_filename = os.path.split(newlist[ind])
        spread = pd.read_csv('./data/spread/spread'+depth_filename[5:])
    else:
        # if not calculate anew
        print('Morphology dataset with the specified filename does not exist.')
        depth, spread, depth_filename = calcmorphology(**kwargs)
    return depth, spread, depth_filename

# Calculate morphology given bodyId list
def calcmorphology(**kwargs):
    # load relevant kwarg
    # Synapse type to use
    if 'synapseType' in kwargs:
        synapseType = kwargs.get('synapseType')
    else:
        synapseType = 'pre'

    # histogram related
    if 'minD' in kwargs:
        minD = kwargs.get('minD')
    if 'maxD' in kwargs:
        maxD = kwargs.get('maxD')
    if 'binSize' in kwargs:
        binSize = kwargs.get('binSize')

    # just making explicit what is being called...
    print('Running calcmorphology...')

    # Get a bodyId list you want to use
    bodyidlist, filename = getbodyids.getbodyids(**kwargs)

    # load PC coefficients and surface model
    pca, modelcoeff, landmarkname = loadlobulamodel(**kwargs)

    # prepare output structure (we build this by appending columns to bodyidlist)
    depth = bodyidlist.copy()
    spread = bodyidlist.copy()

    # define depth histogram edges if not specified already
    if not 'minD' in locals():
        minD    = float(input('Enter minimum depth (in microns): '))
    if not 'maxD' in locals():
        maxD    = float(input('Enter maximum depth (in microns): '))
    if not 'binSize' in locals():
        binSize = float(input('Enter depth bin size (in microns): '))
    binEdges = np.arange(minD,maxD+binSize,binSize)

    # add columns
    for b in range(len(binEdges)-1):
        depth['bin'+str(b)] = np.zeros(len(spread))

    spread['SD1'] = np.zeros(len(spread))
    spread['SD2'] = np.zeros(len(spread))
    spread['SD3'] = np.zeros(len(spread))

    # go through the bodyid list and load synapses
    print('Calculating morphological metrics. This could take a while...')
    for thisId in bodyidlist['bodyId']:
        # load synapses
        synapses = getsynapses.getsynapses(thisId,synapseType)
        # calculate PCs and depth
        rawdepth, PCs = utility.calcrawdepth(pca, modelcoeff, synapses)
        # calculate depth histogram
        for b in range(len(binEdges)-1):
            depth.loc[depth['bodyId']==thisId,'bin'+str(b)] = np.sum(np.logical_and(rawdepth<binEdges[b+1], rawdepth>binEdges[b]))
        spread.loc[spread['bodyId']==thisId,'SD1'] = np.std(PCs[:,0])
        spread.loc[spread['bodyId']==thisId,'SD2'] = np.std(PCs[:,1])
        spread.loc[spread['bodyId']==thisId,'SD3'] = np.std(PCs[:,2])

    # save as csv
    filename_postfix = landmarkname+'_'+synapseType+'_minD'+str(minD)+'_maxD'+str(maxD)+'_bin'+str(binSize)+'_'+filename

    depth.to_csv('./data/depth/depth_'+filename_postfix)
    spread.to_csv('./data/spread/spread_'+filename_postfix)
    return depth, spread, 'depth_'+filename_postfix


# Load (or download) saved synapses, run PCA, and fit a surface
def loadlobulamodel(**kwargs):
    # load relevant kwarg
    if 'landmarkname' in kwargs:
        landmarkname = kwargs.get('landmarkname')
    else:
        landmarkname = ''

    if 'showModel' in kwargs:
        showModel = kwargs.get('showModel')
    else:
        showModel = 1

    # just making explicit what is being called...
    print('Running loadlobulamodel...')

    # See existing "landmark" synapse directories
    landmarklist = glob.glob('./data/landmark/*.csv')
    newlist = [landmark for landmark in landmarklist if landmarkname in landmark]

    if newlist:
        if len(newlist)==1:
            ind = 0
        else:
            print('Found multiple landmark data matching the landmarkname provided')
            utility.print_indexed(landmarklist)
            ind = int(input('Enter which one to use (type -1 if you want to try new cell): '))

        landmark = pd.read_csv(newlist[ind])
        _, filename = os.path.split(newlist[ind])
        landmarkname = filename[:-4] # keep the name of the cell

    # if there is nothing saved or if you want to try a new landmark
    else:
        # in case nothing was saved and nothing was specified, ask here
        if not landmarkname:
            landmarkname = input('Enter the name of cell type you want to use as a landmark: ')

        # Type in the cell type to use
        print('Downloading '+landmarkname+' synapses...')

        # Connect to the neuPrint server
        f = open("authtoken","r")
        tokenstr = f.read()
        c = Client('neuprint.janelia.org', dataset='hemibrain:v1.2.1', token=tokenstr)
        c.fetch_version()

        # define query
        q = """\
            MATCH (a:Neuron)-[:Contains]->(:SynapseSet)-[:Contains]->(s:Synapse)
            WHERE a.type='%s' AND s.type='post' AND s.`LO(R)`
            RETURN DISTINCT s.location.x as x, s.location.y as y, s.location.z as z
            """ % landmarkname

        # run the query
        landmark = c.fetch_custom(q)
        landmark.to_csv('./data/landmark/'+landmarkname+'.csv')

    # Now run PCA on x/y/z in the "landmark"
    # Also convert 8 nm px to microns unit
    X = landmark["x"].to_numpy().flatten()*8/1000
    Y = landmark["y"].to_numpy().flatten()*8/1000
    Z = landmark["z"].to_numpy().flatten()*8/1000
    XYZ = np.array([X,Y,Z]).T
    pca = PCA(n_components=3)
    PCs = pca.fit_transform(XYZ)

    # Do fitting of quadric model
    # This model can be improved
    # fit quadric model
    PC1 = PCs[:,0].flatten()
    PC2 = PCs[:,1].flatten()
    PC3 = PCs[:,2].flatten()

    # predictors
    A = np.array([PC1*0+1, PC1, PC2, PC1**2, PC2**2, PC1*PC2]).T

    # fitting
    modelcoeff, r, rank, s = np.linalg.lstsq(A,PC3,rcond=None)


    # show goodness of fit
    r2 = 1 - r / np.sum(PC3**2)
    print('R2 of the lobula model was: ',r2)

    # visualize this
    if showModel:
        visualize.plotquadricandscatter(pca, modelcoeff, landmark)

    return pca, modelcoeff, landmarkname
