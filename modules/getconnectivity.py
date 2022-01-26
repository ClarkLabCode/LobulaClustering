"""

 Load or calculate connectivity matrix 
 
"""

## Packages
from neuprint import Client
import pandas as pd
import numpy as np
import os
import glob
import modules.getbodyids as getbodyids
import modules.utility as utility

def getconnectivity(**kwargs):
    # load relevant kwarg
    if 'filename' in kwargs:
        filename = kwargs.get('filename')
    else:
        filename = ''
        
    # First, list the existing connectivity matrices
    connectivitylist = glob.glob('.\\data\\connectivity\\connectivity_*.csv')
    # find ones that contain the specified filename
    newlist = [cons for cons in connectivitylist if filename in cons]
    
    # If they exist, ask if we want to read them
    if newlist:
        if len(newlist)==1:
            ind = 0
        else:   
            print('Found multiple saved connectivity dataset with the specified filename:')
            utility.print_indexed(connectivitylist)
            ind = int(input('Which one do you want to load?: '))
        
        connectivity = pd.read_csv(newlist[ind])
        _, filename = os.path.split(newlist[ind])
        
    else: # if saved id list doesn't exist, load them from neuprint
        connectivity, filename = getconnectivityfromserver(**kwargs)
        
    return connectivity, filename

def getconnectivityfromserver(**kwargs):
    print('Newly calculating a connectivity matrix!')
    
    # Connect to the neuPrint server
    f = open("authtoken","r")
    tokenstr = f.read()
    c = Client('neuprint.janelia.org', dataset='hemibrain:v1.2.1', token=tokenstr)
    c.fetch_version()
    
    # Get bodyidlist either from the folder or from neuprint
    bodyidlist, filename = getbodyids.getbodyids(**kwargs)
    
    # we create the output dataframe by appending new columns to the bodyidlist
    # which should already be a pandas dataframe
    connectivity = bodyidlist
    
    # Go through all the bodyids and get connections
    for ii in range(len(bodyidlist)):
        if ii%20==0: print('Working on cell #'+str(ii))
        thisId = bodyidlist.bodyId[ii]
        q = """\
            MATCH (a:Neuron)-[w:ConnectsTo]->(b:Neuron)
            WHERE a.bodyId=%s
            RETURN DISTINCT b.bodyId as bodyId, b.type as type, w.weight as w
            """ % thisId
        df = c.fetch_custom(q)
        # Find unique cell types in the output
        types = df.type.unique()
        for thisType in types:
            # Ignore un-labeled downstream neurons
            if thisType is not None:
                # calculate the total weight
                thisCon = np.sum(df.loc[df["type"]==thisType].w)
                # Register this into the output df 
                if thisType not in connectivity.columns:
                    connectivity[thisType] = np.zeros(len(connectivity))
                connectivity.at[ii, thisType] = thisCon
    newfilename = 'connectivity_'+filename
    connectivity.to_csv('./data/connectivity/'+newfilename)
    return connectivity, newfilename
    
    
    