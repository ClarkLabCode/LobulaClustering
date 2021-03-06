
"""

 Main clustering script

 This function load saved connectivity and morphology (depth/spraed) matrices
 and run hiearchical clustering on them

 If this is the first time you run this, subroutines will
 - download bodyIds of lobula intrinsic neurons (given range of synapse counts)
 - fetch connectivity to downstream neuron types (given a list of bodyIds)
 - Download synapse coordinates of cells (given bodyIds)
 - Do PCA and fit a quadric model to some landmark cells
 - Rotate each synapse into the PC space and calculate their relative depth
 - Calculate spread and depth histogram based on PC-rotated synapse coordinates and relative depth

"""
# Packages
from neuprint import Client
import pandas as pd
import numpy as np
import os
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA
from sklearn.cluster import AgglomerativeClustering
import scipy.cluster.hierarchy as sch
from scipy import signal

# Our own modules
import modules.getbodyids as getbodyids
import modules.getsynapses as getsynapses
import modules.getmorphology as getmorphology
import modules.getconnectivity as getconnectivity
import modules.visualize as visualize
import modules.utility as utility

# just making explicit what is being called...
print('Running lobulaclustering.py...')

## 0. Hard-coded analysis parameters
data_weight = (5,3,1) # how much we trust each dataset (con/dep/spr)
n_cluster = 40

nShow = 30 # this determines how many connectivity features we want to see in the plot (does not affect the analysis itself)

## 1. data preparation

# note: make sure this runs when running the script for the first time

# load connectivity matrix
# or create one if there is none saved
connectivity, con_fn = getconnectivity.getconnectivity()

# load morphology matrix
# or create one if there is none saved
depth, spread, dep_fn = getmorphology.getmorphology()

# Check connectivity and morphology are based on the same bodyidlist
# The assumption is that the order of the bodyId should be the same across these
# three files. This should be true by constructrion (they are created by appending
# new columns to bodyidlist)
if con_fn[con_fn.find('bodyidlist'):] != dep_fn[dep_fn.find('bodyidlist'):]:
    print('Connectivity and morphology matrices are based on different sets of cells. Aborting')
else:
    # Show what hard-coded parameters we are using + which dataset we are using
    print('Connectivity matrix we are using: ', con_fn)
    print('Morphology matrix we are using: ', dep_fn)
    print('Relative weight between connectivity, depth, spread: ',data_weight)
    print('#Cluster requested: ',n_cluster)

    ## Data preparation
    # Find which is the first data column (this should be in principle always 3)
    con_datastart = connectivity.columns.get_loc("bodyId")+1
    dep_datastart = depth.columns.get_loc("bodyId")+1
    spr_datastart = spread.columns.get_loc("bodyId")+1

    # Convert into numpy 2d array
    mat_con = connectivity.iloc[:,con_datastart:].to_numpy()
    mat_dep = depth.iloc[:,dep_datastart:].to_numpy()
    mat_spr = spread.iloc[:,spr_datastart:].to_numpy()

    # normalize by total dispersion (and show dispersion)
    disp_con = np.sum(np.var(mat_con,axis=0))
    disp_dep = np.sum(np.var(mat_dep,axis=0))
    disp_spr = np.sum(np.var(mat_spr,axis=0))
    print('total dispersion of connectivity ',disp_con)
    print('total dispersion of depth ',disp_dep)
    print('total dispersion of spread ',disp_spr)

    # normalize and weight
    norm_mat_con = mat_con / disp_con * data_weight[0]
    norm_mat_dep = mat_dep / disp_dep * data_weight[1]
    norm_mat_spr = mat_spr / disp_spr * data_weight[2]

    # also, get labels for columns (just in case)
    label_con = connectivity.columns[con_datastart:]
    label_dep = depth.columns[dep_datastart:]
    label_str = spread.columns[spr_datastart:]

    # re-order connectivity matrix and its labels by total number of connectivity
    # because we don't care about rare ones (for visualization)
    total_connection = np.sum(mat_con, axis=0)
    important_target_ind = np.argsort(-total_connection)[:nShow]

    # Concatenate the matrices
    mat_all = np.concatenate((norm_mat_con,norm_mat_dep,norm_mat_spr),axis=1)

    ## Actual Clustering
    # Do clustering with ward minimization & show the dendrogram
    fig, ax = plt.subplots()
    linkage = sch.linkage(mat_all, method='ward', metric='euclidean')
    dendrogram = sch.dendrogram(linkage, truncate_mode='lastp', p =n_cluster)
    clabel = sch.fcluster(linkage, n_cluster, criterion='maxclust')
    ax.set_title('Dendrogram of cells of interest (Fig. 3A)')

    ## Visualization and post-processing

    # sort and visualize
    # visualize the connectivity matrix
    fig, ax = visualize.showsortedmatrix(mat_con[:,important_target_ind],clabel,rowlabel=label_con[important_target_ind])
    ax.set_title('Connectivity (Fig. 3B)')
    ax.set_xlabel('Cells of interest')

    fig, ax = visualize.showsortedmatrix(mat_dep,clabel)
    ax.set_title('Innervation Depth (Fig. 3C)')
    ax.set_xlabel('Cells of interest')
    ax.set_ylabel('#Depth bin')

    fig, ax = visualize.showsortedmatrix(mat_spr,clabel)
    ax.set_title('Synapse Spread (um) (Fig. 3D)')
    ax.set_xlabel('Cells of interest')
    ax.set_ylabel('PC axis')

    # visualize the clusters in the PC space
    #visualize.showsortedUMAPscatter(mat_all,clabel,n_components=2)
    fig, ax = visualize.showUMAPscatter2D(mat_all,clabel)
    ax.set_title('UMAP on the concatenated weighted feature matrix (Fig. 3E)')
    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')
    utility.reporttargetpercluster(mat_con, label_con, clabel)

    # visualize mean depth profile for each cluster
    fig, ax = visualize.plotmeanbycluster(mat_dep, clabel)
    ax.set_title('Mean synapse per depth bin for each cluster (Appendix Figs)')
    ax.set_xlabel('#Depth bin')
    ax.set_ylabel('#synapses')

    # visualize mean spread profile for each cluster
    fig, ax = visualize.meanscatterwitherror(mat_spr,clabel)
    ax[0].set_title('synapse spread (Appendix Figs)')
    ax[0].set_xlabel('spread along PC1 (um)')
    ax[0].set_ylabel('spread along PC2 (um)')
    ax[1].set_xlabel('spread along PC2 (um)')
    ax[1].set_ylabel('spread along PC3 (um)')

    plt.show()

    # print morphology parameters
    for cc in np.unique(clabel):
        print('Mean spread of cluster#',cc,np.mean(mat_spr[clabel==cc,:],axis=0))

    # Save results (uncomment for actually saving)
    outdf = depth.bodyId.to_frame()
    outdf.insert(1,"cluster",clabel)
    outfn = 'cluster_N'+str(n_cluster)+dep_fn[5:]
    outdf.to_csv('./data/result/'+outfn)


    ## Additional analysis ##
    # Connectivity from the clusters to LCs
    # limiting this to "classical LCs" up to Wu Nern 2016
    # interested readers can add LC beyond 26
    # list of LC neurons we are going to analyze
    print('Running the additional LC connectivity analysis...')
    LC_list = ('LC4','LC6','LC9','LC10','LC11','LC12','LC13','LC14','LC15',
               'LC16','LC17','LC18','LC20','LC21','LC22','LC24','LC25','LC26',
               'LPLC1','LPLC2','LPLC4')
    LC_list = pd.Index(LC_list)
    LC_index = []
    for LC in LC_list:
        LC_index.append(list(label_con).index(LC))
    # pull out connectivity from cells of interest to LCs
    con_LC = mat_con[:,LC_index]
    # sum within each cluster
    con_LC_byCluster = np.empty([n_cluster, len(LC_list)])
    for cc in np.unique(clabel):
        con_LC_byCluster[cc-1,:] = np.sum(con_LC[clabel==cc,:],axis=0)
    # normalize for each cell type
    norm_con_LC_byCluster = con_LC_byCluster / np.sum(con_LC_byCluster,axis=0)

    # Visualize as a pie chart
    fig, ax = plt.subplots(3,7)
    for i in range(3):
        for j in range(7):
            ax[i,j].set_title(LC_list[i*7+j])
            visualize.showsortedpiechart(norm_con_LC_byCluster[:,i*7+j],
                                         cutoff=0.05,
                                         labels=np.unique(clabel),
                                         ax=ax[i,j])
    fig.suptitle('LP/LPLC inputs by cell of interest clusters (Fig. 5A)')


    # Visualize as a dendrogram
    linkage_reverse = sch.linkage(norm_con_LC_byCluster.T, method='ward', metric='euclidean')
    fig, ax = plt.subplots()
    dendrogram_reverse = sch.dendrogram(linkage_reverse, labels=LC_list)
    out_ind = dendrogram_reverse['leaves']
    ax.set_title('Clustering of LC/LPLCs by their connectivity to cell of interest clusters (Fig. 5B)')

    # visualization
    fig, ax = plt.subplots()
    im = ax.imshow(norm_con_LC_byCluster[:,out_ind].T)
    ax.set_xticks(np.arange(n_cluster))
    ax.set_xticklabels(np.arange(n_cluster)+1)
    ax.set_xlabel('cluster')
    ax.set_yticks(np.arange(len(LC_list)))
    ax.set_yticklabels(LC_list[out_ind])
    ax.set_title('Normalized mean LP/LPLC inputs by cell of interest clusters (Fig. 5B)')
    fig.colorbar(im,ax=ax)
    plt.show()
