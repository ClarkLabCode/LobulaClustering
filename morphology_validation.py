"""

Validate and optimize morphology metrics

Tangential based estimation of synapse depth assumes that lobula layers are fit
well enough with a quadric surface and layers are parallel to each other, which
can introduce systematic error correlated with XY position and can cause over
dividing in the clustering.
This script visualizes postsynapses of LC neurons with known layer innervation
patterns and see if there is anything that can be done on that

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

# cell types to analyze
ctlist = ('LC4','LC6','LC9','LC11','LC12','LC13','LC15','LC16','LC17','LC18',
          'LC20','LC21','LC22','LC24','LC25','LC26','LPLC1','LPLC2')
#ctlist = ('LC4','LPLC2')

# load morphology matrices (calculate if not existing)
label = [] # list to store cell type labels
all_dep = np.empty([0,14])
all_spr = np.empty([0,3])
mean_dep = np.empty([0,14]) # median depth histogram
mean_spr = np.empty([0,3])
sem_spr = np.empty([0,3])

for ct in ctlist:
    # this is necessary only for the first time running the script
    getbodyids.getbodyids_singletype(celltype=ct,synapseType='post')

    # load morphology matrix
    depth, spread, dep_fn = getmorphology.getmorphology(filename=ct,synapseType='post',
                                    landmarkname='LT1',minD=-20,maxD=50,binSize=5)
    # append label
    for i in range(len(depth)):
        label.append(ct)

    dep_datastart = depth.columns.get_loc("bodyId")+1
    spr_datastart = spread.columns.get_loc("bodyId")+1
    
    mat_dep = depth.iloc[:,dep_datastart:].to_numpy()
    mat_spr = spread.iloc[:,spr_datastart:].to_numpy()
    
    all_dep = np.concatenate((all_dep,mat_dep),axis=0)
    all_spr = np.concatenate((all_spr,mat_spr),axis=0)
    
    # calculate normalized mean innervation depth
    this_mean_dep = np.reshape(np.mean(mat_dep,axis=0),[1,14])
    this_mean_dep = this_mean_dep / np.sum(this_mean_dep)    
    mean_dep = np.concatenate((mean_dep,this_mean_dep),axis=0)
    
    # calculate mean synapse spread
    this_mean_spr = np.reshape(np.mean(mat_spr,axis=0),[1,3])
    mean_spr = np.concatenate((mean_spr,this_mean_spr),axis=0)
    this_sem_spr = np.reshape(np.std(mat_spr,axis=0)/np.sqrt(mat_spr.shape[0]),[1,3])
    sem_spr = np.concatenate((sem_spr,this_sem_spr),axis=0)
    

# cast label to np array
label = np.asarray(label)
# need this because some visualization function asks for integer labels
intlabel = [jj for jj in range(len(ctlist)) for ii in range(np.sum(label==ctlist[jj]))]
intlabel = np.asarray(intlabel)
### visualization

## 1. Show mean normalized innervation depth per cell type
fig, ax, im = visualize.showmatrix(mean_dep.T,cmapname='GnBu')
# labels and stuff
ax.set_yticks(np.arange(15)-0.5)
ax.set_yticklabels(np.arange(-20,55,5))
ax.set_xticks(np.arange(len(ctlist)))
ax.set_xticklabels(ctlist,rotation=45,ha='right')
fig.colorbar(im,ax=ax)

## 2. Show cell-averaged spread feature against each other
# with SEM 
fig, ax = plt.subplots(1,2)
for ii in range(2):
    ax[ii].errorbar(mean_spr[:,ii],mean_spr[:,ii+1],xerr=sem_spr[:,ii],yerr=sem_spr[:,ii+1],fmt='none')
    ax[ii].scatter(mean_spr[:,ii],mean_spr[:,ii+1])
    for jj in range(len(ctlist)):
        ax[ii].text(mean_spr[jj,ii],mean_spr[jj,ii+1],ctlist[jj])
        ax[ii].set_xlabel('PC#'+str(ii+1)+' um')
        ax[ii].set_ylabel('PC#'+str(ii+2)+' um')
        

## individual cell data (variable because of clipping etc)
# visualize.showsortedPCscatter(all_dep, intlabel)
# visualize.showsortedscatter(all_spr, intlabel,n_show=3)
