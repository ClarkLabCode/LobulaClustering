"""

 Small utility functions
 
"""
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

# Given pca fit to a landmark, quadric model, and a dataframe with (native) 
# x/y/z synapse locations in it, calculate (for each synapse) lobula layer depth 
# as deviation from predicted PC3 position 
def calcrawdepth(pca,coeff,df):
    # First, flatten the input and convert its unit to microns
    X = df["x"].to_numpy().flatten()*8/1000
    Y = df["y"].to_numpy().flatten()*8/1000
    Z = df["z"].to_numpy().flatten()*8/1000
    XYZ = np.array([X,Y,Z]).T
    
    # apply PCA
    PCs = pca.transform(XYZ)
    
    # flatten PCs again
    PC1 = PCs[:,0].flatten()
    PC2 = PCs[:,1].flatten()
    PC3 = PCs[:,2].flatten()

    # construct predictors
    A = np.array([PC1*0+1, PC1, PC2, PC1**2, PC2**2, PC1*PC2]).T
    
    # apply the model to predict PC3
    predPC3 = np.dot(A,coeff)
    rawdepth = PC3 - predPC3
    
    return rawdepth, PCs

def sortmatrixbylabel(mat,label):
    # Sort along dimension 0 (sorting rows)
    # prepare empty 2d array
    sortedmat = np.empty((0,mat.shape[1]))
    sortedlabel= []
    # go through all the labels and sort
    for thisLabel in np.unique(label):
        thismat = mat[label==thisLabel,:]
        sortedmat = np.concatenate((sortedmat,thismat),axis=0)
        sortedlabel = np.append(sortedlabel,label[label==thisLabel])
    return sortedmat, sortedlabel        

# given a feature matrix, label for features, and cluster vector, report the
# name of the most prominent N features
def reporttargetpercluster(mat,col_label,cluster,**kwargs):
    if 'n_show' in kwargs:
        n_show = kwargs.get('n_show')
    else:
        n_show = 5
        
    for this_cluster in np.unique(cluster):
        this_sum = np.mean(mat[cluster==this_cluster,:], axis=0)
        important_target_ind = np.argsort(-this_sum)[:n_show]
        
        # show results
        print('Cluster #)'+str(this_cluster))
        for ii in range(len(important_target_ind)):
            this_ind = important_target_ind[ii]
            print('- target #'+str(ii)+' '+col_label[this_ind]+' mean synapse count: '+str(round(this_sum[this_ind],2)))
            

def print_indexed(in_list):
    for ii in range(len(in_list)):
        print(ii,' : ',in_list[ii])
        
    