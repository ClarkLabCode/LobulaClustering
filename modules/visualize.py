"""

 A bunch of utility visualization functions

"""
## Packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import umap
from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA
## my modules
import modules.utility as utility


# Show quadric surface and scatter
# Expect coeff to be coefficient for quadric model & PCs to be locations of
# synapses in PC spaces
def plotquadricandscatter(pca,coeff,df,**kwargs):
   # apply PCA and calculate rawdepth
   rawdepth, PCs = utility.calcrawdepth(pca, coeff, df)

   # When there are too many dots, we might want to only show a fraction of them
   if 'fraction' in kwargs:
       F = kwargs.get("fraction")
   else:
       F = 0.1
   showind = np.random.rand(len(PCs))<F

   # Find ranges of PC1 and PC2 for mesh thing
   PC1max = np.max(PCs[:,0])
   PC2max = np.max(PCs[:,1])
   PC1min = np.min(PCs[:,0])
   PC2min = np.min(PCs[:,1])

   # create 5 micron dense mesh
   mesh1,mesh2 = np.meshgrid(np.arange(PC1min,PC1max,5),np.arange(PC2min,PC2max,5))
   mesh1 = mesh1.flatten()
   mesh2 = mesh2.flatten()

   # calculate predicted Z (PC3)
   A = np.array([mesh1*0+1, mesh1, mesh2, mesh1**2, mesh2**2, mesh1*mesh2]).T
   mesh3 = np.dot(A,coeff)

   # Visualize!
   fig, ax = plt.subplots(subplot_kw={"projection": '3d'})
   surf = ax.plot_trisurf(mesh1,mesh2,mesh3,cmap ='cool',alpha=0.5)
   sca  = ax.scatter(PCs[showind,0],PCs[showind,1],PCs[showind,2],c=rawdepth[showind],s=1)


# Given a feature matrix (e.g. connectivity) and a label vector (e.g. clusters)
# sort the matrix and visualize, with cluster boundary as dotted line
# assume the specific structure where 1st dimensions are different samples and
# 2nd dimension is different features
def showsortedmatrix(mat,label,**kwargs):

    # sort label
    sortedmat, sortedlabel = utility.sortmatrixbylabel(mat,label)
    # plot
    fig, ax, im = showmatrix(sortedmat.T, **kwargs)
    # add border lines
    mat_height = sortedmat.shape[1]
    for thisLabel in np.unique(sortedlabel):
        lastind = np.max(np.where(sortedlabel==thisLabel))
        plt.plot([lastind-0.5,lastind-0.5],[-0.5,mat_height],'w--')

    # show row labels (optional)
    if 'rowlabel' in kwargs:
        rowlabel = kwargs.get('rowlabel')
        ax.set_yticks(np.arange(len(rowlabel)))
        ax.set_yticklabels(rowlabel)
        ax.tick_params(labelright=True)

    ax.set_ylim(mat_height-0.5,-0.5)
    # colorbar (saturate at 1/99 percentile to ignore outlyers)
    im.set_clim(np.min(mat),np.percentile(mat,99.9))
    fig.colorbar(im,ax=ax)
    return fig, ax


def showmatrix(mat,**kwargs):
    if 'cmapname' in kwargs:
        cmapname = kwargs.get('cmapname')
    else:
        cmapname = 'viridis'

    fig, ax = plt.subplots()
    im = ax.imshow(mat, aspect='auto', cmap=cm.get_cmap(name=cmapname), interpolation='none')
    return fig, ax, im

# Given a feature matrix and a label vector (from clustering)
# Do PCA on the feature matrix and show scatter plots, label color coded
def showsortedPCscatter(mat,label,**kwargs):
    # You can provide the number of components to show through kwargs
    # Otherwise it is defaulted at 5
    if 'n_components' in kwargs:
        n_components = kwargs.get('n_components')
    else:
        n_components = 5

    # do PCA
    pca = PCA(n_components = n_components)
    PCs = pca.fit_transform(mat)
    fig, ax = showsortedscatter(PCs,label,n_show = n_components)
    return fig, ax


def showsortedUMAPscatter(mat,label,**kwargs):
    # You can provide the number of components to show through kwargs
    # Otherwise it is defaulted at 5
    if 'n_components' in kwargs:
        n_components = kwargs.get('n_components')
    else:
        n_components = 5

    # do PCA
    reducer = umap.UMAP(n_components = n_components,random_state=1)
    embedding = reducer.fit_transform(mat)
    fig, ax = showsortedscatter(embedding,label,n_show = n_components)
    return fig, ax



def showUMAPscatter2D(mat,label):
    # adjust dot size according to the number of samples
    dot_size = np.minimum(np.ceil(1000/mat.shape[0]),10)

    n_cat = len(np.unique(label))
    # do PCA
    reducer = umap.UMAP(n_components = 2,random_state=1)
    embedding = reducer.fit_transform(mat)

    cmap = cm.get_cmap('gist_ncar')

    fig, ax = plt.subplots()
    kk = 0
    for ll in np.unique(label):
        ax.scatter(embedding[label==ll,0],embedding[label==ll,1],
                   s=dot_size, c=np.array([cmap(kk/(n_cat+1))]), edgecolors='none',
                   label=str(ll))
        kk += 1
    ax.legend(ncol=4)
    return fig, ax



# Given a feature matrix and a label vector,
# plot labeled data in the raw features space
def showsortedscatter(mat,label,**kwargs):
    # you can either provide the number of columns to show or the indices of
    # columns to show
    if 'col_to_show' in kwargs:
        col_to_show = kwargs.get('col_to_show')
        n_show = len(col_to_show)
    else:
        if 'n_show' in kwargs:
            n_show = kwargs.get('n_show')
        else:
            n_show = 5
        col_to_show = range(n_show)

    # prepare color vector (assume labels to be integers)
    cmap = cm.get_cmap('gist_ncar')

    # get the number of cluster
    n_cluster = len(np.unique(label))

    # do visualization
    fig = plt.figure()
    for ii in range(n_show):
        # index of this column
        coli = col_to_show[ii]
        # diagonal subplots will show histograms of PCs
        ax = fig.add_subplot(n_show,n_show,ii*n_show+ii+1)
        ax.hist(mat[:,coli])
        ax.set_ylabel('Column #'+str(coli))
        # off-diagonal subplots will show scatters
        for jj in range(ii+1, n_show):
            # index of this column
            colj = col_to_show[jj]
            ax = fig.add_subplot(n_show,n_show,ii*n_show+jj+1)
            if ii == 0:
                ax.set_xlabel('Column #'+str(colj))
                ax.xaxis.set_label_position('top')
            for kk in np.unique(label):
                ax.scatter(mat[label==kk,jj],mat[label==kk,ii],s=0.5,c=np.array([cmap(kk/(n_cluster+1))]))

    # create legend (in a separate dedicated plot on the lower left corner)
    ax = fig.add_subplot(n_show,n_show,n_show*(n_show-1)+1)
    for kk in np.unique(label):
        sc = ax.scatter(kk,1,c=np.array([cmap(kk/(n_cluster+1))]))
    ax.set_yticks([])
    ax.set_xticks(range(1,n_cluster+1))
    ax.set_title('Clusters')
    return fig, ax


def plotmeanbycluster(mat,label,**kwargs):
    # read kwargs
    # transpose if axis=1
    if 'axis' in kwargs:
        if kwargs.get('axis'):
            mat = mat.T
    # x axis variable as vector
    if 'x' in kwargs:
        x = kwargs.get('x')
    else:
        x = np.arange(mat.shape[1])

    uniquelabel = np.unique(label)
    n_cluster = len(uniquelabel)
    cmap = cm.get_cmap('gist_ncar')
    fig, ax = plt.subplots()
    for ii in range(n_cluster):
        thislabel = uniquelabel[ii]
        thisN = sum(label==thislabel)
        thismean = np.mean(mat[label==thislabel,:],axis=0)
        thissem = np.std(mat[label==thislabel,:],axis=0)/np.sqrt(thisN)
        thiscol = cmap(ii/(n_cluster+1))
        ax.plot(x,thismean,color=thiscol,label=thislabel)
        ax.fill_between(x,thismean-thissem,thismean+thissem,color=thiscol,alpha=0.2)
    ax.legend(ncol=4)
    return fig, ax

def meanscatterwitherror(mat,label,**kwargs):
    n_feature = mat.shape[1]
    uniquelabel = np.unique(label)
    n_cluster = len(uniquelabel)
    # prepare mean/sem matrix
    mean_mat = np.empty([0,n_feature])
    sem_mat  = np.empty([0,n_feature])
    # calculate mean/sem
    for ll in uniquelabel:
        this_mean = np.reshape(np.mean(mat[label==ll,:],axis=0),[1,n_feature])
        this_sem  = np.reshape(np.std(mat[label==ll,:],axis=0),[1,n_feature])/np.sqrt(np.sum(label==ll))
        mean_mat = np.concatenate((mean_mat,this_mean),axis=0)
        sem_mat  = np.concatenate((sem_mat,this_sem),axis=0)

    fig, ax = plt.subplots(1,n_feature-1)
    cmap = cm.get_cmap('gist_ncar')
    for ii in range(n_feature-1):
        for jj in range(n_cluster):
            thiscol = cmap(jj/(n_cluster+1))
            ax[ii].errorbar(mean_mat[jj,ii],mean_mat[jj,ii+1],xerr=sem_mat[jj,ii],yerr=sem_mat[jj,ii+1],fmt='none',c=thiscol)
            ax[ii].scatter(mean_mat[jj,ii],mean_mat[jj,ii+1],c=np.array([thiscol]))
            ax[ii].text(mean_mat[jj,ii],mean_mat[jj,ii+1],str(uniquelabel[jj]),c=thiscol)
    return fig, ax


def showsortedpiechart(x, **kwargs):
    # expects a 1D numpy array as an input x
    if 'labels' in kwargs:
        labels = list(kwargs.get("labels"))
    else:
        labels = list(np.arange(len(x)))

    if 'cutoff' in kwargs: # ignore entry with a smaller fraction than this
        cutoff = kwargs.get("cutoff")
    else:
        cutoff = 0.1

    if 'ax' in kwargs:
        ax = kwargs.get("ax")
    else:
        fig, ax = plt.subplots(1,1)

    # normalize x
    x = x / np.sum(x)
    sortind = np.argsort(-x)
    sortedx = x[sortind]
    last_over_cutoff = np.max(np.where(sortedx>cutoff))
    x_plot = sortedx[:last_over_cutoff+1]
    x_plot = np.append(x_plot, 1-np.sum(x_plot)) # merge small entries
    newlabel = [labels[i] for i in sortind]
    shortlabel = [newlabel[i] for i in range(last_over_cutoff+1)]
    shortlabel.append('other')
    ax.pie(x_plot,labels=shortlabel)
    return ax
