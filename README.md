# Lobula Clustering

This repository contains code to run the analysis presented in the following publication:

[Tanaka, R. & Clark, D. A. (2022) [final title TBD] bioRxiv.](https://www.biorxiv.org/content/10.1101/2022.XXX)

The goal of the analysis is to categorize fragmented lobula visual neurons in the [hemibrain](https://www.janelia.org/project-team/flyem/hemibrain) dataset into putative cell types, based on their connectivity to identified cell types as well as their morphology (layer innervation patterns and the spatial extent of their axon terminals).

Note that the results of the analysis presented in the paper (i. e., the IDs and cluster affiliations of the cells analyzed) are provided as a supplementary materials of the paper and you do not need to rerun the whole analysis if you are only interested in taking a closer look at the results.

## Geting started

The code in this repository is written in Python 3. We recommend setting up a virtual environment for this repo using [conda](https://docs.conda.io/en/latest/).
The following packages are required to run the code:

- [neuprint-python](https://github.com/connectome-neuprint/neuprint-python)
- [pandas](https://pandas.pydata.org/)
- [numpy](https://numpy.org/)
- [matplotlib](https://matplotlib.org/)
- [scikit-learn](https://scikit-learn.org/stable/index.html)
- [scipy](https://scipy.org/)
- [umap-learn](https://umap-learn.readthedocs.io/en/latest/)

Install these packages in your virtual environment by running ```conda install -c <channel name> <package name>```, if you are using conda. **neuprint-python** is distributed from **flyem-forge**, and **umap-learn** from **conda-forge**. Everything else should be available from the default channel.

Accesing the hemibrain dataset requires a user authentification token. Visit the [neuPrint+](https://neuprint.janelia.org/) website and sign up (a Google account required). From the **Account** menu, find your Auth Token and copy it into a plain text file named **authtoken** (without extension), place it directly under the root folder of the cloned repository.

## Running the main clustering script

Run ```lobulaclustering.py``` to perform the clustering analysis of the fragmented, unlabeled lobula neurons, as presented in **Figs. 3, 5** in the paper. Please refer to the **Method** section of the [paper](https://www.biorxiv.org/content/10.1101/2022.XXX) for the detailed description of the analysis.

The script will
- load the connectivity matrix
- load the morphology matrix
- perform hiearchical clustering
- visualize the results of the clustering
- perform additional analysis of connectivity from the resultant clusters to LC/LPLC neurons

The results of the clustering will be saved under **data/results**.

Free parameters of this analysis are as follows (see the paper for the details), with the values used in the paper in square brackets:
1. Upper and lower boundaries of the numbers of (pre)synapses of the neurons to be analyzed [500, 50]
2. The cell type to be used to as a landmark to define the layers of lobula [LT1]
3. Upper and lower limits of synaptic positions as well as the bin width to create a histogram of synapse distribution along the layer depth (in microns) [-20, 50, 5]
4. Relative weighting between connectivity, depth, and spread features [5, 3, 1]
5. Number of clusters [40]

When run for the first time, the script will ask you to type in parameters 1 through 3 (through the command line), whereas 4 and 5 are hard-coded (lines 44/45).

When run for the first time, the script downloads the bodyIds of neurons that matches the synapse count criteria, as well as their synapse coordinates and connectivity. It also downloads synapse coordinates of the landmark cell. These process can take long, especially when you are analyzing a large number of neurons. We recommend you to initially set the range of synapse count small (e. g., between 110 and 100), so you can check if the code runs through properly without waiting too long. The list of bodyIds, connectivity, synapse coordinates, and morphological features (i. e., innervation depth and synapse spread) are all saved in the data directory, such that you do not need to repeat the time-consuming process of data download in the subsequent runs.

In subsequent runs, the script will look for the saved connectivity and morphology matrices under **data/connectivity**,  **data/depth**, and **data/spread**. If it finds saved matrices, it will ask you to specify which matrix you want to use. You can also opt to run the analysis using a different set of parameters.


## Running the validation

Run ```morphology_validation.py``` to generate the results shown in **Fig. 2** of the paper.
This script will download the coordinates of postsynapses of specified LC and LPLC neuron types, and calculates morphological summary features.


## Organizations of directory

- ```lobulaclustering.py``` : the main clustering description
- ```morphology_validation.py``` : the script to validate the morphology summary features by analyzing LC/LPLCs with known morphology
- modules : a folder containing modules to download, save, preprocess, and load connectivity and morphology data
- data : a folder containing the end results of the clustering as well as downloaded and preprocessed intermediate data
