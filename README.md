# Lobula Clustering

This repository contains code to run the analysis presented in the following publication:

[Tanaka, R. & Clark, D. A. (2022) [final title TBD] bioRxiv.](https://www.biorxiv.org/content/10.1101/2022.XXX)

The goal of the analysis is to categorize fragmented lobula visual neurons in the [hemibrain](https://www.janelia.org/project-team/flyem/hemibrain) dataset into putative cell types, based on their connectivity to identified cell types as well as their morphology (layer innervation patterns and the spatial extent of their axon terminals).

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



