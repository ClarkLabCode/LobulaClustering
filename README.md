# Lobula Clustering

This repository contains code to run the analysis presented in the following publication:

[Tanaka, R. & Clark, D. A. (2022) [final title TBD] bioRxiv.](https://www.biorxiv.org/content/10.1101/2022.XXX)

The goal of the analysis is to categorize fragmented lobula visual neurons in the [hemibrain](https://www.janelia.org/project-team/flyem/hemibrain) dataset into putative cell types, based on their connectivity to identified cell types as well as their morphology (layer innervation patterns and the spatial extent of their axon terminals).










## Geting started with hemibrain

You can access hemibrain dataset on your browser at [neuPrint](https://neuprint.janelia.org/) website. A google account is required to sign up.

The code in this repository is based on [neuprint-python](https://github.com/connectome-neuprint/neuprint-python) package provided by Janelia. I would recommend going through their [quickstart guide](https://connectome-neuprint.github.io/neuprint-python/docs/quickstart.html) first.

To fetch data from the hemibrain database, **you need authentification token** attached to your google acount, which can obtained from the top-right menu of the neuPrint website (see [quickstart guide](https://connectome-neuprint.github.io/neuprint-python/docs/quickstart.html)). You then need to put your token in a file named **"authtoken"** and put that under the top directory of the repo.

## Writing Cypher query

To fetch data from hemibrain dataset, we pass query to the server written in Cypher Query Language.
See the [manual](https://neuprint.janelia.org/public/neuprintuserguide.pdf) for how to write a query.
