# Lobula Clustering

Extend lobula clustering we have performed as a part of LPLC1 to all neurons in lobula.
Re-creating the code in a seprate folder so it is easy to eventually publish this as a separate repo.


## Geting started with hemibrain

You can access hemibrain dataset on your browser at [neuPrint](https://neuprint.janelia.org/) website. A google account is required to sign up.

The code in this repository is based on [neuprint-python](https://github.com/connectome-neuprint/neuprint-python) package provided by Janelia. I would recommend going through their [quickstart guide](https://connectome-neuprint.github.io/neuprint-python/docs/quickstart.html) first.

To fetch data from the hemibrain database, **you need authentification token** attached to your google acount, which can obtained from the top-right menu of the neuPrint website (see [quickstart guide](https://connectome-neuprint.github.io/neuprint-python/docs/quickstart.html)). You then need to put your token in a file named **"authtoken"** and put that under the top directory of the repo.

## Writing Cypher query

To fetch data from hemibrain dataset, we pass query to the server written in Cypher Query Language.
See the [manual](https://neuprint.janelia.org/public/neuprintuserguide.pdf) for how to write a query.
