"""

 Given a list of bodyId of lobula neurons of interest, download all the synapses
 and save them so we don't need to redo this everytime we want to re-calculate
 morphological metrics

 Because synapses belong to each cell and does not depend on which cell you are
 analyzing (e.g. different synapse count UB/LBs) let's have one big folder that
 stores all the synapse csv

"""

## Packages
from neuprint import Client
import pandas as pd
import numpy as np
import os
import glob


# This will go through the "synapselist" folder and download synapse if necessary
def getsynapses(bodyid,synapseType):

    # refer to different directory depending on which synapse type you are using
    if synapseType=='pre':
        synapseDir = '.\\data\\synapselist\\'
    else:
        synapseDir = '.\\data\\postsynapselist\\'

    # just making explicit what is being called...
    print('Running getsynapses...')

    # First, check if synapses of this neuron has been already saved
    thisSynapseFile = glob.glob(synapseDir+str(bodyid)+'.csv')

    # if it does not exist, download
    if not thisSynapseFile:
        print('Downloading the '+synapseType+'synapses of cell#'+str(bodyid))
        # First, connect to the neuPrint server
        f = open("authtoken","r")
        tokenstr = f.read()
        c = Client('neuprint.janelia.org', dataset='hemibrain:v1.2.1', token=tokenstr)
        c.fetch_version()
        # Prepare query
        q = """\
            MATCH (a:Neuron)-[:Contains]->(:SynapseSet)-[:Contains]->(s:Synapse)
            WHERE a.bodyId=%s AND s.type = '%s' AND s.`LO(R)`
            RETURN DISTINCT s.location.x as x, s.location.y as y, s.location.z as z
            """ % (bodyid,synapseType)
        df = c.fetch_custom(q)
        # save it
        df.to_csv(synapseDir+str(bodyid)+'.csv')
    else:
        # load it otherwise
        df = pd.read_csv(synapseDir+str(bodyid)+'.csv')
    return df
