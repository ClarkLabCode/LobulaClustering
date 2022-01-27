"""

 Get bodyId of lobula neurons of interest by

 - reading csv if it has already been saved
 - running neuprint query if not

"""

## Packages
from neuprint import Client
import pandas as pd
import numpy as np
import os
import glob
import modules.utility as utility

def getbodyids(**kwargs):
     # load relevant kwarg
     # provide some non-existent filename when downloading new bodyid
    if 'filename' in kwargs:
        filename = kwargs.get('filename')
    else:
        filename = ''

    # just making explicit what is being called...
    print('Running getbodyids...')

    # Check if we already have fetched bodyids from hemibrain
    bodyidlistlist = glob.glob('.\\data\\bodyidlist\\*.csv')
    newlist = [listname for listname in bodyidlistlist if filename in listname]

    # If they exist, ask if we want to read them
    if newlist:
         if len(newlist)==1:
             ind = 0
         else:
             print('Found multiple lists of bodyIds with the provided filename:')
             utility.print_indexed(newlist)
             ind = int(input('Which one do you want to load?: '))
         # load the list
         bodyidlist = pd.read_csv(newlist[ind])
         _, filename = os.path.split(newlist[ind])
    else:
        # if saved id list doesn't exist, load them from neuprint
        # We only care about un-labeled neurons
        bodyidlist, filename = getbodyidsfromservertypenull(**kwargs)
    return bodyidlist, filename


def getbodyidsfromservertypenull(**kwargs):
    # check if upper/lower bounds of synapse numbers are provided as kwargs
    if 'ub' in kwargs:
        ub = kwargs.get('ub')
    if 'lb' in kwargs:
        lb = kwargs.get('lb')

    # just making explicit what is being called...
    print('Running getbodyidsfromservertypenull...')

    # Run neuPrint query to get bodyIds of lobula small terminals
    # which we will analyze by running clustering
    print('Will fetch bodyids from neuPrint')
    # Connect to the server
    f = open("authtoken","r")
    tokenstr = f.read()
    c = Client('neuprint.janelia.org', dataset='hemibrain:v1.2.1', token=tokenstr)
    c.fetch_version()

    # ask upper/lower bounds of the synapse counts (if not provided)
    if not 'ub' in locals():
        ub = input('Enter upper bound of total synapse count: ')
    if not 'lb' in locals():
        lb = input('Enter lower bound of total synapse count: ')

    # Define query
    q = """\
        MATCH (a:Neuron)
        WHERE a.type IS NULL AND apoc.convert.fromJsonMap(a.roiInfo)["LO(R)"].pre=a.pre AND apoc.convert.fromJsonMap(a.roiInfo)["LO(R)"].post=a.post AND a.pre+a.post<%s AND a.pre+a.post>%s
        RETURN DISTINCT a.bodyId as bodyId
        """ % (ub,lb)

    bodyidlist = c.fetch_custom(q)
    print('Found',len(bodyidlist),'cells without labels. Saving...')
    filename = 'typenullbodyidlist_ub'+ub+'_lb'+lb+'.csv'
    bodyidlist.to_csv('./data/bodyidlist/'+filename);
    return bodyidlist, filename

# This is only for validation and will be called directly from scripts
def getbodyids_singletype(**kwargs):
    if 'celltype' in kwargs:
        celltype = kwargs.get('celltype')
    else:
        celltype = input('Enter the cell type of interest: ')

    # just making explicit what is being called...
    print('Running getbodyids_singletype...')

    # we need this to avoid cells with no synapse of requested type
    if 'synapseType' in kwargs:
        synapseType = kwargs.get('synapseType')
    else:
        synapseType = 'post'

    # Fetch and save body Ids of a set cell type
    # for validation of depth extraction
    print('Fetching bodyIds fron the server...')

    # Connect to the server
    f = open("authtoken","r")
    tokenstr = f.read()
    c = Client('neuprint.janelia.org', dataset='hemibrain:v1.2.1', token=tokenstr)
    c.fetch_version()

    # Define query
    q = """\
        MATCH (a:Neuron)-[:Contains]->(:SynapseSet)-[:Contains]->(s:Synapse)
        WHERE a.type = '%s' AND s.`LO(R)` AND s.type = '%s'
        RETURN DISTINCT a.bodyId as bodyId
        """ % (celltype,synapseType)

    bodyidlist = c.fetch_custom(q)
    print('Found',len(bodyidlist),celltype,'. Saving...')
    filename = 'bodyidlist_'+celltype+'.csv'
    bodyidlist.to_csv('./data/bodyidlist/'+filename);
    return bodyidlist, filename
