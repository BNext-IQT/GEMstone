"""A database of genes, represented as bloom filters. Able to load the
addgene-plasmids-sequences data set, encode using bloom filters and search
the data using IOU comparisons.
"""

import json
from collections import defaultdict, namedtuple
import time
import sys, os
import numpy as np
import pickle as p
from joblib import Parallel, delayed
from Bio import SeqIO
from phe import paillier
from p_bloom_filter import encode

data_directory = None
num_cores = 48 # Number of cores for parellel processing
seq_len = 20000

####################
# Search for a query in a "database"
####################
def search(query, data_dir):
    """Searches the database for the 'best match' to the given query. Returns
    relevent information to find the IOU scores of all the genes in the database
    in order to determine the 'best match'

    Args:
        query: The encrypted bloom filter (array) of the gene being searched for.

    Returns:
        A dictionary where each ID maps to the encrypted intersection of the
        gene and query and the magnitude of the gene.
    """
    global num_cores
    global data_directory
    
    data_directory = data_dir
        
    data = os.listdir(data_directory)
    
    print('\nFound %s entries in database\n' % str(len(data)))
    
    data = data[:7000]
    print('Using %s entries from database\n' % str(len(data)))
    
    scores = Parallel(n_jobs=num_cores)(delayed(gen_scores)(id_, query) for id_ in data)
    
    return scores


####################
# Calculate the dot product and magnitude of result based on a sequence ID
####################
def gen_scores(id_, query):
    global data_directory
    global seq_len
    
    seq_file = os.path.join(data_directory, id_)
    
    entry_seq = ''

    with open(seq_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            entry_seq += str(record.seq)
    
    entry_seq = entry_seq[:seq_len]
    
    entry_bloom = encode(entry_seq)
    
    seq_code = entry_seq
    
    return (dotproduct(entry_bloom, query), magnitude(entry_bloom), seq_code)
    
    
####################
# Calculate the dot product between a binary vector and an encrypted vector
####################
def dotproduct(v1, v2):
    """Finds the dot product of two vectors (arrays). The first vector must
    be binary.

    Args:
        v1: A binary vector (array).
        v2: A vector (array).

    Returns:
        The dot product of the two vectors.
    """
    v1_array = np.asarray(v1)
    v2_array = np.asarray(v2)
    #dot = sum([v2[i] for i,_ in enumerate(v1) if v1[i] == 1])
    dot = np.dot(v1_array, v2_array)
    
    return dot


####################
# Calculate the magnitude of a binary vector
####################
def magnitude(v):
    """Finds the magnitude of a binary vector (array).

    Args:
        v: A binary vector (array).

    Returns:
        The magnitude of the vector.
    """
    
    return sum(v)                