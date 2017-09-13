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

from phe import paillier
from p_bloom_filter import encode

# Type of sequence used in database
SEQUENCE_TYPE = "public_addgene_full_sequences"

# Database of genes
data = None

# Gene data structure holding the name, principle investigator, a partial
# sequence and the corresponding bloom filter of a gene.
Gene = namedtuple("Gene", "name, pi, sequence, bloom")

def search(query):
    """Searches the database for the 'best match' to the given query. Returns
    relevent information to find the IOU scores of all the genes in the database
    in order to determine the 'best match'

    Args:
        query: The encrypted bloom filter (array) of the gene being searched for.

    Returns:
        A dictionary where each ID maps to the encrypted intersection of the
        gene and query and the magnitude of the gene.
    """
    global data
    global Gene

    # Encode data
    if data:
        print("data already encoded")
    else:
        print("endcoding data...")
        #print("Loading data...")
        encode_data()

        #data = p.load(open('../encoded_addgene.p', 'rb'))
        print("...database complete")
        #print("...database loaded")

    scores = {}

    for id_ in data:
        d = data[id_].bloom
        scores[id_]=(dotproduct(d, query), magnitude(d))

    return scores

def dotproduct(v1, v2):
    """Finds the dot product of two vectors (arrays). The first vector must
    be binary.

    Args:
        v1: A binary vector (array).
        v2: A vector (array).

    Returns:
        The dot product of the two vectors.
    """

    dot = 0
    for i in range(0, len(v1)):
        if v1[i] == 1:
            dot += v2[i]

    return dot


def magnitude(v):
    """Finds the magnitude of a binary vector (array).

    Args:
        v: A binary vector (array).

    Returns:
        The magnitude of the vector.
    """
    sum = 0
    for x in v:
        sum += x
    return sum

def get_gene(id_):
    """Returns gene with given ID.

    Args:
        id_: ID corresponding to desired gene.

    Returns:
        The corresponding gene.
    """
    if data:
        return data[id_]
    return None

def encode_data():
    """Reads in the addgene-plasmids-sequences data from a json file and stores
    the important information in an internal data structure. Also genenerates
    and stores the bloom filter for each gene in the database.
    """
    global data
    # Reads in data from JSON file
    try:
	   gene_data = os.path.join('../data', 'addgene-plasmids-sequences.json')
	   with open(gene_data) as data_file:
           raw_data = json.load(data_file)
       data_file.close()

    except FileNotFoundError:
        print("\nFileNotFoundError: [Errno 2] No such file or directory: ",
            "'addgene-plasmids-sequences.json'\n")
        print("File should be in current folder. Encode failed.")
        sys.exit(2)

    plasmids = raw_data['plasmids']

    # Initialization of internal data structures for set.
    data = defaultdict(Gene)
    id_ = 0

    # Read in and organize important fields of data from the raw data.
    for k in plasmids:
        name = k['name']
        if k['pi']:
            pi = k['pi'][0]
            if k['sequences'][SEQUENCE_TYPE]:
                sequence = k['sequences'][SEQUENCE_TYPE][0]
                bf = encode(sequence)
                gene = Gene(name = name, pi = pi, sequence = sequence, bloom = bf)
                data[id_] = gene
                id_ += 1
