"""A database of genes, represented as bloom filters. Able to load the
addgene-plasmids-sequences data set, encode using bloom filters and search
the data using IOU comparisons.
"""

import json
from collections import defaultdict, namedtuple
import time
import sys

from bloom_filter import encode

# Type of sequence used in database
SEQUENCE_TYPE = "public_addgene_full_sequences"

# Database of genes
data = None

# Gene data structure holding the name, principle investigator, a partial
# sequence and the corresponding bloom filter of a gene.
Gene = namedtuple("Gene", "name, pi, sequence, bloom")

def search(query):
    """Searches the database for the 'best match' to the given query. Uses
    intersection-over-union (IOU) comparison of the bloom filters to determine
    closest match.

    Args:
        query: The bloom filter (array) of the gene being searched for.

    Returns:
        The 'Gene' that is the 'best match' to the query.

        The IOU for the 'best match' and the query.
    """
    global data

    # Encode data
    if data:
        print("data already encoded")
    else:
        print("endcoding data...")
        encode_data()
        print("...database complete")

    max_iou = 0
    best_id = 0

    query_mag = magnitude(query)

    for id_ in data:
        score = iou(data[id_].bloom, query, query_mag)
        if score > max_iou:
            max_iou = score
            best_id = id_

    return data[best_id], max_iou

def iou(data, query, query_mag):
    """Finds the IOU for two bloom filters.

    Args:
        data: The bloom filter (array) of a particular gene from the database.
        query: The bloom filter (array) of the gene being searched for.
        query_mag: The magnitude of the gene being searched for. (Since this
            remains constant in all IOU calcutlations in the search)

    Returns:
        The IOU for the two genes.
    """
    intersection = dotproduct(data, query)
    union = magnitude(data) + query_mag - intersection
    return intersection/union

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

def encode_data():
    """Reads in the addgene-plasmids-sequences data from a json file and stores
    the important information in an internal data structure. Also genenerates
    and stores the bloom filter for each gene in the database.
    """
    global data
    # Reads in data from JSON file
    try:
        with open("../addgene-plasmids-sequences.json") as data_file:
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
            if k['sequences']["public_addgene_full_sequences"]:
                sequence = k['sequences'][SEQUENCE_TYPE][0]
                bf = encode(sequence)
                gene = Gene(name = name, pi = pi, sequence = sequence, bloom = bf)
                data[id_] = gene
                id_ += 1
