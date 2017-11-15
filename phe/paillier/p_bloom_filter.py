"""Bloom filter for genetic sequences. Stored as unsigned character array."""

import sys
from array import array
from collections import defaultdict

#BLOOM FILTER DEFAULTS
# Information calculated from https://krisives.github.io/bloom-calculator/
K = 8
H = hash
#HASH_MAX = sys.maxsize + 1
SIZE = 500

def encode(gene, size=SIZE, k=K, h=H, HASH_MAX=sys.maxsize + 1):
    """Creates a bloom filter. Used to encode a genetic sequence.

    Args:
        gene: A string holding all or part of a DNA sequence.
        size: The size of the bloom filter. Set to the default size if
            no size is given.
        k: The size of the k-mer in the filter (how many nucleotides
            [characters] are encoded at once). Set to the default k if no
            k is given.
        h: The hash used to encode each k-mer entered in the bloom filter.
            Set to the default hash if no hash is given.

    Returns:
        The corresponding bloom filter. An array where each each entry a hashed
        k-mer maps to is a one and all other entries are zero.

        The number of of unique k-mers in the gene.
    """
    bf = initialize_bloom_filter(size)
    gene = gene.upper()                         # Make gene all uppercase.

    # Loop through all k-mers for gene.
    for n in range(0, len(gene)-k + 1):
        # Get k-mer of length k and hash it.
        k_mer = gene[n:n + k]                   # TODO: ignore case, 'N's?
        k_hash = (h(k_mer) + HASH_MAX) % size   # Make hash positive and within
                                                # size of filter.

        # Set entry in the bloom filter corresponding to the hashed k-mer to one.
        bf[k_hash] = 1

    return bf


def initialize_bloom_filter(size=SIZE):
    """ Creates empty bloom filter.

    Args:
        size: The size of the bloom filter. Set to the default size if
            no size is given.

    Returns:
        An unsigned character array of the given size where each entry is 0.
        The empty bloom filter.
    """
    bf = array('b',[0])
    bf = bf * size
    return bf

def tostring(bf):
    """ Prints bloom filter on one line.

    Args:
        bf: The bloom filter.
    """
    for x in bf:
        print(x, end='')
    print()
