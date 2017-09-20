"""Queries the database for specific genes."""

import sys
import time

from joblib import Parallel, delayed
from phe import paillier
from p_bloom_filter import encode
from p_database import search, get_gene, Gene
from optimize_invert import invert

paillier.invert = invert

num_cores = 32 # Number of cores for parellel processing


####################
#
####################
def main(f=None):
    """Reads in queries from a file and searches for them. If no file present,
    reads in quieries from the standard input and searches for them.

    Note: queries from standard in have a maximum of 1023 characters.

    Args:
        f: A text file containing queries. Each query should be on its own line.
    """
    start = time.time()
    print('Start time: ' + str(start) + '\n')

    print('generating key pair...')
    public_key, private_key = paillier.generate_paillier_keypair()
    print('...key pair complete\n')

    
    if f:
        try:
            with open(f, 'r') as queries_file:
                text = queries_file.read()
                queries = text.splitlines()
            queries_file.close()
        except:
            print("Requires a text file with the queries.")
            sys.exit(2)

            
        for query_sequence in queries:
            if query_sequence:
                q_start = time.time()
                print("Query: ", query_sequence.upper(), "\n")
                gene, iou = query(query_sequence, public_key, private_key)
                print("Best IOU: ", iou, "\n")
                print("Sequence: ", gene.sequence)
                print("---------------------------------------------\n")
                q_end = time.time()
                q_elapsed = q_end - q_start
                print('Query run time: ' + str(q_elapsed) + '\n')
    else:
        print("Enter query: ")
        for query_sequence in sys.stdin:
            gene, iou = query(query_sequence, public_key, private_key)
            print("Best IOU: ", iou)
            print("Sequence: ", gene.sequence, "\n")
            print("Enter query: ")

            
    end = time.time()
    print('End time: ' + str(end))
    elapsed = end - start
    print('Time elapsed: ' + str(elapsed))


####################
#
####################
def query(query, public_key, private_key):
    """Encodes a query and searches for it in the data base.

    Args:
        query: A genetic sequence (string) to be searched for.
        public_key: The public key for the paillier encryption.
        private_key: The private key for the paillier encryption.

    Returns:
        The 'Gene' that is the 'best match' to the query.

        The IOU for the 'best match' and the query.
    """
    global num_cores

    print("encoding query...")
    query = encode(query)
    print("...encode complete")

    query_mag = magnitude(query)

    print("encrypting query...")
    #query = [public_key.encrypt(x) for x in query]
    query = Parallel(n_jobs=num_cores)(delayed(public_key.encrypt)(x) for x in query)
    
    print("...encrypt complete")
    print("generating scores...")
    
    scores = search(query)
    
    print("...scores complete")

    print("performing search...")
    max_iou = 0
    best_id = 0

    for id_ in scores:
        intersection = private_key.decrypt(scores[id_][0])
        score = iou(intersection, scores[id_][1], query_mag)
        if score > max_iou:
            max_iou = score
            best_id = id_

    gene = get_gene(best_id)

    print("...search complete \n")

    return gene, max_iou


####################
#
####################
def iou(intersection, data_mag, query_mag):
    """Finds the IOU for two bloom filters.

    Args:
        intersection: The intersection of the two genes.
        data_mag: The magnitude of the gene being compared to.
        query_mag: The magnitude of the gene being searched for.

    Returns:
        The IOU for the two genes.
    """
    union = (data_mag + query_mag) - intersection
    return intersection/union


####################
#
####################
def magnitude(v):
    """Finds the magnitude of a binary vector (array).

    Args:
        v: A binary vector (array).

    Returns:
        The magnitude of the vector.
    """
    #sum = 0
    #for x in v:
    #    sum += x
    #return sum
    return sum(v)


####################
# Main
####################
if __name__ == '__main__':
    if(len(sys.argv) > 1):
        main(sys.argv[1])
    else:
        main()
    sys.exit(0)
