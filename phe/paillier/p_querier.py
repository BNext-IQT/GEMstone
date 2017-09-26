"""Queries the database for specific genes."""

import sys
import time

from joblib import Parallel, delayed
from phe import paillier
from p_bloom_filter import encode
from p_database import search, get_gene, Gene, magnitude
from optimize_invert import invert
from Bio import SeqIO

paillier.invert = invert

num_cores = 32 # Number of cores for parellel processing


####################
# Main function to run pipeline
####################
def main(f=None):
    """Reads in queries from a file and searches for them. If no file present,
    reads in quieries from the standard input and searches for them.

    Note: queries from standard in have a maximum of 1023 characters.

    Args:
        f: A text file containing queries. Each query should be on its own line.
    """
    start = time.time()
    print('\nStart time: ' + str(start) + '\n')

    print('generating key pair...')
    public_key, private_key = paillier.generate_paillier_keypair()
    print('...key pair complete\n')

    
    if f:
        try:
        #if True:
            #with open(f, 'r') as queries_file:
            #    text = queries_file.read()
            #    queries = text.splitlines()
            #queries_file.close()
            
            #queries = list(SeqIO.parse(f, "fasta"))
            
            queries = ''

            with open(f, "r") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    queries += str(record.seq)
        
        except:
            print("Requires a text file with the queries.")
            sys.exit(2)

            
        #for query_sequence in queries:
        #if query_sequence:
        if queries:
                q_start = time.time()
                #seq = str(query_sequence.seq)
                seq = queries
                print("Query: ", seq.upper()[:1000], "\n")
                
                #gene, max_iou, max_ioLquery, max_ioLresult = query(seq, public_key, private_key)
                max_iou, max_ioLquery, max_ioLresult = query(seq, public_key, private_key)
                                
                q_end = time.time()
                q_elapsed = q_end - q_start
                
                print('Query run time: ' + str(q_elapsed) + '\n')
                print("Best IoU: ", max_iou, '\n')
                print("Best IoLenQuery: ", max_ioLquery, '\n')
                print("Best IoLenResult: ", max_ioLresult, '\n') 
                #print("Sequence: ", gene.sequence, "\n")
                
                print("---------------------------------------------\n")    
    else:
        print("Enter query: ")
        for query_sequence in sys.stdin:
            gene, max_iou, max_ioLquery, max_ioLresult = query(query_sequence, public_key, private_key)
            print("Best IoU: ", max_iou, '\n')
            print("Best IoLenQuery: ", max_ioLquery, '\n')
            print("Best IoLenResult: ", max_ioLresult, '\n')            
            print("Sequence: ", gene.sequence, "\n")
            print("Enter query: ")

            
    end = time.time()
    print('End time: ' + str(end))
    elapsed = end - start
    print('Time elapsed: ' + str(elapsed))


####################
# Query a database with a query and public key and decrypt using a private key
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
    print('Length of query BF: ' + str(len(query)))
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
    
    result_scores = Parallel(n_jobs=num_cores)(delayed(calc_iou)(id_, private_key, query_mag) for id_ in scores)
    
    for score_set in result_scores:
        if score_set[0] > max_iou:
            #max_iou = score_set[0]
            #best_id = score_set[1]
            #max_ioLquery = score_set[2]
            #max_ioLresult = score_set[3]   
            max_iou = score_set[0]
            max_ioLquery = score_set[1]
            max_ioLresult = score_set[2]  
    
    #gene = get_gene(best_id)

    print("...search complete \n")

    #return gene, max_iou, max_ioLquery, max_ioLresult
    return max_iou, max_ioLquery, max_ioLresult


####################
# Calculate the Intersection over Union
####################
def calc_iou(id_, private_key, query_mag):
    intersection = private_key.decrypt(id_[0])
    result_iou, max_ioLquery, max_ioLresult = iou(intersection, id_[1], query_mag)
    
    #return result_iou, id_[2], max_ioLquery, max_ioLresult
    return result_iou, max_ioLquery, max_ioLresult


####################
# Calculate the Intersection over Union
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
    
    iou = intersection/union
    max_ioLquery = intersection/query_mag
    max_ioLresult = intersection/data_mag
    
    return iou, max_ioLquery, max_ioLresult


####################
# Main
####################
if __name__ == '__main__':
    if(len(sys.argv) > 1):
        main(sys.argv[1])
    else:
        main()
    sys.exit(0)
