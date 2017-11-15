"""This is the primary code document for the execution of this workflow.
Queries the database for specific genes."""

# Load our packages for the environment
import sys
import time
from joblib import Parallel, delayed
from phe import paillier
from p_bloom_filter import encode
from p_database import search, magnitude
from optimize_invert import invert
from Bio import SeqIO

paillier.invert = invert
num_cores = 48 # Number of cores for parellel processing
query_len = 100

####################
# Main function to run pipeline
####################
def main(f, d, dev = False):
    """Reads in queries from a file and searches for them. If no file present,
    reads in quieries from the standard input and searches for them.

    Note: queries from standard in have a maximum of 1023 characters.

    Args:
        f: A FASTA file containing a query. 
        d: A path do a directory with FASTA files to act as the database to search
        dev: a boolean indicator, if True, runs the pipeline on only the first 
             entry in a data set rather than the whole set. 
    """
    global query_len
    
    start = time.time()
    
    print('\nStart time: ' + str(start) + '\n')
    print('generating key pair...')
    
    # Create the encryption public and private key pair
    public_key, private_key = paillier.generate_paillier_keypair()
    
    print('...key pair complete\n')
    
    
    # Build the query by concatenating the entries in a FASTA file together
    queries = ''
    
    with open(f, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            queries += str(record.seq)
        
        
    q_start = time.time()
    seq = queries
    seq = seq[:query_len]
    
    print("Query: ", seq.upper()[:1000], "\n")

    
    max_iou, max_ioLquery, max_ioLresult, best_seq, best_mag = query(seq, public_key, private_key, dev = dev, data_dir = d)
    
    
    q_end = time.time()
    q_elapsed = q_end - q_start

    
    print('Query run time: ' + str(q_elapsed) + '\n')
    print('Length of result: %s' % str(len(best_seq)))
    print('Length of query: %s' % str(len(seq)), '\n')
    print("Best IoU: ", max_iou)
    print("Best IoLenQuery: ", max_ioLquery)
    print("Best IoLenResult: ", max_ioLresult, '\n') 
    
    print("Sequence: ", best_seq[:1000])
    print("---------------------------------------------\n")    
            
    end = time.time()
    print('End time: ' + str(end))
    elapsed = end - start
    print('Time elapsed (min): ' + str(float(elapsed)/60))
    
    out = str([max_iou, max_ioLquery, max_ioLresult, q_elapsed, seq[:2000], best_seq][:2000])+'\n'
    with open("results.txt", "a") as myfile:
        myfile.write(out)

####################
# Query a database with a query and public key and decrypt using a private key
####################
def query(query, public_key, private_key, dev, data_dir):
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
    print("...encode complete\n")

    query_mag = magnitude(query)

    
    # Encrypting query
    print("encrypting query...")
    
    encrypt_start = time.time()
    
    query = Parallel(n_jobs=num_cores)(delayed(public_key.encrypt)(x) for x in query)

    encrypt_end = time.time()
    
    print("...encrypt complete: Encrypt time (min) = %s" % str(float(encrypt_end - encrypt_start)/60))
    print("generating scores...")
    
    scores = search(query, data_dir = data_dir)
    
    print("...scores complete")
    print("performing search...")
    
    max_iou = 0
    best_id = 0
    best_seq = ''
    result_scores = Parallel(n_jobs=num_cores)(delayed(calc_iou)(id_, private_key, query_mag) for id_ in scores)
    
    for score_set in result_scores:
        if score_set[0] >= max_iou: 
            max_iou = score_set[0]
            max_ioLquery = score_set[1]
            max_ioLresult = score_set[2]  
            best_seq = score_set[3]
            result_mag = score_set[4]
            
    print("...search complete \n")
    
    return max_iou, max_ioLquery, max_ioLresult, best_seq, result_mag


####################
# Calculate the Intersection over Union
####################
def calc_iou(id_, private_key, query_mag):
    intersection = private_key.decrypt(id_[0])
    Iou, IoLquery, IoLresult = iou(intersection, id_[1], query_mag)
    
    return Iou, IoLquery, IoLresult, id_[2], id_[1]


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
    main(sys.argv[1], sys.argv[2])
    sys.exit(0)