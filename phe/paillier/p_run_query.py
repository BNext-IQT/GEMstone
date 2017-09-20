"""Queries the database for specific genes."""

import sys
import time
import pickle as p
import numpy as np
from phe import paillier
from p_querier import query_database, iou, get_gene, magnitude
from p_database import encrypt_bloom_filter
from p_bloom_filter import build_bloom_filter
from joblib import Parallel, delayed

### Main function to run our query system ###
def main(f):
    """Reads in pre-prepared queries from a pickle file and searches for them. 

    Args:
        f: A pickle file containing queries. 
    """
    start = time.time()
    print('Start time: ' + str(start) + '\n')
    
    print('Initializing database...')
    data = p.load(open('prepared_encoded_database.p', 'rb'))
    print("...database loaded")
    
    print('Importing query and encryption information...\n')
    
    # Load the data
    with open(f, 'r') as queries_file:
        text = queries_file.read()
        queries = text.splitlines()
    queries_file.close()
    
    qs = []
    SIZE = 12000
    K = 16
    H = hash

    # Build query bloom filters
    for query in queries:
        qs.append(build_bloom_filter(query, SIZE, K, H))
    
    
    if True:
    #try:
        '''
        information = p.load(open(f, 'rb'))

        public_key = information['public_key']
        private_key = information['private_key']
        query_list = information['query_list']
        query_blooms = information['query_blooms']
        encrypted_queries = information['encrypted_blooms']
        '''
        #for (query, query_bloom, enc_query) in zip(query_list, query_blooms, encrypted_queries):
        
        for query_bloom,query in zip(qs,queries):
            q_start = time.time()
            print("Query: ", query.upper(), "\n")
                
            print('Encrypting query...\n')
            
            public_key, private_key = paillier.generate_paillier_keypair()
            enc_query = encrypt_bloom_filter(query_bloom, public_key)
            query_mag = sum(query_bloom)
            
            print('...query encrypted\n')
            
            #query_mag = sum(query_bloom)

            scores = query_database(enc_query, data)
            
            max_iou = 0
            best_id = 0

            for id_ in scores:
                intersection = private_key.decrypt(scores[id_][0])
                score = iou(intersection, scores[id_][1], query_mag)
                if score > max_iou:
                    max_iou = score
                    best_id = id_

            gene = get_gene(best_id, data)

            print("...search complete \n")

            
            
            #unencrypted_inter = Parallel(n_jobs=32)(delayed(private_key.decrypt)(entry[0]) for entry in results)
            #magnitudes = [entry[1] for entry in results]
            
            #IoUs = Parallel(n_jobs=32)(delayed(iou)(inter, mag, query_mag) for inter,mag in zip(unencrypted_inter, magnitudes))
            
            #IoU_max = max(IoUs)
          
            print("Best IOU: ", max_iou, "\n")
            print("Sequence: ", gene.sequence, '\n')
            
          
            q_end = time.time()
            q_elapsed = q_end - q_start
          
            print('Query run time: ' + str(q_elapsed) + '\n')
            print("---------------------------------------------\n")
            return 0 # This is just to stop the loop for debugging purposes
          
    #except:
    #    print("Requires a pickle file with the queries.")
    #    sys.exit(2)

              
    end = time.time()
    print('End time: ' + str(end))
    elapsed = end - start
    print('Time elapsed: ' + str(elapsed))
          
          
if __name__ == '__main__':
    if(len(sys.argv) > 0):
        main(sys.argv[1])
    else:
        print('Please included a query information file')
    sys.exit(0)
