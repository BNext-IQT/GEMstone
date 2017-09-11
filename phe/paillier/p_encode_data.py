import json
import pickle as p
from collections import defaultdict, namedtuple
from p_bloom_filter import encode
import sys

# Gene data structure holding the name, principle investigator, a partial
# sequence and the corresponding bloom filter of a gene.
Gene = namedtuple("Gene", "name, pi, sequence, bloom")

# Type of sequence used in database
SEQUENCE_TYPE = "public_addgene_full_sequences"

def encode_data():
    """Reads in the addgene-plasmids-sequences data from a json file and stores
    the important information in an internal data structure. Also genenerates
    and stores the bloom filter for each gene in the database.
    """
    # Reads in data from JSON file
    try:
        with open("../addgene-plasmids-sequences.json") as data_file:
           print('encoding database...')
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

    p.dump(data, open('../encoded_addgene.p', 'wb'))
    print('...database encoding complete')
    
def main():
    encode_data()

if __name__ == '__main__':
    main()
    sys.exit(0)
