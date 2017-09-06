"""Queries the database for specific genes."""

import sys

from bloom_filter import encode
import database

def main(f=None):
    """Reads in queries from a file and searches for them. If no file present,
    reads in quieries from the standard input and searches for them.

    Note: queries from standard in have a maximum of 1023 characters.

    Args:
        f: A text file containing queries. Each query should be on its own line.
    """

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
                print("Query: ", query_sequence.upper(), "\n")
                gene, iou = query(query_sequence)
                print("Best IOU: ", iou, "\n")
                print("Sequence: ", gene.sequence)
                print("---------------------------------------------\n")
    else:
        print("Enter query: ")
        for line in sys.stdin:
            query_sequence = line
            gene, iou = query(query_sequence)
            print("Best IOU: ", iou)
            print("Sequence: ", gene.sequence, "\n")
            print("Enter query: ")


def query(query_sequence):
    """Encodes a query and searches for it in the data base.

    Args:
        query: A genetic sequence (string) to be searched for.

    Returns:
        The 'Gene' that is the 'best match' to the query.

        The IOU for the 'best match' and the query.
    """

    print("encoding query...")
    query = encode(query_sequence)
    print("...query complete")

    print("performing search...")
    gene, iou = database.search(query)
    print("...search complete \n")

    return gene, iou

if __name__ == '__main__':
    if(len(sys.argv) > 1):
        main(sys.argv[1])
    else:
        main()
    sys.exit(0)
