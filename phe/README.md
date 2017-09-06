# BioML

## Unencryted
The unencrypted directory demonstrates searching of a genetic database using bloom filters and IOU comparison in an unencrypted space. There are two separate modules, the *querier* and the *database*. Both rely on a third module which encodes genes as bloom-filters. This bloom-filter module allows for the hash, k-mer length and bloom-filter size to be adjusted based on specs of the data. The querier module encodes a genetic sequence query as a bloom-filter and sends it to the database. The database then uses an intersection over union (IOU) comparison to determine the gene in the database that is the best match and returns that gene, as well as the corresponding IOU score.

The database is currently set up to encode and query the addgene-plasmids-sequences data set. This data set is too large to load in Git Hub but must be in the BioML directory to use the project.

The query module can either read queries from a text file or the command line. If a text file is included, it will query each separate line in the file. A sample text file is included in the repository. Otherwise, the user will be prompted to enter queries directly into the commmand line.

To run searches from the command line:
```shell
python querier.py
```

To run searches from a text file:
```shell
python querier.py sample.txt
```

Note that the builtin python hash function chooses a random seed each run so results with low certainty may change due to false positives. To keep the seed consistent use:
```shell
PYTHONHASHSEED=0 python querier.py sample.txt
```

## Paillier
The paillier directory attempts to replicate the unencrypted directory in an encrypted space using the partially homomorphic encryption algorithm paillier. Paillier is homomorphic over addition, which does not allow for a direct computation of the IOU score for two genes in the encrypted space. However, the components - the intersection and the union of the two genes - can be computed in the encrypted space. Thus, this implementation instead returns the encrypted intersection and union of the query and each of the genes in the database back to the querier where these components can be decrypted and used to find the best IOU score and corresponding gene.

This implementation relies on the python package *phe* which is a very basic implementation of paillier encryption which is not particularly efficient. Some parts of the implementation have been monkey-patched in the *optimize_invert* module to improve performance.

This can be run the same way as the unencrypted search. Note the querier module is called *p_callier.py* so the following command would be used.

```shell
PYTHONHASHSEED=0 python p_querier.py sample.txt
```

## Seal
The seal directory begins to explore a fully homomorphic encryption algorithm implemented by Microsoft. The implementation has been offered to run on MacOS and there is a small sample demonstrating encryption and decryption as well as the calculation of a one-bit max in the *GeneEncryption* directory. The make in this directory produces an executable in the *bin* directory called *gene*.

The next step in exploring this library would be to create a larger max function that could compare to IOU scores. Included below is sudo-code. However, I could not figure out how to separate the encrypted numbers into bits. While it is possible to divide by a constant, I could not figure out how to do integer division, although it is likely possible.

``` cpp
Ciphertext max(Ciphertext x, Ciphertext y, int n)
{
  if (n == 1)
  {
    return bit_max(x, y)
  }
  else
  {
    // x_n, y_n are the nth digit of x, y
    // x_rest, y_rest are the (n-1)th...0th digits of x, y
    return 2^n * bit_max(x_n, y_n)
        + (x_n * y_n) * max(x_rest, y_rest)
        + x_n * (1 - y_n) * x_rest
        + (1- x_n) * y_n * y_rest
        + (1 - x_n) * (1 - y_n) * max(x_rest, y_rest);
  }
}
```

## docker
A docker file, as well as a build and run script, are included to test both the unencrypted and paillier encrypted search.

## TODO:
* return top *n* IOU scores instead of just top 1 (use heap)
* serialize database so it does not have to be encoded each run
* continue SEAL research
