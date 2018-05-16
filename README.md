# GEMstone - SIG-DB

SIG-DB is an algorithm developed to search privately held genomic databases using homomorphic encryption and locality-sensitive hashing. The primary implementation uses partially homomorphic encryption via the Paillier cryposystem and the fully homomorphic systems implements an Python wrapper version of the SEAL encryption library, developed by Microsoft Research. 

The algorithm code is hosted in the Python notebook named [main_define_classes.ipynb](https://github.com/BNext-IQT/GEMstone/blob/master/phe/paillier/main_define_classes.ipynb). There are two code blocks that execute the PHE and FHE versions of the algorithms, and classes defined for each algorithm below. 

If you have any questions, [contact us for questions](https://www.bnext.org/contact/).

**Manuscript Abstract:** Genomic data are becoming increasingly valuable as we develop methods to utilize the information at scale and gain a greater understanding of how genetic information relates to biological function. Advances in synthetic biology and the decreased cost of sequencing are increasing the amount of privately held genomic data. As the quantity and value of private genomic data grows, so does the incentive to acquire and protect such data, which creates a need to store and process these data securely. We present an algorithm for the Secure Interrogation of Genomic DataBases (SIG-DB). The SIG-DB algorithm enables databases of genomic sequences to be searched with an encrypted query sequence without revealing the query sequence to the Database Owner or any of the database sequences to the Querier. SIG-DB is the first application of its kind to take advantage of locality-sensitive hashing and homomorphic encryption to allow generalized sequence-to-sequence comparisons of genomic data. 

- [Blog Description and Early Results](https://medium.com/bioquest/gemstone-series-secure-information-sharing-for-genetic-queries-early-results-a14e79bb4f57)
- [SIG-DB on arXiv](https://arxiv.org/abs/1803.09565)
- [PySEAL on arXiv](https://arxiv.org/abs/1803.01891)
