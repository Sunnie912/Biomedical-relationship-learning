# Biomedical-relationship-learning
A python implementation of the Ensemble Biclustering for Classification (EBC) algorithm. Although having "biclustering" in its name, EBC is a co-clustering algorithm that allows you to perform co-clustering on very large sparse N-dimensional matrix. For details and examples of using EBC please reference [this paper](http://www.ncbi.nlm.nih.gov/pubmed/26219079).

## Files
- 'Data': all the data we got from the Medline database. This file also includes before/after cleaning data.
- 'EBC': implementations of EBC algorithm using the sparse matrix, 2d NumPy dense array, and training/tests code with parsed data from 'Data' file.
- 'pubmed_parser': data cleaning using Pubmed Parser.
- 'standford_parser': data processing using Standford Parser.
- 'README.md' file.
