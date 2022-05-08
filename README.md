# RNAFolding

This project is a method to find the best way to fold an RNA, given its sequence.

## Nussinov algorithm

The method used in this project is the Nussinov algorithm, written in [Nussinov.py](Nussinov.py)

### Simple Nussinov

The simple Nussinov algorithm is well explained [here](https://en.wikipedia.org/wiki/Nussinov_algorithm)

### Multi Nussinov

Th' nussinov algorithm can be extended to be computed on a bunch of sequences instead of 2 sequences. It is written in [Nussinov.py](Nussinov.py), class `ACC_Computation`.

You can find 2 differents ways to do so : Maximum Expected Accuracy or Maximum total accurracy.