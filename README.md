# BLAST_gapalign

This is a Rust port of the BLAST ALIGN_EX function for gapped alignment 
extension. It is based on the original C code found in the 2.14.1 release 
of NCBI BLAST with minimal modifications to simplify the data structures.
This is **not** a full port of the BLAST algorithm, which involves many more
components (seed finding, ungapped extension, statistics generation, sequence
management etc). However, it provides the same means of generating a greedy
alignment of two subsequences given an anchor (seed), a scoring system 
(matrix, and affine gap penalties) and a x-drop threshold for termining the 
extension of the alignment.


## Features
- A simple API for greedy seeded alignment extension.
- Support for custom scoring matrices and affine gap penalties.

## Example


## References

Original C code: <https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.14.1/ncbi-blast-2.14.1+-src.tar.gz>
file: ncbi-blast-2.14.1+-src/c++/src/algo/blast/core/blast_gapalign.c

Zhang, Zheng, et al. 
"*A greedy algorithm for aligning DNA sequences.*" 
Journal of Computational biology 7.1-2 (2000): 203-214.


