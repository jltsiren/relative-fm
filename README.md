# Relative FM-index

This is an implementation of the Relative FM-index (RFM) and the Relative (Compressed) Suffix Tree (RST, RCST). [The early implementatation](https://jltsiren.kapsi.fi/relative-fm) is available elsewhere.

The implementation uses OpenMP, libstdc++ parallel mode, and C++11. It should compile with g++ 4.9 or later.

# Dependencies

The implementation uses [SDSL](https://github.com/simongog/sdsl-lite) extensively. Set `SDSL_DIR` in the makefile to point to the SDSL directory.

The LCP array is compressed using the [lcp_library branch of RLZAP](https://github.com/farruggia/rlzap/tree/lcp_library). Set `RLZ_USE_EXTERNAL_SDSL` to `ON` in the `CMakeLists.txt` of RLZAP, and set `RLZAP_DIR` in the makefile of this project to point to the RLZAP directory.


## Compilation options

The maximum resident size reported by `getrusage()` is in kilobytes in Linux and in bytes in OS X. By default, the implementation assumes Linux-like behavior. To get the correct memory usage reports in OS X, uncomment the line `RUSAGE_FLAGS=-DRUSAGE_IN_BYTES` in the makefile.

By default, the test programs will produce minimal output. To change this, add `-DVERBOSE_OUTPUT` and/or `-DVERBOSE_STATUS_INFO` to `VERBOSE_FLAGS` in the makefile. The first option will add more detailed output to `stdout`, while the second option prints additional status information to `stderr`.

Other flags may or may not work correctly if changed from the defaults.

## Test programs

The following programs generally depend on the output produced by the earlier ones in the list:

* `mutate` creates synthetic document collections, as described in the paper. Example: `mutate source target 0.001`.
* `build_bwt` builds the BWT and samples the suffix array. It will also build the LCP array with option `-l`.
* `align_bwts` builds the RFM index. Option `-i` can be used to build the full RFM supporting locate/extract queries.
* `build_rlcp` builds the RLCP.
* `query_test` queries one or more indexes with the given pattern file, where each non-empty line is used as a pattern. The most relevant index types are `sp` (SSA) and `rpp` (RFM). Tags `L` switches to from find queries to locate queries.
* `verify` can be used for benchmarking the basic queries.
* `cst_traverse` and `cst_compare` were used for the CST experiments (traversal and maximal matches) in the paper. These programs will reuse the SDSL compressed suffix trees and relative select structures built during earlier executions.

**Base name** is the name of the input file, **ref** is the base name of the reference sequence, **seq** or **target** is the base name of a target sequence, and **paper** refers to the RCST paper in subdirectory `rcst`. When the structures are stored on disk, extensions that depend on the type of the structure are appended to the base name.

## Datasets

The following datasets have been used in the articles based on the relative FM-index:

* The 1000 Genomes Project version of the GRCh37 assembly of the human reference genome: [ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/](ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/)
* The genome of a Han Chinese individual from the YanHuang project: [ftp://public.genomics.org.cn/BGI/yanhuang/fa/](ftp://public.genomics.org.cn/BGI/yanhuang/fa/)
* The genome of the 1000 Genomes Project individual NA12878: [http://sv.gersteinlab.org/NA12878_diploid/](http://sv.gersteinlab.org/NA12878_diploid/)
* Read set `ERR019904_1` of the 1000 Genomes Project individual HG00122: [ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/HG00122/sequence_read/](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/HG00122/sequence_read/)

## Future development

* Optimizations:
  * In `RelativeFM::inverse()`, position *i-1* can provide faster access to *ISA[i]*, if it is an lcs-position.
  * In `increasingSubsequence()`, we could use binary search (or an index) in the left/right match arrays, instead of storing the values in arrays `smallest` and `previous`.
  * `RelativeCST::select_leaf()`: We can usually retrieve `LCP[i]` and `LCP[i+1]` using the same iterator.
  * Left and right matches do not have to be adjacent to the current suffix in the mutual suffix array. The original definition used the nearest suffix (of the target sequence) on one side and the adjacent suffix on the other side. Could we use the nearest suffix on both sides?
* Implementation:
  * Add copy constructors, assignment operators etc. to the classes.
  * Switch to SDSL structures as references.
* Alternative solutions:
  * `hyb_vector` could work well in LCS bitvectors, once it supports select queries.
  * Gap encoded / run-length encoded bitvectors with Huffman coded gap lengths could be useful somewhere (see Nicola Prezza: A Compressed-Gap Data-Aware Measure for Indexable Dictionaries).
  * Finding bwt-invariant subsequences using SA instead of CSA.
  * A basic RFM index for *S* can be build directly by backward searching for *S* in the reference index.
  * RFM index relative to GCSA.
  * Relative CST-Sada: RLZ bitvectors for suffix tree topology and PLCP.
  * A relative version of the range min-max tree used in the CST-npr for repetitive collections (see Andrés Abeliuk, Rodrigo Cánovas, and Gonzalo Navarro: Practical Compressed Suffix Trees).

## References

Djamal Belazzougui, Travis Gagie, Simon Gog, Giovanni Manzini, and Jouni Sirén: **Relative FM-indexes**.
Proc. SPIRE 2014, Springer LNCS 8799, pp. 52-65, Ouro Preto, Brazil, October 20-22, 2014.
[DOI: 10.1007/978-3-319-11918-2_6](https://doi.org/10.1007/978-3-319-11918-2_6)

Christina Boucher, Alexander Bowe, Travis Gagie, Giovanni Manzini, and Jouni Sirén: **Relative Select**.
Proc. SPIRE 2015, Springer LNCS 9309, pp. 149-155, London, UK, September 1-3, 2015.
[DOI: 10.1007/978-3-319-23826-5_15](https://doi.org/10.1007/978-3-319-23826-5_15)

Anthony J. Cox, Andrea Farruggia, Travis Gagie, Simon J. Puglisi, and Jouni Sirén: **RLZAP: Relative Lempel-Ziv with Adaptive Pointers**.
Proc. SPIRE 2016, Springer LNCS 9954, pp. 1-14, Beppu, Japan, October 18-20, 2016.
[DOI: 10.1007/978-3-319-46049-9_1](https://doi.org/10.1007/978-3-319-46049-9_1)

Andrea Farruggia, Travis Gagie, Gonzalo Navarro, Simon J. Puglisi, and Jouni Sirén: **Relative Suffix Trees**.
Accepted to The Computer Journal, 2017.
[arXiv:1508.02550](https://arxiv.org/abs/1508.02550)
