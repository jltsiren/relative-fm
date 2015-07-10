# Relative FM-index

This is an implementation of the Relative FM-index (RFM) and the Relative Compressed Suffix Tree (RCST). [The early implementatation](http://jltsiren.kapsi.fi/relative-fm) is available elsewhere.

To compile, install the [Succinct Data Structures Library 2.0](https://github.com/simongog/sdsl-lite) and set `SDSL_DIR` in the makefile to point to the SDSL directory. The implementation uses OpenMP, libstdc++ parallel mode, and C++11, so g++ 4.7 or newer is required.

## Compilation options

The maximum resident size reported by `getrusage()` is in kilobytes in Linux and in bytes in OS X. By default, the implementation assumes Linux-like behavior. To get the correct memory usage reports in OS X, uncomment the line `RUSAGE_FLAGS=-DRUSAGE_IN_BYTES` in the makefile.

By default, the test programs will produce minimal output. To change this, add `-DVERBOSE_OUTPUT` and/or `-DVERBOSE_STATUS_INFO` to `VERBOSE_FLAGS` in the makefile. The first option will add more detailed output to `stdout`, while the second option prints additional status information to `stderr`.

Other flags may or may not work correctly if changed from the defaults.

## Test programs

The following programs generally depend on the output produced by the earlier ones in the list:

* `mutate` creates synthetic document collections, as described in the paper. Example: `mutate source target 0.001`.
* `build_bwt` builds the BWT for one or more input files. It can also be used to build the LCP array and to sample SA/ISA.
* `index_dlcp` builds an index for the DLCP array of the reference, given its LCP array. The index is used when building the RLCP array.
* `align_bwts` builds RFM and RCST. The relevant options are `-i` (build a full RFM instead of a basic RFM) and `-L` (also build the RLCP array).
* `query_test` queries one or more indexes with the given pattern file, where the patterns are listed one per line. The most relevant index types are `sp` (SSA) and `rpp` (RFM index built with `align_bwts`). Tag `b` allows using other RFM types, while tags `L` and `L V` change the type of queries executed.
* `cst_traverse` and `cst_compare` were used for the CST experiments (traversal and matching statistics) in the paper. These programs will reuse the SDSL compressed suffix trees and relative select structures built during earlier executions.

These programs are technical and/or obsolete:

* `test_rlz` and `verify` can be used for verifying and/or timing various RLZ options and queries.
* `build_rlzfm` and `lcs` have been abandoned and may not work at all.

"Base name" is the name of the input file, "ref" is the base name of the reference sequence, "seq" is the base name of a target sequence, and "paper" refers to the RCST paper in subdirectory `rcst`. When different structures are stored on disk, extensions that depend on the type of the structure are appended to the base name.

## Todo

* RLZ pointer encoding:
  * Absolute pointers require a lot of space.
  * Run-length encoded relative pointers do not take any less space in the RLCP array.
  * Differential encoding the relative pointers could work, especially if we find the longest increasing subsequence of the absolute pointers and store the rest separately.
  * Non-greedy parsing could also help.
* Optimizations:
  * In `RelativeFM::inverse()`, position *i-1* can provide faster access to *ISA[i]*, if it is an lcs-position.
  * In `increasingSubsequence()`, we could use binary search (or an index) in the left/right match arrays, instead of storing the values in arrays `smallest` and `previous`.
  * RLCP construction may not need separate LCP and DLCP arrays.
  * Left and right matches do not have to be adjacent to the current suffix in the mutual suffix array. The original definition used the nearest suffix (of the target sequence) on one side and the adjacent suffix on the other side. Could we use the nearest suffix on both sides?
  * Specialized `getComplement()` for `RLSequence`.
  * Could it be more space-efficient to store the leaf values in the RLCP minimum tree relative to their parents?
  * RLZ parsing that streams the text from disk.
* Alternate solutions:
  * RLCP based on LCP instead of DLCP can be smaller.
  * `hyb_vector` could work well in LCS bitvectors, once it supports select queries.
  * Finding bwt-invariant subsequences using SA instead of CSA.
* Implementation:
  * Use class-specific `size_type` instead of `uint64_t` when possible.
  * Add copy constructors, assignment operators etc. to the classes.
* RLZFM:
  * Construction from any `SimpleFM`.
  * Mismatch vector could store comp values instead of character values.
  * Can the SA samples be compressed with RLZ? Maybe if we use text order sampling.
* RLZ bitvector:
  * Not all structures are required for all queries. Try moving them to rank/select support structures.
  * Implement full decompression.

## References

Djamal Belazzougui, Travis Gagie, Simon Gog, Giovanni Manzini, and Jouni Sirén: **Relative FM-indexes**.
Proc. SPIRE 2014, Springer LNCS 8799, pp. 52-65, Ouro Preto, Brazil, October 20-22, 2014.
[DOI: 10.1007/978-3-319-11918-2_6](http://dx.doi.org/10.1007/978-3-319-11918-2_6)

Christina Boucher, Alexander Bowe, Travis Gagie, Giovanni Manzini, and Jouni Sirén: **Relative Select**.
Accepted to SPIRE 2015.
