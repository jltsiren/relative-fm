# Relative FM-index

This is an implementation of the Relative FM-index and the Relative Compressed Suffix Tree. [The early implementatation](http://jltsiren.kapsi.fi/relative-fm) is available elsewhere.

To compile, install the [Succinct Data Structures Library 2.0](https://github.com/simongog/sdsl-lite) and set `SDSL_DIR` in the makefile to point to the SDSL directory. The implementation uses OpenMP, libstdc++ parallel mode, and C++11, so g++ 4.7 or newer is required.

More documentation will follow.

## References

Djamal Belazzougui, Travis Gagie, Simon Gog, Giovanni Manzini, and Jouni Sirén: **Relative FM-indexes**.
Proc. SPIRE 2014, Springer LNCS 8799, pp. 52-65, Ouro Preto, Brazil, October 20-22, 2014.
[DOI: 10.1007/978-3-319-11918-2_6](http://dx.doi.org/10.1007/978-3-319-11918-2_6)
