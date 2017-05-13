# Experiment logs for the RCST paper

* **Note: These logs are for the old arXiv version of the paper from 2015.**
* "human" is the human reference genome.
* "female" is the human reference genome without chromosome Y.
* "maternal" and "paternal" are the corresponding haplotypes of NA12878.
* Prefixes and numbers indicate copies of the same files.

## Index construction and index sizes

* Basic RFM: `old_rfm_*.log`
* SSA, Full RFM: `build_rfm_*.log`
* RCST: `test_rfm_*.log`

## Component sizes

* Basic RFM: `old_rfm_*.log`
* Full RFM, RLCP: `test_rfm_*.log`

## Query tests

* LF/Psi in SSA, RFM: `verify_psi_*.log`
* LCP, RLCP: `verify_lcp_*.log`
* Locate queries: `locate_test_*.log`

## Synthetic datasets

* Varying mutation rate: `index_synth_*.log`
* Comparison to GCT: `index_multi_*.log`, `index_multi_*.err`
  * Use `cst_size.py` to parse the `.log` files.

## Compressed suffix trees

* Full traversal: `cst_traversal_*.log`
* Matching statistics: `cst_compare_*.log`
