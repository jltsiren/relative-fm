# Experiment logs for the RST paper

* **Note: These logs are for the submitted 2017 version of the paper.**
* "human" is the human reference genome.
* "female" is the human reference genome without chromosome Y.
* "maternal" and "paternal" are the corresponding haplotypes of NA12878.
* Numbers indicate copies of the same files.

## Index construction

* SSA, Full RFM: `build_*.log`
* Basic RFM: `old_*.log`, `old_*.err`
* RLCP: `rlcp_*.log`

## Component sizes

* RFM: Same as for construction
* RLCP: `breakdown_*.log`

## Benchmarks

* Basic queries: `verify_*.log`
* Find queries: `query_*.log`
* Locate queries: `locate_*.log`
  * The last number indicates SA sample interval.
* CST comparison: `cst_*.log`

## Synthetic collections

* `multi_*.log`, where the last number indicates mutation rate.
* Use `cst_size.py` to parse the logs.
* The comparison to GCT comes from the old benchmarks, as the construction of each GCT takes around 2.5 days.
