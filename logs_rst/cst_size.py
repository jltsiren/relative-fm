#! /usr/bin/env python
# -*- coding: utf-8 -*-

# cst_size.py logfile1 [logfile2 ...]

# This computes the total size of a collection of Relative CSTs, including
# the size of the reference.

import sys, re

def printSize(type, megabytes, bytes):
  print (type + ":"), megabytes, "MB,", (megabytes * 1048576.0 * 8 / bytes), "bpc"

def main():

  if len(sys.argv) < 2:
    return

  repet_sizes = {}
  repet_sizes["0.0001"] = 101453735
  repet_sizes["0.0003"] = 127063602
  repet_sizes["0.001"] = 158264566
  repet_sizes["0.003"] = 199825193
  repet_sizes["0.01"] = 240325131
  repet_sizes["0.03"] = 250474638
  repet_sizes["0.1"] = 230936675

  rate_pattern = "Mutation rate:\s(.*)"
  rate = re.compile(rate_pattern)

  size_pattern = "Target size:\s*(.*)"
  size = re.compile(size_pattern)

  fm_pattern = "Simple FM:\s*(.*) MB"
  fm = re.compile(fm_pattern)

  lcp_pattern = "LCP array:\s*(.*) MB"
  lcp = re.compile(lcp_pattern)

  rfm_pattern = "Relative FM:\s*(.*) MB"
  rfm = re.compile(rfm_pattern)

  rlcp_pattern = "Relative LCP:\s*(.*) MB"
  rlcp = re.compile(rlcp_pattern)

  repet_pattern = "Total:\s*(.*)"
  repet = re.compile(repet_pattern)

  for arg in range(1, len(sys.argv)):
    fm_size = 0.0
    lcp_size = 0.0
    rfm_size = 0.0
    rlcp_size = 0.0
    total_size = 0.0
    repet_size = 0.0
    bytes = 0
    infile = open(sys.argv[arg], "rb")
    state = 0

    for line in infile:
      if state == 0:   # Mutation rate
        res = rate.search(line)
        if res:
          if res.group(1) in repet_sizes:
            repet_size = repet_sizes[res.group(1)] / 1048576.0
          state = 1
      elif state == 1:   # Collection size
        res = size.search(line)
        if res:
          bytes += int(res.group(1))
        res = fm.search(line)
        if res:
          fm_size = float(res.group(1))
          state = 2
      elif state == 2: # RFM size
        res = rfm.search(line)
        if res:
          rfm_size += float(res.group(1))
        res = lcp.search(line)
        if res:
          lcp_size = float(res.group(1))
          state = 3
      elif state == 3: # RLCP size
        res = rlcp.search(line)
        if res:
          rlcp_size += float(res.group(1))

    infile.close()
    total_size = fm_size + lcp_size + rfm_size + rlcp_size

    print "File:", sys.argv[arg]
    print "Original:", (bytes / 1048576.0), "MB"
    printSize("Reference", fm_size + lcp_size, bytes)
    printSize("RFM", rfm_size, bytes)
    printSize("RLCP", rlcp_size, bytes)
    printSize("RFM+RLCP", rfm_size + rlcp_size, bytes)
    printSize("Relative CST", total_size, bytes)
    printSize("Repetitive CST", repet_size, bytes)
    print

if __name__ == "__main__":
  main()
