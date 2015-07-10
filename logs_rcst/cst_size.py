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
      if state == 0:
        res = size.search(line)
        if res:
          bytes += int(res.group(1))
        res = fm.search(line)
        if res:
          fm_size = float(res.group(1))
          state = 1
      elif state == 1:
        res = lcp.search(line)
        if res:
          lcp_size = float(res.group(1))
          state = 2
      else:
        res = rfm.search(line)
        if res:
          rfm_size += float(res.group(1))
        res = rlcp.search(line)
        if res:
          rlcp_size += float(res.group(1))
        res = repet.search(line)
        if res:
          repet_size = int(res.group(1)) / 1048576.0

    infile.close()
    total_size = fm_size + lcp_size + rfm_size + rlcp_size

    print "File:", sys.argv[arg]
    print "Original:", (bytes / 1048576.0), "MB"
    printSize("Reference", fm_size + lcp_size, bytes)
    printSize("Sequences", rfm_size + rlcp_size, bytes)
    printSize("Relative CST", total_size, bytes)
    printSize("Repetitive CST", repet_size, bytes)
    print

if __name__ == "__main__":
  main()
