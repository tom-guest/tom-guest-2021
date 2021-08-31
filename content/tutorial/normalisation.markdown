+++
title = "Normalising read counts"

date = 2019-04-19T00:00:00
# lastmod = 2019-04-19T00:00:00

draft = false  # Is this a draft? true/false
toc = true  # Show table of contents? true/false
type = "docs"  # Do not modify.

# Add menu entry to sidebar.
linktitle = "Normalising read counts"
[menu.tutorial]
  parent = "ChIP-seq with R"
  weight = 5
+++

## Introduction

In this tutorial we will use the normalise read counts. 

To do this you will need to import the BAM files from your alignments and use the mapping statistics and the average read length from the sequencing run.

