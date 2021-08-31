+++
title = "Getting started with Bioconductor"

date = 2019-04-21T00:00:00
# lastmod = 2019-04-19T00:00:00

draft = false  # Is this a draft? true/false
toc = true  # Show table of contents? true/false
type = "docs"  # Do not modify.

# Add menu entry to sidebar.
linktitle = "Introduction to Bioconductor"
[menu.tutorial]
  parent = "ChIP-seq with R"
  weight = 2
+++

# Introduction

In this tutorial you will learn how to install and use bioconductor. 

[Bioconductor](https://www.bioconductor.org) is an opensource collection of R packages that provides a framework for doing bioinformatics in R. 

Install bioconductor with this:


```r
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install()
```

Afterwards, installing bioconductor (BioC) packages is a little different from other R packages and makes use of the BiocManager::install() function. For example, to install QuasR:


```r
BiocManager::install("QuasR", version = "3.8")
```

Bioconductor releases contain a number of R packages that have been designed to perform different tasks or provide the data structures required for interacting with genomic data.  

## Some examples of BioC packages

| Function  | Example packages  |
|-----------|-------------------|
Data structures | **IRanges, GenomicRanges, Biostrings, BSgenome**
Input of data   | **ShortRead**, Rsamtools, GenomicAlignments, **rtracklayer**
Annotation  | GenomicFeatures, **BSgenome, TxDb**
Alignment | **Rbowtie, QuasR**, Biostrings
ChIP-seq  | **ChIPseeker**, chipseq, ChIPseqR, **ChIPpeakAnno**, DiffBind, BayesPeak
De-novo motif discovery | rGADEM, MotifDb, SeqLogo, **ChIPpeakAnno**
RNA-seq | EdgeR

Packages in **bold** are ones that are used later.

# Creating and using GRanges objects

One of the most useful data structures is GRanges. These are essentially a list of genomic intervals that could be anything from genes to transcription factor binding sites.

## Creating simple GRanges


```r
library(GenomicRanges)
```

```
## Loading required package: stats4
```

```
## Loading required package: BiocGenerics
```

```
## Warning: package 'BiocGenerics' was built under R version 4.0.5
```

```
## Loading required package: parallel
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
##     union, unique, unsplit, which.max, which.min
```

```
## Loading required package: S4Vectors
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following object is masked from 'package:base':
## 
##     expand.grid
```

```
## Loading required package: IRanges
```

```
## Loading required package: GenomeInfoDb
```

```
## Warning: package 'GenomeInfoDb' was built under R version 4.0.5
```


To start, here is a GRanges object which has 3 genes, all on chromosome 1. The first gene runs from position 1-3, and are 3 nucleotides long.


```r
gr1 <- GRanges(seqnames = "chr1", strand = c("+", "-", "+"),
              ranges = IRanges(start = c(1,3,5), width = 3))
gr1
```

```
## GRanges object with 3 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chr1       1-3      +
##   [2]     chr1       3-5      -
##   [3]     chr1       5-7      +
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

To create a GRange for a single chromosome.


```r
chrI <- GRanges(seqnames = "chrI",
              ranges = IRanges(start = 1, width = 3000000))
chrI
```

```
## GRanges object with 1 range and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chrI 1-3000000      *
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

To create a GRange with ChIP-seq peaks, all on the same chromosome you could assign the start and end co-ordinates to vectors and create a GRange from that.


```r
peak_start <- c(100,220,450,767,899,1040)
peak_end <- c(140,260,490,800,945,1100)
peaks_gr <- GRanges(seqnames = "chrI",
                    ranges = IRanges(start=peak_start, end=peak_end))
peaks_gr
```

```
## GRanges object with 6 ranges and 0 metadata columns:
##       seqnames    ranges strand
##          <Rle> <IRanges>  <Rle>
##   [1]     chrI   100-140      *
##   [2]     chrI   220-260      *
##   [3]     chrI   450-490      *
##   [4]     chrI   767-800      *
##   [5]     chrI   899-945      *
##   [6]     chrI 1040-1100      *
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

## GRanges from dataframes

It is more likely that you will want to create GRanges objects from other data structures. For example, if you import peaks from a BED file. Bioconductor packages provide import functions for different filetypes and then this can be coerced into a dataframe and easily converted into a GRanges object with the `makeGRangesFromDataFrame` function.


