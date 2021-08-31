---
date: "2021-01-24"
linkTitle: ChIP-seq Analysis using R
summary: Some guidance for using R to do ChIP-seq analysis.
title: "ChIP-seq Analysis using R"
type: book
---

This tutorial aims to help analyse data from ChIP-seq experiments in bacteria that investigate the binding of transcription factors.

The tutorial will cover:

[__Introduction to Bioconductor__](/tutorial/bioconductor_intro/)  
- Installing Bioconductor and bioconductor packages  
- Accessing and creating genome packages  
- Interacting with genome data  

[__Aligning to a reference genome__](/tutorial/alignment/)  
- Aligning fastq files to your reference genome  
- Checking alignment stats and quality  

[__Normalising read counts__](/tutorial/normalisation/)  

[__Peak calling__ ](/tutorial/peakcalling/)  
- Calling peaks using MACS  
- Exporting and visualising peak data  

[__De-novo motif discovery__](/tutorial/motif/)  
- Finding transcription factor motifs from called peaks  

[__Visualisation__](/tutorial/visualisation/)  
- [Overview plots](/tutorial/circleplot/)  
- Distance from TSS  
- Functional analysis  
- Coverage plots  

## About

This was built in R markdown which is an easy-to-use, human readable markdown language. It is particularly useful because you can include R code whenever you like and throughout you will see sections of code so that you can a) see what code has been run and b) copy and use it yourself. 

I decided to bring together the tools that I have used to analyse ChIP-seq data because there is an abundance of information online for ChIP-seq analysis of eukaryotes, but not prokaryotes.