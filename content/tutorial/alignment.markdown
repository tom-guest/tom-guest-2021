+++
title = "Aligning to a reference genome"

date = 2019-04-19T00:00:00
# lastmod = 2019-04-19T00:00:00

draft = false  # Is this a draft? true/false
toc = true  # Show table of contents? true/false
type = "docs"  # Do not modify.

# Add menu entry to sidebar.
linktitle = "Aligning to a reference genome"
[menu.tutorial]
  parent = "ChIP-seq with R"
  weight = 4
+++

## Introduction

In this tutorial we will use the [QuasR](https://bioconductor.org/packages/release/bioc/html/QuasR.html) package to align our sequencing data to a reference genome.
 
QuasR uses the alignment program `Bowtie`, and will produce `bam` files for each alignment. There are more advanced options for including auxiliary genomes (which are used to align 'leftover' unmapped sequences) which are useful to check for contaminating DNA and for spiked experiments. 

### Sequencing reads

In ChIP-seq the sequencing reads are typically short single end dsDNA. This means that the 5’ end will be sequenced on “+” strand and the 3’ end will be on “-” strand. "+” reads extend only in positive direction and “-” reads in negative direction which results in the typical bimodal peak at transcription factor binding sites.

![](/tutorial/alignment_files/reads.png)

## Install the QuasR package


```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("QuasR", version = "3.8")
```

## Folder structure

I recommend having a folder for each new project. Within this create separate folders for `raw_fastq`, `alignments`, `genomes`, etc. You can then keep things tidy and avoid making any changes or accidentally deleting your raw fastq sequencing data.

You should set the project root as working directory, (or even better create a RStudio project).


```r
setwd("path/to/myChIPseqProject")
```

## Preparation

Download and save the genome `fasta` file(s) to your genomes folder. Then direct R to their location(s).


```r
chromosomeI <- "genomes/chromosome_I.fasta"
chromosomeII <- "genomes/chromosome_II.fasta"
```

You will also need to create a tabulated `sampleFile.txt` which should contain two columns `FileName` and `SampleName`, list the names of the files (exactly as they are) and a name for each sample.

If you multiple replicates you can either name them as `replicate_**` or if you give them all the same SampleName, R will know to treat them as replicates, but you will need to decide if that is how you want R to behave. If you are unsure stick with the former approach.

| FileName  | SampleName  |
|-----------|-------------|
| sample1.fastq  | sample1  |
| sample2.fastq  | sample2  |

Then assign this to `sampleFile` in R:


```r
sampleFile <- "raw_fastq/sampleFile.txt"
```

## Making the alignment(s)

If you have multiple chromosomes you can assign the alignment to each as a separate QuasR project.

Note:
You could choose to align using a BSgenome object instead, see the [QuasR documentation](https://bioconductor.org/packages/release/bioc/vignettes/QuasR/inst/doc/QuasR.html#617_using_a_bsgenome_package_as_reference_genome) if you prefer to do that.

Before you run the `qAlign()` function, make sure that you have created the folders to store your alignments (or you will get an error saying the directory does not exist).

To see a full list of other alignment parameters, use `?qAlign` to see the helpfile.


```r
proj_chrI <- qAlign(sampleFile, genome=chromosomeI, projectName = "chrI", alignmentsDir= "alignments/chrI/")
proj_chrII <- qAlign(sampleFile, genome=chromosomeII, projectName = "chrII", alignmentsDir= "alignments/chrII/")
```

If you run the same code again and you already have the output `.bam` and `.txt` files it will not repeat the alignment. If you alter any alignment parameters it will do a new alignment.

You can type the project name into the console to see which output files correspond to each alignment.


```r
proj_chrI
proj_chrII
```

### Auxiliary alignments

You can specify auxiliary genomes to be used to align unmapped sequences from the core genome. For example, if you have a spiked experiment or to check for contaminating sequences.

To do this supply a text file with a list of additional genome files like this:

| FileName  | AuxName |
|-----------|---------|
| NC_001422.1.fa  | phiX174 |

Then assign it in R:


```r
additionalGenomes <- "genomes/additionalGenomes.txt"
```

And then when running qAlign, add `auxiliaryFile = additionalGenomes` to the arguments list.


```r
proj_chrI <- qAlign(sampleFile, genome=chromosomeI, projectName = "chrI", alignmentsDir= "alignments/chrI/",
                    auxiliaryFile = additionalGenomes)
```


## Alignment stats and quality reports

You can use the `alignmentStats()` function to find out the number of mapped/unmapped reads.


```r
alignmentStats(proj_chrI)
alignmentStats(proj_chrII)
```

You can produce PDF quality reports for each alignment using the `qQCReport()` function.


```r
qQCReport(proj_chrI, pdfFilename = "quality_reports/chrI_quality.pdf")
qQCReport(proj_chrII, pdfFilename = "quality_reports/chrII_quality.pdf")
```
