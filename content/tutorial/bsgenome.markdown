+++
title = "Forging a BSgenome package"

date = 2018-09-09T00:00:00
# lastmod = 2018-09-09T00:00:00

draft = false  # Is this a draft? true/false
toc = true  # Show table of contents? true/false
type = "docs"  # Do not modify.

# Add menu entry to sidebar.
linktitle = "Forging a BSgenome package"
[menu.tutorial]
  parent = "ChIP-seq with R"
  weight = 3
+++

## Introduction

In this tutorial we will use the BSgenome package to create an R package that contains the genome for _Vibrio cholerae_, but you can replace it with __your favourite organism__. 

The [BSgenome](https://bioconductor.org/packages/release/bioc/html/BSgenome.html) package provides a framework for interacting with genome information in R.

Alert:
There may already be a BSgenome package for your favourite organism. Check the [list](https://kasperdanielhansen.github.io/genbioconductor/html/BSgenome.html) of available genomes first.

Once you have _forged_ (created) and installed the package, you will be able to load the genome as you would any other R package, like this:


```r
library(BSgenome.Vcholerae.NCBI.N16961)
```

If _V. cholerae_ happens to be __your favourite organism__ too, and want to save yourself some time, you can [download](/tutorial/files/BSgenome.Vcholerae.NCBI.N16961_1.0.0.tar.gz) the package I created, and skip to the 'Install your genome package' section to get started.

## Install the BSgenome package


```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("BSgenome", version = "3.8")
```

## Download the genome sequence

You will need to download the fasta files for the genome you want to forge. Since _V. cholerae_ has two chromosomes, these are downloaded separately. You can use whatever source you like (EMBO, NCBI etc..) but make sure you get it in __fasta__ file type.  

For the _V.cholerae_ reference genome (N16961):  

 - [chromosome I](https://www.ncbi.nlm.nih.gov/nuccore/NC_002505.1)  
 - [chromesome II](https://www.ncbi.nlm.nih.gov/nuccore/NC_002506.1)   

At the NCBI website, click send to file, then select FASTA.    

## Prepare your files

I created a new folder on my Desktop and set this as my working directory:


```r
setwd("~/Desktop/genomepackage")
```

Within this create another named `seqs_srcdir`, move the fasta files here.  

Make sure the files are appropriately named (i.e.. as-is from source) for _V. cholerae_ the two chromosome fasta files are named: `NC_002505.1.fa` and  `NC_002506.1.fa`. It is tempting to name them something more readable like "chromosome 1" but this can cause problems later.

Ensure the file extensions are `.fa` if they are not already, on Mac double check by right-clicking and choosing 'Get Info' because it could still be `.fasta` - if so, amend it.

### Seed file

The seed file contains all the relevant metadata for the BSgenome package, so it is worth supplying as much information as you can. The easiest way to make a seed file is to edit one that already exists, so you can download my [seed file](/tutorial/files/BSgenome.Vcholerae.NCBI.N16961-seed) and use it as a template.

You will need to use a text editor such as TextEdit on Mac (right-click and select open with > TextEdit) or RStudio. Use NCBI to populate the relevant information:

>Package: BSgenome.Vcholerae.NCBI.N16961  
Title: Full genome sequence for Vibrio cholerae O1 biovar El Tor str N16961  
Description: Full genome sequence for the two chromosomes of Vibrio cholerae El Tor N16961 provided by NCBI  
Version: 1.0.0  
organism: Vibrio cholerae  
common_name: V. cholerae  
provider: NCBI  
provider_version: ASM674v1  
release_date: 2014/02  
release_name: N16961  
source_url: https://www.ncbi.nlm.nih.gov/genome/?term=Vibrio%20cholerae  
organism_biocview: Vibrio_cholerae  
BSgenomeObjname: Vcholerae  
seqnames: c("NC_002505.1","NC_002506.1")  
seqs_srcdir: /User/Desktop/genomepackage/seqs_srcdir  

For genomes with multiple chromosomes, list them as a vector (see above example). Otherwise `seqnames: chromosomenameFileName` is sufficient (you can remove the c()). 

The `BSgenomeObjname` is important because this is the name you will use to access the package in R once it has been installed.  

Save it as is, and then edit the file name. Be careful with the file extensions, double check using 'Get Info' to ensure it has not been changed to `.txt` or anything else. 

## Forge the package

The package is forged using the `forgeBSgenomeDataPkg` function.  

Simply use the name of the seed file as the only argument and it will create your package files to the same directory.  

Alert:  
Double check the sequence files are `.fa` file types and that the details in the seed are correct before running.


```r
forgeBSgenomeDataPkg("BSgenome.Vcholerae.NCBI.N16961-seed")
```

If you need to run the function again, delete the previous package files first. 

## Install your genome package

To install the genome package you will need to use the Mac command line (Terminal). 

- Close R 
- Open Terminal
- In Terminal navigate to your working directory:
    + use `ls` to see list of files in the current directory
    + use `cd` to move to a directory (i.e.. `cd Desktop`)
- Run `R CMD build BSgenome.Vcholerae.EBI.N16961` to compile the package
- Run `R CMD check BSgenome.Vcholerae.EBI.N16961.tar.gz` to check it
- Run `R CMD INSTALL BSgenome.Vcholerae.NCBI.N16961_1.0.0.tar.gz`

Alert:  
If you have downloaded my _V. cholerae N16961_ [package](/tutorial/files/BSgenome.Vcholerae.NCBI.N16961_1.0.0.tar.gz) you will need to navigate to wherever you have saved the file.  Then run:  
`R CMD INSTALL BSgenome.Vcholerae.NCBI.N16961_1.0.0.tar.gz`

And you're done! It should now be ready to use.  

## Loading and accessing the genome in R

To use the genome in R you will need to load the package using the `library()` function.


```r
library(BSgenome.Vcholerae.NCBI.N16961)
```

Enter the `BSgenomeObjname` (in this case Vcholerae) to print some general information about the genome to the console.


```r
Vcholerae
```

```
## V. cholerae genome:
## # organism: Vibrio cholerae (V. cholerae)
## # genome: ASM674v1
## # provider: NCBI
## # release date: 2014/02
## # 2 sequences:
## #   NC_002505.1 NC_002506.1                                                
## # (use 'seqnames()' to see all the sequence names, use the '$' or '[[' operator
## # to access a given sequence)
```

`length()` tells you how many chromosomes there are


```r
length(Vcholerae)
```

```
## [1] 2
```

`Vcholerae$NC_002505.1` tells you the length of chromosome I and a bit of its sequence


```r
Vcholerae$NC_002505.1
```

```
## 2961149-letter DNAString object
## seq: AGGGTCATTAAATATATATAAAGATCTATATAGAGA...GGCTAGAAAATCGCTTTCCTGTTTTTTCGATCAAGG
```

`alphabetFrequency(Vcholerae$NC_002505.1)` shows you the ACGT content of chromosome I


```r
alphabetFrequency(Vcholerae$NC_002505.1)
```

```
##      A      C      G      T      M      R      W      S      Y      K      V 
## 769234 703384 708931 779567      0      0      0      0      0      0      0 
##      H      D      B      N      -      +      . 
##      0      0      0     33      0      0      0
```

You can extract sequence from specific co-ordinates, for example to select the sequence from position 45 to 65 on chromosome I:


```r
Vcholerae$NC_002505.1[45:65]
```

```
## 21-letter DNAString object
## seq: TTAGATCTACTATTAAGGAGC
```

You could then store this in a data frame, export the data frame as fasta file (with multiple sequence etc..) and use it somewhere else...
