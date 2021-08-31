+++
title = "Chromosome overview plots"

date = 2020-01-27T00:00:00
# lastmod = 2019-04-19T00:00:00

draft = false  # Is this a draft? true/false
toc = true  # Show table of contents? true/false
type = "docs"  # Do not modify.
+++



## Introduction

In this tutorial we will use the [ggbio](https://bioconductor.org/packages/release/bioc/html/ggbio.html) package to build overview plots. 

These are good for showing a simple representation of ChIP-seq and other datasets.

![Example of a simple overview plot. In this plot genes of the positive (red) and negative (blue) strands of a single bacterial chromosome are shown as rectangles on two tracks.](/tutorial/circle_files/example.png)

In the example above only two tracks are shown, which are a representation of the all the genes encoded on the chromosome. Inner tracks can then be added to represent over types of data that have been aligned to the genome, such as sequencing reads.

## Required packages and data

Install (if necessary) and then load the following packages:

- tidyverse  
- data.table  
  
These packages are bioconductor packages:  

- ggbio  
- AnnotationDbi  
- GenomicRanges  
- Biostrings  
- GenomicFeatures  

You will also need to import annotation data for the chromosome(s) as a TxDB. You can use GFF to make a TxDB or download a TxDB that already exists for the organism.

### Create a TxDB object

A TxDB object contains genome annotation information. It is very simple to make, once you have done you can save and reload it whenever you wish to use it. 

- Download a GFF or GTF annotation file (from NCBI)


```r
library(GenomicFeatures)
library(AnnotationDbi)
<genomeName>TxDB <- makeTxDbFromGFF("annotation<genomeName>.gff")
saveDb(<genomeName>TxDB, file="<genomeName>TxDB.sqlite")
```

- Loading the TxDB file (replacing <genome>, with a name for the genome you are using):


```r
<genome>Annotation <- loadDb("<genomeName>TxDB.sqlite")
```

## Creating a gene track

This track shows the genes on the positive and negative strand of the chromosome. To do this in a single track you need to collapse it, because there are often overlapping genes in different reading frames.  


```r
p_strand <- as_tibble(cds(<genome>Annotation, filter = list(tx_strand = "+"))) %>%
  dplyr::select(c('seqnames', 'start', 'end')) %>%
  arrange(start) %>%
  dplyr::group_by(g = cumsum(cummax(lag(end, default = first(end))) < start)) %>% 
  dplyr::summarise(start = first(start), end = max(end)) %>%
  dplyr::select(c('start', 'end'))

p_strand$seqnames <- "NC_000000.0"

n_strand <- as_tibble(cds(chrIAnnotation, filter = list(tx_strand = "-"))) %>%
  dplyr::select(c('seqnames', 'start', 'end')) %>%
  dplyr::arrange(start) %>%
  dplyr::group_by(g = cumsum(cummax(lag(end, default = first(end))) < start)) %>% 
  dplyr::summarise(start = first(start), end = max(end)) %>%
  dplyr::select(c('start', 'end'))

n_strand$seqnames <- "NC_000000.0"
```

First, import the coding sequences from the annotation file and filter for the positive strand only.  

Then arrange data frame by the start position of each gene.   

The group_by and summarise functions are used to 'collapse' the data, by using the start and end positions of genes to merge overlapping coding regions in the different reading frames.  

`lag(end)` adds a column that tracks the previous end values, `cummax` stores the highest of previous end values which is then compared with the end. When a start value is higher than the maximum of the previous stop, it creates a new group. `cumsum` tracks the number of TRUEs and gives each different group (i.e.. overlapping genes form groups) a unique number. `summarise` aggregates based on the group id using the first start value of each group and the largest end value.  

Finally need to add back in the `seqname` column. This should be the RefSeq value from NCBI and should match up the chromosome name used in the reference genome. This is repeated for the negative strand.  

Next you create GRanges objects with the `p_strand` and `n_strand` information using the `makeGRangesFromDataFrame` function from the `GenomicRanges` package.


```r
pos_genes <- makeGRangesFromDataFrame(p_strand)
neg_genes <- makeGRangesFromDataFrame(n_strand)
```

## Creating the plot

The plot is created using the [ggbio](https://bioconductor.org/packages/release/bioc/html/ggbio.html) package. It works similarly to ggplot, so you can add multiple elements by using '+' and various geoms (such as point, line, link, etc..) For the gene track we will use the rectangle (`rect`) geom.


```r
plot <- ggbio(buffer = 0.2, radius = 10) +
  #circle(avg_reads_track, geom = "line", trackWidth= 3,
  #      color = "dodgerblue2", aes(y=average), buffer=35, radius=6, space.skip=0) +
  circle(neg_genes, geom = "rect", color=NA, 
        fill='blue', trackWidth=0.5, space.skip=0) +
  circle(pos_genes, geom = "rect", color=NA, 
         fill='red', trackWidth=0.5, space.skip=0) 
plot
```

The plot is initialised using a call to the `ggbio` function. The two parameters:  

- buffer: the space between each track  
- radius: the size of the circle plot (might need to increase this depending on the number of tracks you want to add)  

The `circle()` function is called for creating a track in the plot, there a number of parameters that can be set to alter the aesthetics of the track. You can view the help file and all the parameters that can be altered by entering `?circle` in the R console. 

![Help page for the ggbio circle layout function.](/tutorial/circle_files/circle_help.png)

The `space.skip` parameter set to zero means that there is no space between the beginning and end of the chromosome. This is required because the default behaviour is designed to have space between multiple chromosomes. 





