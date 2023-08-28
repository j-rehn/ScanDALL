# ScanDALL
Scan RNA-seq data for Deletions in Acute Lymphoblastic Leukaemia (ALL)

## Overview
Jacqueline Rehn : jacqueline.rehn@sahmri.com  Jimmy Breen : Jimmy.Breen@telethonkids.org.com

This algorithm uses splice-junction information from STAR-aligned RNA-seq data to predict the presence of common gene deletions observed in acute lymphoblastic leukaemia (ALL). 

## Instructions for use
ScanDALL requires BEDTools to run, specifically the bedtools intersect program from this suite of tools. Please ensure this is installed and available in your path before running ScanDALL. ScanDALL also requires a BED file, containing exon locations for the canonical transcript for the gene(s) of interest, and a text file containing known splice junction information for annotated transcripts from the relevant reference annotation. R scripts are provided to generate your own exon.bed and annotated.SJ files for genes of interest. Alternatively, files have been pre-generated for 8 ALL associated genes (IKZF1, PAX5, ETV6, RB1, ERG, LEF1, PTEN and STIL) using ensembl GRCh37 and UCSC hg19 reference annotations.

For installation of ScanDALL, clone this github repository:

```
git clone https://github.com/j-rehn/ScanALL.git
```

## Quick Start
