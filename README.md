# ScanDALL
Scan RNA-seq data for Deletions in Acute Lymphoblastic Leukaemia (ALL)

## Overview
Jacqueline Rehn : jacqueline.rehn@sahmri.com  Jimmy Breen : Jimmy.Breen@telethonkids.org.com

This algorithm uses splice-junction information from STAR-aligned RNA-seq data to predict the presence of common gene deletions observed in acute lymphoblastic leukaemia (ALL). 

## Instructions for use
ScanDALL requires BEDTools to run, specifically the bedtools intersect program from this suite of tools. Please ensure this is installed and available in your path before running ScanDALL. ScanDALL also requires a BED file, containing exon locations for the canonical transcript for the gene(s) of interest, and a text file containing known splice junction information for annotated transcripts from the relevant reference annotation. R scripts are provided to generate your own exon.bed and annotated.SJ files for genes of interest. Alternatively, files have been pre-generated for 8 ALL associated genes (IKZF1, PAX5, ETV6, RB1, ERG, LEF1, PTEN and STIL) using ensembl GRCh37 and UCSC hg19 reference annotations.

For installation of ScanDALL, clone this github repository and change into the created directory:

```
git clone https://github.com/j-rehn/ScanALL.git
cd ScanDALL
```

## Quick Start
The script scandall.sh requires the SJ.out.tab file generated by STAR during alignment for the sample of interest. Additionally, the user must supply the path to the relevant BED file and annotated SJ text file. Optionally, users may select to use RRS or PRS as the normalisation score. If not supplied, RRS is used as default.

1. To view the required arguments for scandall.sh:
```
bash scandall.sh

Usage: scanDALL.sh [SJ.out.tab] [exon.bed] [anno.SJ.txt] [norm.score]
Incorrect number of arguments
- Required: SJ.out.tab = SJ.out.tab file
- Required: exon.bed = bed file containing exon locations for genes of interest
- Required: anno.SJ.tx = text file containing known SJ locations for genes of interest
- Optional: norm.score = normalisation score to use for filtering. Options [RRS] (default) or [PRS]
```
Full path to each of these files needs to be provided.

2. To run ScanDALL for detetion of ALL associated gene deletions on your own sample (e.g. mysample.SJ.out.tab) using the UCSC hg19 reference annotation:
```
bash scandall.sh path/to/mysample.SJ.out.tab bed_files/hg19.keyGenes.exon.bed bed_files/known.SJ.hg19.txt
```
Multiple temporary files will be generated at the location where the script is run. These are removed once ScanDALL is complete. There are 2 final output files: 
  a. mysample.exonDel.bed containing information about splice-junctions (SJs) which spanned multiple exons for the genes of interest. Many of these may not have passed subsequent filtering steps. This information includes chr:position of SJ start, chr:position SJ end, gene, strand, exons the SJ spans, number of reads supportin the SJ. E.g.
```
12:12022904	12:12043874	ETV6	+	6-7	12
21:39763638	21:39774478	ERG	-	8-9	16
21:39764367	21:39775427	ERG	-	7-8	13
4:108969908	4:108991818	LEF1	-	10-11	7
4:108969908	4:109004511	LEF1	-	6-11	22
4:108991927	4:109088710	LEF1	-	2-8	8
4:109010414	4:109088710	LEF1	-	2-3	6
7:50343956	7:50444230	IKZF1	+	1-3	6
7:50358698	7:50467615	IKZF1	+	3-7	8
7:50367354	7:50467615	IKZF1	+	4-7	628
7:50367354	7:50529860	IKZF1	+	4-8	16
```
  
  b. mysample.exonDel.anno containing information on predicted deletions which passed filters. Data includes Chr, Gene, Strand, exons deleted, number of distinct splice junctions supporting this deletion, number of reads supporting the SJ, total number of reads mapping to the acceptor and donor sites of the deletion SJ, PRS value, median number of reads mapping to the annotated SJ's in the deletion window, RRS value, confidence level, SJ start - SJ end positions. E.g.
 E.g.
```
7	IKZF1	+	4-7	1	628	1302	0.482335	158	0.798982	high	50367354-50467615
```  

## Application of ScanDALL to alternate genes or reference genomes

