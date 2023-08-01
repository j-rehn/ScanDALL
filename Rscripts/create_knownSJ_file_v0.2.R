#!/apps/bin/env Rscript

#############################################
#                                           #
#   Generate intron.bed files for use       #
#         with bedtools intersect           #
#                                           #
#############################################

# execute with Rscript --vanilla create_intron_bed.R -r gtf -e exclude.transcriptID.txt -o output.bed"

# load libraries
suppressMessages(library(tidyverse))
suppressMessages(library(rtracklayer))
suppressMessages(library(optparse))

# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

#making options
option_list = list(
  make_option(c("-r", "--ref"), type="character", default=NULL,
              help="gtf file", metavar="character"),
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="text file containing genes for inclusion in bed file"),
  make_option(c("-e","--exclude"), type="character", default=NULL,
              help="text file containing transcript.ID for exclusion", metavar="character"),
  make_option(c("-o","--output"), type="character", default=NULL,
              help="output bed file", metavar="character")

);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# error checking
if (is.null(opt$ref)){
  print_help(opt_parser)
  stop("gtf file must be supplied", sep="\n", call.=FALSE)
}

if (is.null(opt$output)){
  cat("Output file not supplied. File will be written to known.SJ.txt", sep="\n")
}

if (is.null(opt$input)){
  cat("Genes of interest not supplied. Splice junctions for all genes in gtf will be included.")
}

if (is.null(opt$exclude)){
  cat("transcriptID's for exclusion not supplied. Splice junctions for all transcripts in gtf will be included in known.SJ.txt file.", sep="\n")
}


# import UCSC_hg19 gtf file
cat("Importing gtf file", sep="\n")
gtf_obj <- import(opt$ref)

# import transcriptID's to exclude from known SJ file
if (!is.null(opt$exclude)) {

  # import list of transcriptID's to exclude from known.SJ file
  exclude <- read_delim(opt$exclude, delim = "\t", col_names = "transcript_id", show_col_types = FALSE)
  cat("Excluding SJ's for the following transcripts:", sep="\n")
  write.table(exclude, col.names = FALSE, row.names = FALSE, quote = FALSE)

} else {

  # create empty data.frame
  exclude <- data.frame(transcript_id = c(""))

}

# import genes's to include
if (!is.null(opt$input)) {
  
  # import list of genes's to include
  include <- read_delim(opt$input, delim = "\t", col_names = "gene_name", show_col_types = FALSE)
  cat("Including SJ's for the following genes:", sep="\n")
  write.table(include, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
} else {
  
  # create empty data.frame
  include <- data.frame(gene_name = c(""))

}

# specify output file name
if (is.null(opt$output)) {

  # not output file specified. Provide generic name.
  output_file <- "known.SJ.txt"

} else {

  # specify output file location/name as indicated
  output_file <- opt$output
}


###### Determine start and end positions of exons (represent the SJ coordinates) #######

# extract all data about gene of interest
data <- data.frame(gtf_obj[gtf_obj$gene_name %in% include$gene_name & gtf_obj$type == "exon", ])

cat("Determining SJ coordinates based on known transcripts.", sep="\n")

# start with positive strand genes
plusStrand_introns <- data %>%
  # select for positive strand genes
  filter(strand == "+") %>%
  # exclude transcripts in the exclude list (e.g. IKZF1 transcripts representing the exon 4-7 deletion)
  filter(!transcript_id %in% exclude$transcript_id) %>%
  # select required columns
  dplyr::select(seqnames,start,end,width,strand,
                gene_name,transcript_id,exon_number) %>%
  # group_by transcript_id to extract intron locations for each transcript separately
  group_by(transcript_id) %>%
  # calculate start and end position of intron for genes on positive strand
  mutate(intron_start = end + 1,
         intron_end = lead(start) - 1) %>%
  # drop rows without an intron_end
  filter(!is.na(intron_end)) %>%
  ungroup()


# extract minus strand intron information
minusStrand_introns <- data %>%
  # select for negative strand
  filter(strand == "-") %>%
  # exclude transcripts in the exclude list (e.g. IKZF1 transcripts representing the exon 4-7 deletion)
  filter(!transcript_id %in% exclude$transcript_id) %>%
  # select required columns
  dplyr::select(seqnames,start,end,width,strand,
                gene_name,transcript_id,exon_number) %>%
  # group_by transcript_id to extract intron locations for each transcript separately
  group_by(transcript_id) %>%
  # calculate start and end position of intron for genes on negative strand
  arrange(seqnames,transcript_id,start) %>%
  mutate(intron_start = end + 1,
         intron_end = lead(start)-1) %>%
  # drop rows without an intron_end
  filter(!is.na(intron_end)) %>%
  ungroup() %>%
  arrange(seqnames,transcript_id,desc(start))

cat("Writing known SJ to file.", sep="\n")

# Combine positive and negative strand intron information and sort the bed files
plusStrand_introns %>%
  bind_rows(minusStrand_introns) %>% 
  arrange(seqnames,intron_start,intron_end) %>%
  # reformat SJ coordinates for use with ScanDel
  mutate(SJ = paste(seqnames,intron_start,intron_end, sep = ":")) %>%
  # remove duplicate splicing events
  distinct(SJ) %>%
  # write to file
  write_delim(file = output_file,
              delim = "\t",
              col_names = FALSE)



