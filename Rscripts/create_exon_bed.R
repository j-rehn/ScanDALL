#!/apps/bin/env Rscript

# execute with Rscript --vanilla create_exon_bed.R -r gtf -i transcriptID.txt -o output.bed"

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
              help="text file containing transcriptID's for inclusion in bed file"),
  make_option(c("-o","--output"), type="character", default=NULL,
              help="output bed file", metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# error checking
if (is.null(opt$ref)){
  print_help(opt_parser)
  stop("gtf file must be supplied", call.=FALSE)
}

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("trascriptID's for genes to include in bed file must be supplied", call.=FALSE)
}

if (is.null(opt$output)){
  cat("Output file not supplied. Bed file will be written to canonical.exon.bed", call.=FALSE)
}


# import UCSC_hg19 gtf file
gtf_obj <- import(opt$ref)

# specify output file name
if (is.null(opt$output)) {
  
  # not output file specified. Provide generic name.
  output_file <- "canonical.exon.bed"
  
} else {
  
  # specify output file location/name as indicated
  output_file <- opt$output
}


########## Generate exon bed containing transcripts of interest ##########

# import transcriptIDs for transcripts to include in bed file
transcripts <- read_delim(opt$input, delim = "\t", col_names = FALSE)
  
# extract exon locations in bed file format for for indicated transcript
data.frame(gtf_obj[gtf_obj$type == "exon", ]) %>% 
  # select for key transcripts
  filter(transcript_id %in% transcripts$X1) %>% 
  # calculate total number of features (exons) for each transcript (gene)
  add_count(transcript_id) %>% 
  # reformat as bed file
  dplyr::select(seqnames,start,end,gene_name,strand,n) %>% 
  ### coordinates for start position of all transcripts off-by-one (GTF 1-based but bed 0-based for start position)
  # adjust start coordinates to generate 0-based position
  mutate(start = (start - 1)) %>% 
  # add exon_numbers (different order depending on whether gene appears on '+' or '-' strand)
  group_by(gene_name) %>% 
  mutate(counter = row_number()) %>% 
  arrange(desc(start)) %>% 
  mutate(reverse_counter = row_number()) %>% 
  arrange(seqnames, start) %>% 
  # select the counter that reflects correct exon order for strand
  mutate(exon_number = ifelse(strand == "+", counter, reverse_counter)) %>% 
  dplyr::select(-counter,-reverse_counter) %>% 
  # write this to a tab separated file without header
  write_delim(file = output_file,
              delim = "\t", col_names = FALSE)



