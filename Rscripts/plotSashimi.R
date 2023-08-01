
#########################
#########################
###                   ###
###    plotSashimi    ###
###                   ###
#########################
#########################

## Create a sashimi plot from SJ.out.tab file for a specified gene/region

# load libraries
library(tidyverse)
library(rtracklayer)

## required parameters

# gtf : imported gtf_obj
# SJfile : location of SJ.out.tab for sample of interest
# gene2plot : gene wish to create sashimi plot for
# transcript2plot : transcript_id for gene transcript to use for plotting

## optional parameters

# excludeTranscripts : vector of transcript id to exclude before annotating splice junctions
# minReads : minimum uniq reads supporting splice junction
# SJ_anno : types of splice junctions to include in plot
# propCov :

plotSashimi <- function(gtf, ref, SJfile, gene2plot, transcript2plot,
                        excludeTranscripts = NULL,
                        minReads = 5,
                        SJ_anno = NULL, ...){

  # extract all data about gene of interest
  data <- data.frame(gtf[gtf$gene_name == gene2plot])

  # remove transcripts annotated as exclude
  if(!is.null(excludeTranscripts)){
    data <- filter(data, !transcript_id %in% excludeTranscripts)
  }

  ## subset data for basic gene information:
  ## requires different code depending of gtf file
  if(ref == "ensembl"){

    gene_info <- filter(data, type == "gene") %>%
      dplyr::select(seqnames,start,end,width,strand,source,gene_id,
                    gene_name,gene_version,gene_biotype)

  }else if(ref == "UCSC"){

    gene_info <- filter(data, type == "transcript" & transcript_id == transcript2plot) %>%
      dplyr::select(seqnames,start,end,width,strand,source,gene_name) %>%
      slice_max(start) %>%
      slice_max(end)

  }


  # extract exon data for transcript of interest
  exon_info <- filter(data, type == "exon", transcript_id == transcript2plot) %>%
    dplyr::select(seqnames,start,end,width,exon_number,strand,source,
                  transcript_id,
                  # transcript_version,
                  gene_name)

  ## exons in UCSC gtf file are labelled according to genomic location regardless of strand.
  ## need to reverse the order of exons for genes on minus strand if using UCSC

  if(gene_info$strand == "-"){

    exon_info <- exon_info %>%
      dplyr::select(-exon_number) %>%
      arrange(transcript_id, desc(start)) %>%
      mutate(exon_number = row_number())

  }

  # import SJ.out.tab and extract SJ data for gene of interest
  SJ_out <- read_delim(file = SJfile,
                       col_names = c("seqnames","SJ_start","SJ_end","strand",
                                     "intron_motif","annotation","uniq_reads",
                                     "multimap_reads","max_overhang"),
                       col_types = cols(seqnames = col_character(),
                                        SJ_start = col_integer(),
                                        SJ_end = col_integer(),
                                        strand = col_integer(),
                                        intron_motif = col_integer(),
                                        annotation = col_integer(),
                                        uniq_reads = col_integer(),
                                        multimap_reads = col_integer(),
                                        max_overhang = col_integer())) %>%
    # subset for region involving gene of interest
    filter(seqnames %in% gene_info$seqnames,
           (SJ_start >= gene_info$start & SJ_start <= gene_info$end) | (SJ_end >= gene_info$start & SJ_end <= gene_info$end)) %>%
           # SJ_start >= gene_info$start,
           # SJ_end <= gene_info$end) %>%
    # filter to exclude SJ with more multimap_reads than uniq_reads & fewer than 5 uniq_reads
    filter(uniq_reads > multimap_reads,
           uniq_reads >= minReads) %>%
    arrange(SJ_start) %>%
    dplyr::select(SJ_start,SJ_end,uniq_reads)

  ### Scale transcript for plotting ###

  # first convert start and end location to coordinate positions beginning at 1
  exon_plot_data <- exon_info %>%
    arrange(start) %>%
    # convert loc to coordinate beginning at 0
    mutate(plot_start = start - start[1]) %>%
    mutate(plot_end = end - start[1]) %>%
    ## scale exons & introns to enable plotting within a specified range (0-10)
    # determine scaling factor
    mutate(scale_factor = 10/max(plot_end)) %>%
    # apply scaling factor to each start and stop position of exons
    mutate(plot_start = plot_start*scale_factor,
           plot_end = plot_end*scale_factor)

  ## reverse coordinates for plotting if gene is on the minus strand
  if(gene_info$strand == "-"){
    exon_plot_data <- exon_plot_data %>%
      mutate(plot_start = (plot_start*-1)+10,
             plot_end = (plot_end*-1)+10) %>%
      # switch plot start and plot end
      dplyr::rename(plot_start = plot_end,
                    plot_end = plot_start)
  }

  # apply same scaling factor to SJ data
  SJ_out <- SJ_out %>%
    # adjust coordinate position
    mutate(plot_SJ_start = SJ_start - exon_plot_data$start[1]) %>%
    mutate(plot_SJ_end = SJ_end - exon_plot_data$start[1]) %>%
    # apply scaling factor
    mutate(plot_SJ_start = plot_SJ_start*(exon_plot_data$scale_factor[1]),
           plot_SJ_end = plot_SJ_end*(exon_plot_data$scale_factor[1]))

  ### annotate detected SJ as being known or novel ###

  # extract recognised SJ from gtf file
  known_introns <- data %>%
    filter(type == "exon") %>%
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
    # calculate total number of introns per transcript + intron number
    arrange(seqnames,transcript_id,desc(start)) %>%
    mutate(intron_number = row_number()) %>%
    add_count(transcript_id) %>%
    # reformat as bed file
    arrange(seqnames,transcript_id,start) %>%
    ungroup() %>%
    # reformat as bed file
    dplyr::select(intron_start,intron_end,
                  transcript_id,intron_number) %>%
    # exclude duplicate introns from separate transcripts
    distinct(intron_start,intron_end,.keep_all = TRUE)

  ## add annotation info to SJ_out and use to specify colours for plotting SJ

  # must account for strand
  if(gene_info$strand == "+"){

    SJ_out <- SJ_out %>%
      left_join(known_introns, by = c("SJ_start" = "intron_start", "SJ_end" = "intron_end")) %>%
      # annotate as known donor, known acceptor, known donor & acceptor or novel donor & acceptor
      mutate(anno = ifelse(SJ_start %in% known_introns$intron_start & SJ_end %in% known_introns$intron_end & is.na(transcript_id), "NDA",
                           ifelse(SJ_start %in% known_introns$intron_start & !SJ_end %in% known_introns$intron_end & is.na(transcript_id), "D",
                                  ifelse(!SJ_start %in% known_introns$intron_start & SJ_end %in% known_introns$intron_end & is.na(transcript_id), "A",
                                         ifelse(!SJ_start %in% known_introns$intron_start & !SJ_end %in% known_introns$intron_end & is.na(transcript_id), "N","DA"))))) %>%
      # also adjust the thickness of the curved line
      mutate(prop = uniq_reads/sum(uniq_reads)) %>%
      mutate(line = cut(prop, breaks = 3, labels = c("low","med","high")))

  } else if(gene_info$strand == "-") {

    SJ_out <- SJ_out %>%
      left_join(known_introns, by = c("SJ_start" = "intron_start", "SJ_end" = "intron_end")) %>%
      # annotate as known donor, known acceptor, known donor & acceptor or novel donor & acceptor
      mutate(anno = ifelse(SJ_start %in% known_introns$intron_start & SJ_end %in% known_introns$intron_end & is.na(transcript_id), "NDA",
                           ifelse(SJ_start %in% known_introns$intron_start & !SJ_end %in% known_introns$intron_end & is.na(transcript_id), "A",
                                  ifelse(!SJ_start %in% known_introns$intron_start & SJ_end %in% known_introns$intron_end & is.na(transcript_id), "D",
                                         ifelse(!SJ_start %in% known_introns$intron_start & !SJ_end %in% known_introns$intron_end & is.na(transcript_id), "N","DA"))))) %>%
      # also adjust the thickness of the curved line
      mutate(prop = uniq_reads/sum(uniq_reads)) %>%
      mutate(line = cut(prop, breaks = 3, labels = c("low","med","high"))) %>%
      # and reverse coordinates for plotting
      mutate(plot_SJ_start = (plot_SJ_start*-1)+10,
             plot_SJ_end = (plot_SJ_end*-1)+10) %>%
      # switch plot start and plot end
      dplyr::rename(plot_SJ_start = plot_SJ_end,
                    plot_SJ_end = plot_SJ_start)

  }

  # plot only requested SJ annotations
  if(!is.null(SJ_anno)){
    SJ_out <- filter(SJ_out, anno %in% SJ_anno)
  }


  # create plot annotation (gene name, transcript_id, scale)
  plot_anno <- exon_plot_data %>%
    distinct(scale_factor,gene_name,strand,transcript_id) %>%
    mutate(scale_seg_length = 10000*scale_factor) %>%
    # determine position
    mutate(line_end = 10, line_start = 10-scale_seg_length)

  # determine plot breaks on x-axis
  f <- seq(floor(min(SJ_out$plot_SJ_start)),ceiling(max(SJ_out$plot_SJ_end)))

  # create base plot
  p <- exon_plot_data %>%
    ggplot() +
    # theme_bw() +
    theme_void() + 
    geom_segment(x=0, y=1, xend=10, yend=1, colour = "grey50") +
    geom_rect(aes(xmin = plot_start, xmax = plot_end, ymin = 0.75, ymax = 1.25),
              fill = "grey70", colour = "grey50") +
    # # number exons on the plot
    # geom_text(aes(x = plot_start + ((plot_end - plot_start)/2), y = 0.7, label = exon_number),
    #           check_overlap = TRUE) +
    geom_segment(aes(x = line_start, xend = 10, y = 0.5, yend = 0.5),
                 data = plot_anno, size = 1, lineend = "round") +
    geom_text(aes(x = line_start + ((line_end - line_start)/2), y = 0.45, label = "10 Kb"),
              data = plot_anno) +
    # add to plot gene name and strand (changing to white to keep plot dimension but moving labelbelow exons)
    geom_text(aes(x = 0, y = 2, label = paste(gene_name, " ( ",strand," strand )",sep = "")),
              data = plot_anno, hjust = "left", size = 5, colour = "white") +
    # add to plot transcript ID used for exons
    geom_text(aes(x = 0, y = 1.9, label = transcript_id),
              data = plot_anno, hjust = "left", size = 5, colour = "white") +
    ## re-add gene name & transcripitID
    geom_text(aes(x = 0, y = 0.65, label = gene_name),
              data = plot_anno, hjust = "left", size = 7, fontface = "italic") +
    geom_text(aes(x = 0, y = 0.53, label = transcript_id),
              data = plot_anno, hjust = "left", size = 5.5) +
    scale_color_manual(name = "Splice junction annotation",
                       values = c("DA" = "#377EB8",
                                  "NDA" = "#E41A1C",
                                  "D" = "#4DAF4A",
                                  "A" = "#A65628",
                                  "N" = "#984EA3"),
                       labels = c("Known donor-acceptor combo (DA)",
                                  "Novel donor-acceptor combo (NDA)",
                                  "Known donor, novel acceptor (D)",
                                  "Known acceptor, novel donor (A)",
                                  "Novel donor and novel acceptor (N)")) +
    scale_x_continuous(breaks = f) +
    # scale_y_continuous(limits = c(0,2)) +
    theme(legend.position = "none",
          # legend.direction="vertical",
          axis.title = element_blank(),
          # axis.text = element_text(size = 12),
          # legend.title = element_text(size = 12, face = "bold"),
          # legend.text = element_text(size = 12)
          axis.text = element_blank()
          ) +
    guides(color=guide_legend(ncol=3))

  # add SJ's to plot as separate layers
  if(nrow(SJ_out[SJ_out$line == "high",]) > 0) {
    p <- p + geom_curve(aes(x = plot_SJ_start, xend = plot_SJ_end, y = 1, yend = 1, colour = anno),
                        data = SJ_out[SJ_out$line == "high",], curvature = -0.5, size = 2)
  }
  if(nrow(SJ_out[SJ_out$line == "med",]) > 0){
    p <- p+ geom_curve(aes(x = plot_SJ_start, xend = plot_SJ_end, y = 1, yend = 1, colour = anno),
                       data = SJ_out[SJ_out$line == "med",], curvature = -0.5, size = 1)
  }
  if(nrow(SJ_out[SJ_out$line == "low",]) > 0){
    p <- p + geom_curve(aes(x = plot_SJ_start, xend = plot_SJ_end, y = 1, yend = 1, colour = anno),
                        data = SJ_out[SJ_out$line == "low",], curvature = -0.5, size = 0.5)
  }

  # return final plot
  return(p)

}



