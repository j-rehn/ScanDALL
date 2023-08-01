#!/bin/bash

# help page
if [ "$#" != "3" ]; then
	echo "Usage: `basename $0` [SJ.out.tab] [exon.bed] [anno.SJ.txt] [norm.score]"
  	echo "Incorrect number of arguments"
	echo "- Required: SJ.out.tab = SJ.out.tab file"
	echo "- Required: exon.bed = bed file containing exon locations for genes of interest"
	echo "- Required: anno.SJ.tx = text file containing known SJ locations for genes of interest"
	echo "- Optional: norm.score = normalisation score to use for filtering. Options [RRS] (default) or [PRS]"
	exit 0
fi

#### Required: bedtools v2.26.0 or higher ###

# supplied variables
sample=`basename $1 .SJ.out.tab`
exon_bed=$2
knownSJ=$3
norm=${4:-RRS}

# start time
start_time=$SECONDS

#################################################
#################################################
###                                           ###
###   Step01: extract exon skipping events    ###
###                                           ###
#################################################
#################################################

echo -e "Running Deletion Detection on ${sample}"

# run bedtools intersect to extract SJ spanning 2 or more canonical exons in key genes
bedtools intersect -a $1 -b ${exon_bed} -F 1 -c | \
awk '{if($10 > 1 && $7 > 5 && $7 > $8){print $0}}' | \
bedtools intersect -a stdin -b ${exon_bed} -wo > ${sample}.unfilt.exonDel.bed

# perform check and exit program at this stage if no deletions detected
if [ ! -s ${sample}.unfilt.exonDel.bed ]
then
  echo "No deletions detected"
  rm ${sample}.unfilt.exonDel.bed
  touch ${sample}.exonDel.bed
  echo -e "chr\tgene\tstrand\tdeletion\tn_SJ\tsupport_reads\ttotal_reads\tPRS\tmedian_reads\tRRS\tconfidence\tSJ" > ${sample}.exonDel.anno
  exit
else
  echo -e "SJ's spanning multiple exons identified"
fi

#################################################
#################################################
###                                           ###
###     Step02: identify exons involved       ###
###                                           ###
#################################################
#################################################

# select required fields from exonDel.bed file
awk -vOFS="\t" '{print $14, $15, $1, $2, $3, $7, $17}' ${sample}.unfilt.exonDel.bed | \
# combine gene identifiers and SJ start and stop to generate a unique deletion identifier
awk -vOFS="\t" '{print $8=$1":"$2":"$3":"$4":"$5":"$6,$7}' | \
# array to select the max (and min) exon number thus represent the deleted exons
awk -vOFS="\t" '{Del[$1]++; min[$1]=Del[$1]==1||min[$1]>$2?$2:min[$1]; max[$1]=max[$1]<$2?$2:max[$1]} \
  END { for (var in Del) print var, min[var], max[var]}' | \
# split combined deletion identifier
awk -F '[:\t]' '{ print $1,$2,$3,$4,$5,$6,$7,$8}' OFS="\t" | \
# combine chr:SJstart:SJend into single identifier, followed by geneName, strand, uniq_map and deleted exons (e.g. 4-8)
awk -vOFS="\t" '{print $3":"$4":"$5,$1,$2,$6,$7"-"$8}' > ${sample}.tmp.exonDel.bed

#################################################
#################################################
###                                           ###
###       Step03: exclude known SJ's          ###
###                                           ###
#################################################
#################################################

echo -e "Excluding annotated SJ's"

# use known.SJ.txt to remove known SJ from tmp.exonDel.bed
awk 'NR == FNR {a[$1]; next} !($1 in a)' ${knownSJ} ${sample}.tmp.exonDel.bed | \
# edit SJ identifier to be formatted as chr:SJstart, chr:SJend
awk -F '[:\t]' '{ print $1":"$2,$1":"$3,$4,$5,$7,$6}' OFS="\t" | \
# sort the output by genomic location
sort -t$'\t' -k1,1 > ${sample}.tmp2.exonDel.bed

# perform check and exit program at this stage if no deletions detected
if [ ! -s ${sample}.tmp2.exonDel.bed ]
then
  echo "Identified SJ's are present in annotated transcripts. No deletions detected"
  rm ${sample}.*
  touch ${sample}.exonDel.bed
  echo -e "chr\tgene\tstrand\tdeletion\tn_SJ\tsupport_reads\ttotal_reads\tPRS\tmedian_reads\tRRS\tconfidence\tSJ" > ${sample}.exonDel.anno
  exit
else
  echo -e "Excluding SJ's with novel acceptor & donor site"
fi

#################################################
#################################################
###                                           ###
###       Step04: exclude novel SJ's          ###
###                                           ###
#################################################
#################################################

# generate tmp file containing known SJ acceptor and donor sites
cat <(awk -F '[:\t]' '{print $1":"$2}' OFS="\t" ${knownSJ}) <(awk -F '[:\t]' '{print $1":"$3}' OFS="\t" ${knownSJ}) > ${sample}.tmp.annot.SJ.txt
# use grep to identify only lines containing either a known acceptor or donor site
grep -f ${sample}.tmp.annot.SJ.txt ${sample}.tmp2.exonDel.bed > ${sample}.exonDel.bed

# perform check and exit program at this stage if no deletions detected
if [ ! -s ${sample}.exonDel.bed ]
then
  echo "Only SJ with novel acceptor-donor detected. No deletions to report"
  rm ${sample}.*
  touch ${sample}.exonDel.bed
  echo -e "chr\tgene\tstrand\tdeletion\tn_SJ\tsupport_reads\ttotal_reads\tPRS\tmedian_reads\tRRS\tconfidence\tSJ" > ${sample}.exonDel.anno
  exit
else
  echo -e "Calculating read support"
fi

#####################################
#####################################
###                               ###
###   Step05: sum read support    ###
###                               ###
#####################################
#####################################

## Generate a separate file containing read counts for SJ's involved in deletion
cat ${sample}.unfilt.exonDel.bed | cut -f 1-3 | uniq > ${sample}.tmp.bed

# search original SJ file for any SJ_start matching SJ_start in deletion file
awk 'FNR==NR{a[$2]=$1; next}; $2 in a {print $0;}' ${sample}.tmp.bed $1 | awk '{if($7 > 5 && $7 > $8){print $0}}' > ${sample}.unfilt.count.bed
# repeat for SJ_end
awk 'FNR==NR{a[$3]=$1; next}; $3 in a {print $0;}' ${sample}.tmp.bed $1 | awk '{if($7 > 5 && $7 > $8){print $0}}' >> ${sample}.unfilt.count.bed

## take count.bed file, extract cols representing chr:SJstart, chr:SJend and uniq_map
# create cols containing chr:SJstart and chr:SJend columns
awk -vOFS="\t" '{$10=$1":"$2; $11=$1":"$3}1' ${sample}.unfilt.count.bed | \
# select only the chr:SJstart, chr:SJend and uniq_map cols
awk -vOFS="\t" '{print $10,$11,$7}' | \
# remove duplicated lines & write to new file
awk '!visited[$0]++' > ${sample}.uniq.count.bed

# sum counts for each unique SJstart and SJend
awk -vOFS="\t" '{start[$1]++; sum1[$1]+=$3} END{ for (var in start) print var, sum1[var]}' ${sample}.uniq.count.bed | \
# sort the output by genomic location
# awk -F '[:\t]' '{print $1,$2,$3}' OFS="\t" | sort -t'\t' -nk1,1 -nk2,2 | awk -vOFS="\t" '{print $1":"$2,$3}'
sort -t$'\t' -k1,1 > ${sample}.sum.count.SJstart.bed
# repeat above with SJend
awk -vOFS="\t" '{end[$2]++; sum[$2]+=$3} END{ for (var in end) print var, sum[var]}' ${sample}.uniq.count.bed | \
sort -t$'\t' -k1,1 > ${sample}.sum.count.SJend.bed

#################################################
#################################################
###                                           ###
###   Step05: collate SJ for same deletion    ###
###                                           ###
#################################################
#################################################

echo -e "Collating predictions"

## collate SJ uniq_map representing the same deletion event (e.g. IKZF1 ex 4-8)
awk -F '[:\t]' '{print $5":"$6":"$7,$8}' OFS="\t" ${sample}.exonDel.bed | \
awk -vOFS="\t" '{Del[$1]++; sum[$1]+=$2} END {for (var in Del) print var, Del[var], sum[var]}' | \
# sort for future joining
sort -t$'\t' -k1,1 > ${sample}.tmp1.anno

###################################################
###################################################
###                                             ###
### Step06: identify SJ range for each deletion ###
###                                             ###
###################################################
###################################################

# create unique identifier for each deletion
cat ${sample}.exonDel.bed | \
awk -F '[:\t]' '{ print $1";"$5":"$6":"$7,$2,$4 }' OFS="\t" | \
# use an array to select the min/max start/end of SJ for single deletion event
awk -vOFS="\t" '{Del[$1]++; min[$1]=Del[$1]==1||min[$1]>$2?$2:min[$1]; max[$1]=max[$1]<$3?$3:max[$1]} \
  END { for (var in Del) print var, min[var], max[var]}' | \
# sort deletions by genomic coordinates
awk -F'[;\t]' '{print $2,$1,$3,$4}' OFS="\t" | sort -n -k2 -k3 -k4 > ${sample}.tmp.count.bed

###################################################
###################################################
###                                             ###
###   Step07: calculate median support reads    ###
###                                             ###
###################################################
###################################################

echo -e "Calculating RRS"

### create a subset of SJ.out.tab containing only known SJ

# combine chr:SJstart:SJend into single identifier
  awk -vOFS="\t" '{if($7 > 5 && $7 > $8){print $1":"$2":"$3,$7}}' $1 > ${sample}.tmp.SJ.out.tab
  # exclude unknown/novel SJ
  awk 'NR == FNR {a[$1]; next} ($1 in a)' ${knownSJ} ${sample}.tmp.SJ.out.tab | \
  # expand chr:start:end into separate columns
  awk -F'[:\t]' '{print $1,$2,$3,$4}' OFS="\t" > ${sample}.filt.SJ.out.tab

# create empty count file
echo -n "" > ${sample}.median.count.bed

while read line; do

  # extract variables
  DEL=$(echo $line | awk '{ print $1 }')
  CHR=$(echo $line | awk '{ print $2 }')
  START=$(echo $line | awk '{ print $3 }')
  END=$(echo $line | awk '{ print $4 }')

  # subset SJ.out.tab for SJ's that fall within the range of the potential deletion
  MEDIAN=$(awk -v chr=$CHR -v start=$START -v end=$END -F'\t' '{if($1 == chr && $2 >= start && $3 <= end){print $4}}' ${sample}.filt.SJ.out.tab | \
  # sort support read col
  sort -n | \
  # compute median
  awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }')
  
  # print DEL and MEDIAN to count.bed file
  echo -e "${DEL}\t${MEDIAN}" >> ${sample}.median.count.bed
  # sort median.count.bed for future join
  sort -t$'\t' -k1,1 ${sample}.median.count.bed > ${sample}.sorted.median.count.bed

done < ${sample}.tmp.count.bed

################################
################################
###                          ###
###  Step08: calculate PRS   ###
###                          ###
################################
################################

echo -e "Calculating PRS"

# sort exonDel bed
awk -vOFS="\t" '{print $1,$3,$4,$5}' ${sample}.exonDel.bed | sort -t$'\t' -k1,1 > ${sample}.resort.exonDel.bed
# join sorted exonDel.bed and sum.SJstart
join -t $'\t' -1 1 -2 1 ${sample}.resort.exonDel.bed ${sample}.sum.count.SJstart.bed | \
# remove duplicate totals representing the same start location
awk '!visited[$0]++' | \
# sum read counts supporting the same deletion event
awk -vOFS="\t" '{a[$2"\t"$3"\t"$4]++; sum[$2"\t"$3"\t"$4]+=$5} END {for (var in a) print var, sum[var]}' | \
# sort by gene name and deletion fields
sort -t$'\t' -k1,1 -k3,3 | \
# create a uniq idenfifier for joining
awk -vOFS="\t" '{print $1":"$2":"$3, $4}' | \
# sort for future joining
sort -t$'\t' -k1,1 > ${sample}.tmp2.anno

# repeat for SJend, note must resort sample.exonDel.bed according to SJend location
awk -vOFS="\t" '{print $2,$3,$4,$5}' ${sample}.exonDel.bed | sort -t$'\t' -k1,1 > ${sample}.resort.exonDel.bed
# join resorted exonDel bed with sum.counts.SJend
join -t $'\t' -1 1 -2 1 ${sample}.resort.exonDel.bed ${sample}.sum.count.SJend.bed | \
# remove duplicates
awk '!visited[$0]++' | \
# sum read counts supporting the same deletion event
awk -vOFS="\t" '{a[$2"\t"$3"\t"$4]++; sum[$2"\t"$3"\t"$4]+=$5} END {for (var in a) print var, sum[var]}' | \
# sort by gene name and deletion fields
sort -t$'\t' -k1,1 -k3,3 | \
# create a uniq idenfifier for joining
awk -vOFS="\t" '{print $1":"$2":"$3, $4}' | \
# sort for future joining
sort -t$'\t' -k1,1 > ${sample}.tmp3.anno

################################
################################
###                          ###
###  Step09: calculate RRS   ###
###                          ###
################################
################################

echo -e "Calculating RRS"

# join tmp1.anno with tmp2.anno and then tmp3.anno
join -t $'\t' -1 1 -2 1 ${sample}.tmp1.anno ${sample}.tmp2.anno | join -1 1 -2 1 - ${sample}.tmp3.anno | \
# calculate total reads mapping to SJ of deletion
awk -vOFS="\t" '{print $1,$2,$3,$4,$5,$6=($4+$5-$3)}' | \
# calculate proportion of reads for del SJ and total reads with those SJ
awk -vOFS="\t" '{print $1,$2,$3,$6,$7=($3/$6)}' | \
# join with median read count
join -1 1 -2 1 - ${sample}.sorted.median.count.bed | \
# calculate RRS
awk -vOFS="\t" '{print $1,$2,$3,$4,$5,$6,$7=($3/($6+$3))}' > ${sample}.tmp4.anno

################################################
################################################
###                                          ###
###  Step10: filter & annotate predictions   ###
###                                          ###
################################################
################################################

if [ $norm = "RRS" ]; then

  echo "Filtering/confidence based on RRS"

  # exclude deletions with RRS value < 0.05
  awk -vOFS="\t" '{if($7 >= 0.05){print $0}}' ${sample}.tmp4.anno | \
  # add confidence level based on RRS
  awk -vOFS="\t" '{if($7 < 0.1 || $3 < 11) conf="low"; else if($7 >= 0.2) conf="high"; else conf="moderate"; print $0, conf}' > ${sample}.tmp4.filtered.anno

  else

  echo "Filtering/confidence based on PRS"

  # exclude deletions with PRS value < 0.025
  awk -vOFS="\t" '{if($5 >= 0.025){print $0}}' ${sample}.tmp4.anno | \
  # add confidence level based on PRS
  awk -vOFS="\t" '{if($5 < 0.05 || $3 < 11) conf="low"; else if($5 >= 0.1) conf="high"; else conf="moderate"; print $0, conf}' > ${sample}.tmp4.filtered.anno

fi


# Extract SJ locations from exon.bed file and add these to exonDel.anno
awk -F '[:\t]' '{print $1,$2,$3,$4,$5,$6,$7,$8}' OFS="\t" ${sample}.exonDel.bed | \
awk -vOFS="\t" '{print $9=($5":"$6":"$7";"$1),$10=($2"-"$4)}' OFS="," | \
awk -F, '{a[$1]=(a[$1]?a[$1]FS$2:$2)} END {for (i in a) print i, a[i]}' OFS="\t" | \
awk -F '[;\t]' '{print $1,$2,$3}' OFS="\t" | \
sort -t$'\t' -k1,1 > ${sample}.tmp5.anno

# add SJ locations to sample annotation file
join -t $'\t' -1 1 -2 1 ${sample}.tmp4.filtered.anno ${sample}.tmp5.anno | \
# tidy output and print to final file
awk -F '[:\t]' '{print $11,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$12}' OFS="\t" > ${sample}.exonDel.anno

###########################
###########################
###                     ###
###  Step11: clean-up   ###
###                     ###
###########################
###########################

## Perform a check to confirm that there are lines present in exonDel.anno
if [ ! -s ${sample}.exonDel.anno ]
then
  echo "No deletions detected"
  echo -e "chr\tgene\tstrand\tdeletion\tn_SJ\tsupport_reads\ttotal_reads\tPRS\tmedian_reads\tRRS\tconfidence\tSJ" >> ${sample}.exonDel.anno
else
  echo -e "Deletion detection complete"
  sed -i $'1i chr\tgene\tstrand\tdeletion\tn_SJ\tsupport_reads\ttotal_reads\tPRS\tmedian_reads\tRRS\tconfidence\tSJ' ${sample}.exonDel.anno
fi

## clean up tmp files
rm ${sample}.tmp*
rm ${sample}.unfilt.*
rm ${sample}.resort.exonDel.bed
rm ${sample}.filt.SJ.out.tab
rm ${sample}.median.count.bed
rm ${sample}.sorted.median.count.bed
rm ${sample}.uniq.count.bed
rm ${sample}.sum.count*

# calculate elapsed time and print
elapsed=$(( SECONDS - start_time ))
# report elapsed time
echo -e "Elapsed time (s) for ${sample} : ${elapsed}"
