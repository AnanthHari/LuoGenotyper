#!/bin/bash

# Do for each gene

# inputs
gene=$1
sample_name=$2
igblast_output_prefix=$3
spades_output_prefix=$4
bowtie2_output_prefix=$5
gatk_output_prefix=$6
hapcut2_output_prefix=$7

# run bowtie2 against the contigs
if [ ! -d $bowtie2_output_prefix ]
then
  mkdir -p $bowtie2_output_prefix
fi
out_prefix="$bowtie2_output_prefix"/"$gene"
ref="$spades_output_prefix"/"$gene"/contigs.fasta
index="$out_prefix"/index
bowtie2_bam="$out_prefix"/"$gene"-contigs.bam
mkdir -p $out_prefix
# create bowtie2 index
bowtie2-build -f $ref $index
# map filtered reads
bowtie2 --local --score-min G,20,30 -x $index -U "$igblast_output_prefix"/"$gene".fq | samtools view -b -o $bowtie2_bam


# GATK
if [ ! -d $gatk_output_prefix ]
then
  mkdir -p $gatk_output_prefix
fi
gatk_bam=$(echo $bowtie2_bam | sed 's/bowtie2/gatk/')
out_prefix="$gatk_output_prefix"/"$gene"
vcf="$out_prefix"/"$gene"-vs-ref.vcf
mkdir -p $out_prefix
# get .fai and .dict files
samtools faidx $ref
gatk CreateSequenceDictionary -R $ref
# prep the bam file
gatk AddOrReplaceReadGroups \
  --INPUT $bowtie2_bam \
  --OUTPUT $gatk_bam \
  --RGID $sample_name \
  --RGLB lib1 \
  --RGPL ILLUMINA \
  --RGPU unit1 \
  --RGSM $sample_name \
  --SORT_ORDER coordinate
samtools index $gatk_bam
# run gatk
gatk HaplotypeCaller -ploidy 2 --stand-call-conf 30 -I $gatk_bam -R $ref -O $vcf

# HapCUT2
if [ ! -d $gatk_output_prefix ]
then
  mkdir -p $gatk_output_prefix
fi
out_prefix="$hapcut2_output_prefix"/"$gene"
frag_file="$out_prefix"/"$gene"-fragment.txt
hap_file="$out_prefix"/"$gene"-haplotype.txt
mkdir -p $out_prefix
# build fragment file
extractHAIRS --bam $gatk_bam --VCF $vcf --out $frag_file
# get haplotypes
HAPCUT2 --fragments $frag_file --VCF $vcf --output $hap_file
