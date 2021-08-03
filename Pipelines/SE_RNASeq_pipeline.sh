#!/bin/bash
############################################################################################################################### 
#input parameters
############################################################################################################################### *
input_path=$1
secuencer=$2  #Nextera , TrueSeq2, TrueSeq3
genomeGTFPath=$3
genomeStarIndexPath=$4
thread=$5
kallistoIndexPath=$6
trim=$7
############################################################################################################################### 
source activate Aligment
############################################################################################################################### 
#Folders
############################################################################################################################### 
fastq_path=${input_path}/fastq
fastqc_path=${input_path}/fastQC
trimm_path=${input_path}/trimmed
trimmed_fastqc_path=${input_path}/trimmedFastQc
mapping_path=${input_path}/mapping
sorted_path=${input_path}/sorted
htseq_path=${input_path}/htseq
kallisto_path=${input_path}/kallisto
dexseq_path=${input_path}/dexseq

mkdir $fastqc_path
mkdir $trimm_path
mkdir $trimmed_fastqc_path
mkdir $mapping_path
mkdir $sorted_path
mkdir $htseq_path
mkdir $kallisto_path
mkdir $dexseq_path
############################################################################################################################### 
#Adapter file
############################################################################################################################### 
adapter=""
if [ "$secuencer" = "TruSeq2" ]; then
	adapter="TruSeq2-SE.fa"
elif [ "$secuencer" = "TruSeq3" ]; then
	adapter="TruSeq3-SE.fa"
else
echo " no adaptar specified or found"
fi

adapter="./anaconda3/pkgs/trimmomatic-0.38-1/share/trimmomatic-0.38-1/adapters/"${adapter}
echo ${adapter}


############################################################################################################################### 
#Settings
############################################################################################################################### 
maxMemory=15000000000
trimmingParameters="2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30"
############################################################################################################################### 
#                                                Single end sequencing
###############################################################################################################################
echo "Procesing GTF file for DEXSeq"
#dexseq_prepare_annotation.py ${genomeGTFPath} ${genomeGTFPath}.DEXSeq.chr.gff

for i in $(ls ${fastq_path}/*.fastq)
do
#sample_name=`basename $i|awk -F"_" '{print $1}'`
sample_name=`basename $i|sed 's/.fastq//'`
echo "Processing "${sample_name} 
echo "Running Quality control on raw samples (fastqc)"
fastqc ${fastq_path}/${sample_name}.fastq -t ${thread} -outdir ${fastqc_path}

if [ "$trim" = "TRUE" ]; then 
echo "Removing low quality reads and trimming adapters (trimmomatic)"
trimmomatic SE -threads ${thread} ${fastq_path}/${sample_name}.fastq ${trimm_path}/${sample_name}.fastq ILLUMINACLIP:${adapter}:${trimmingParameters}

echo "Running Quality control after trimming (fastqc)"
fastqc ${trimm_path}/${sample_name}.fastq -t ${thread} -o ${trimmed_fastqc_path}
else
trimm_path=${fastq_path} 
fi

echo "Mapping reads to genome (star)"
mkdir ${mapping_path}/${sample_name}
STAR --runMode alignReads --genomeDir ${genomeStarIndexPath} --runThreadN ${thread} --readFilesIn ${trimm_path}/${sample_name}.fastq --limitBAMsortRAM maxMemory --outSAMtype BAM Unsorted --outFileNamePrefix ${mapping_path}/${sample_name}/${sample_name}

echo "Sorting bam file  output (samtools)"
samtools sort -@ ${thread} -n -o ${sorted_path}/${sample_name}.bam  ${mapping_path}/${sample_name}/${sample_name}Aligned.out.bam

echo "Gene expression quantification (HTSeq)"
htseq-count -f bam -t gene -m union  ${sorted_path}/${sample_name}.bam ${genomeGTFPath} >${htseq_path}/${sample_name}.txt

echo "Transcript expression quantification (Kalisto)"
kallisto quant -i ${kallistoIndexPath} -o ${kallisto_path}/${sample_name} -t ${thread} --single -l 75 -s 20 ${trimm_path}/${sample_name}.fastq

echo "Transcript expression quantification (DEXSeq)"
dexseq_count.py ${genomeGTFPath}.DEXSeq.chr.gff ${sorted_path}/${sample_name}.bam ${dexseq_path}/${sample_name}.txt

done

echo "RNASeq finished"

#####################################################################################################

