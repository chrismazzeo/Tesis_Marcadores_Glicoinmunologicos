#!/bin/bash
############################################################################################################################### 
#input parameters
############################################################################################################################### 
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
#Folders/
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
mkdir $sorted_path/Users/christianmazzeo
mkdir $htseq_path
mkdir $kallisto_path
mkdir $dexseq_path
############################################################################################################################### 
#Adapter file
############################################################################################################################### 
adapter=""
if [ "$secuencer" = "Nextera" ]; then 
	adapter="NexteraPE-PE.fa"
elif [ "$secuencer" = "TruSeq2" ]; then
	adapter="TruSeq2-PE.fa"
elif [ "$secuencer" = "TruSeq3" ]; then
	adapter="TruSeq3-PE.fa"
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
#                                                Paired end sequencing
###############################################################################################################################
for i in $(ls ${fastq_path}/*_1.fastq)
do

##sample_name=`basename $i|awk -F"_" '{print $1}'`
sample_name=`basename $i|sed 's/_1.fastq//'`

echo "Processing "${sample_name} 
echo "Running Quality control on raw samples (fastqc)"
fastqc ${fastq_path}/${sample_name}_1.fastq -t ${thread} -o ${fastqc_path}
fastqc ${fastq_path}/${sample_name}_2.fastq -t ${thread} -o ${fastqc_path}

if [ "$trim" = "TRUE" ]; then 
echo "Removing low quality reads and trimming adapters (trimmomatic)"
trimmomatic PE -threads ${thread} ${fastq_path}/${sample_name}_1.fastq ${fastq_path}/${sample_name}_2.fastq  ${trimm_path}/${sample_name}_1.fastq ${trimm_path}/${sample_name}_unpaired_1.fastq ${trimm_path}/${sample_name}_2.fastq ${trimm_path}/${sample_name}_unpaired_2.fastq ILLUMINACLIP:${adapter}:${trimmingParameters}

echo "Running Quality control after trimming (fastqc)"
fastqc ${trimm_path}/${sample_name}_1.fastq -t ${thread} -o ${trimmed_fastqc_path}
fastqc ${trimm_path}/${sample_name}_2.fastq -t ${thread} -o ${trimmed_fastqc_path}
else
trimm_path=${fastq_path} 
fi


echo "Mapping reads to genome (star)"
mkdir ${mapping_path}/${sample_name}
STAR --runMode alignReads --genomeDir ${genomeStarIndexPath} --runThreadN ${thread} --readFilesIn ${trimm_path}/${sample_name}_1.fastq 
${trimm_path}/${sample_name}_2.fastq --limitBAMsortRAM maxMemory --outSAMtype BAM Unsorted --outFileNamePrefix ${mapping_path}/${sample_name}/${sample_name}


echo "Sorting bam file  output (samtools)"
samtools sort -@ ${thread} -n -o ${sorted_path}/${sample_name}.bam ${mapping_path}/${sample_name}/${sample_name}Aligned.out.bam

echo "Gene expression quantification (HTSeq)"
htseq-count -f bam -t gene -m union -s $strand ${sorted_path}/${sample_name}.bam ${genomeGTFPath} >${htseq_path}/${sample_name}.txt

echo "Transcript expression quantification (Kalisto)"
kallisto quant -i ${kallistoIndexPath} -o ${kallisto_path}/${sample_name} -t ${thread}  ${trimm_path}/${sample_name}_1.fastq ${trimm_path}/${sample_name}_2.fastq


echo "Transcript expression quantification (DEXSeq)"
dexseq_count.py ${genomeGTFPath}.DEXSeq.chr.gff ${sorted_path}/${sample_name}.bam ${dexseq_path}/${sample_name}.txt

done
echo "RNASeq finished"


#####################################################################################################

