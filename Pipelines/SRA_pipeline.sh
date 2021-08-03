
#########################
#SRA
#########################
input_path=$1
output_path=$2

#fastq-dump también lo puedo hacer el retrieve del SRA

#SRA -->fastq

for i in $(ls ${input_path}/*.sra)
do
#sample_name=`basename $i|awk -F"_" '{print $1}'`
sample_name=`basename $i|sed 's/.sra//'`
echo "Processing "${sample_name}

fastq-dump --split-3 ${input_path}/${sample_name}.sra -O ${output_path}
done

echo "SRA finished"
