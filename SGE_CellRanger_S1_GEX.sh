#!/bin/bash
#$ -N S1_GEX
# #$ -V # Job has the same environment variables as the submission shell
#$ -l h=highmem
#$ -pe peomics 8 # -pe threads
#$ -l mem_free=8G # h_vmem=4G hard limit of the maximum amount of vitual memory
#$ -o "logs/S1_GEX.stdout"
#$ -e "logs/S1_GEX.stderr"
#$ -wd "/data/CellRanger/"
# #$ -cwd

# #$ -v "REF=hg38,FQ=test.fq,SAM=test.sam" # no space!!!
# #$ -l h_rt=00:02:00
# #$ -t 1-3 # task array

# usage: qsub SGE_CellRanger_S1_GEX.sh
# OR
# usage: qsub -N S1_GEX -V -pe peomics 8 -l mem_free=8G -o "logs/S1_GEX.stdout" -e "logs/S1_GEX.stderr" SGE_CellRanger_S1_GEX.sh

echo 'Start'
cellranger count --id=S1_GEX \
		 --fastqs=/data/FASTQ \
		 --transcriptome=/home/rli3/bin/refdata-cellranger-GRCh38-3.0.0 \
		 --sample=S1_GEX \
		 --localcores=8

echo 'Done'


#:<<SKIP

### TASK ARRAY
#N=$SGE_TASK_ID

#echo $CPU

CPU=$NSLOTS

#FILE=`ls raw/test*\.fastq.gz | grep _1.fastq.gz | head -n $N | tail -n 1`
#SAMPLE=${FILE%_1.fastq.gz}
#SAMPLE=${SAMPLE#raw/}
#fq1=$FILE
#fq2=${FILE/_1/_2}


#echo 'Start'


#echo ${REF}
#echo $FQ
#echo $SAM

#echo $fq1
#echo $fq2

#echo 'Done'

#SKIP

