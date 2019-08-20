#!/bin/bash
#$ -N S5_S1_GEX
# #$ -V # Job has the same environment variables as the submission shell
#$ -l h=fccccappdevn03
#$ -pe threads 8 # -pe peomics 2 ??? # -pe threads
#$ -l mem_free=8G # h_vmem=4G hard limit of the maximum amount of vitual memory
#$ -o "logs/S5_S1_GEX.stdout"
#$ -e "logs/S5_S1_GEX.stderr"
#$ -wd "/home/rli3/Projects/HPK1_Single_Cell_RNAseq_20190812/data/CellRanger/"
# #$ -cwd

# #$ -v "REF=hg38,FQ=test.fq,SAM=test.sam" # no space!!!
# #$ -l h_rt=00:02:00
# #$ -wd "/home/rli3/Documents"
# #$ -t 1-3 # task array

# usage: qsub SGE_qsub_template.sh
# OR
# usage: qsub -N run_test -V -pe threads 2 -l h_vmem=4G -o "logs/test.stdout" -e "logs/test.stderr" SGE_qsub_template.sh

CPU=$NSLOTS

echo 'Start'
#echo $CPU

#echo ${REF}
#echo $FQ
#echo $SAM

cellranger count --id=S1_GEX \
		 --fastqs=/resbioinfo/data/MedGenome/P2000413_07292019/FASTQ \
		 --transcriptome=/home/rli3/bin/refdata-cellranger-GRCh38-3.0.0 \
		 --sample=S1_GEX \
		 --localcores=8

echo 'Done'



#:<<SKIP

### TASK ARRAY
#N=$SGE_TASK_ID

#FILE=`ls raw/test*\.fastq.gz | grep _1.fastq.gz | head -n $N | tail -n 1`
#SAMPLE=${FILE%_1.fastq.gz}
#SAMPLE=${SAMPLE#raw/}
#fq1=$FILE
#fq2=${FILE/_1/_2}


#echo 'Start'

#echo $fq1
#echo $fq2

#echo 'Done'

#SKIP

