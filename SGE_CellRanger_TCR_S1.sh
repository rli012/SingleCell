#!/bin/bash
#$ -N S1_TCR
# #$ -V # Job has the same environment variables as the submission shell
#$ -l h=fccccappdevn02
#$ -pe peomics 8 # -pe threads
#$ -l mem_free=8G # h_vmem=4G hard limit of the maximum amount of vitual memory # mem_free=8G, memory needed
#$ -o "logs/S1_TCR.stdout"
#$ -e "logs/S1_TCR.stderr"
#$ -wd "/data/CellRanger/"
# #$ -cwd

# #$ -v "REF=hg38,FQ=test.fq,SAM=test.sam" # no space!!!
# #$ -l h_rt=00:02:00
# #$ -t 1-3 # task array

# usage: qsub SGE_CellRanger_TCR_S1.sh
# OR
# usage: qsub -N S1_TCR -V -pe peomics 8 -l mem_Free=8G -o "logs/S1_TCR.stdout" -e "logs/S1_TCR.stderr" SGE_CellRanger_TCR_S1.sh

echo 'Start'

cellranger vdj --id=S1_TCR \
	       --fastqs=data/FASTQ \
	       --reference=/home/rli3/bin/refdata-cellranger-vdj-GRCh38-alts-ensembl-3.1.0 \
	       --sample=S1_TCR \
	       --localcores=8

echo 'Done'


#:<<SKIP

### TASK ARRAY
#N=$SGE_TASK_ID
#CPU=$NSLOTS

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

