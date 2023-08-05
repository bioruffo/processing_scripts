#!/bin/bash

# module load STAR
# module load HTSeq

set -eu -o pipefail


# SETUP VARIABLES
SRADIR="./SRA"
GENOMEDIR="/home/common-data/star/hg38_noalts"
GTF="/home/common-data/ucsc/hg38/gencode.v38.annotation.gtf"
THREADS=16
STRAND=2
# 1 =  first read in pair is in feature strand, '++ --'
# 2 = second read in pair is in feature strand, '+- -+'


# Some quick notes
echo ""
echo "Did you check with a subsample, that the forward (R1) reads are on the same strand as the features?"
echo "https://htseq.readthedocs.io/en/release_0.11.1/count.html#cmdoption-htseq-count-s"
echo "Use infer_experiment.py from the rseqc conda env!"
echo "conda activate rseqc"
echo "infer_experiment.py -r /home/common-data/ucsc/hg38/gencode.v38.annotation.bed -i <SRRname>.Aligned.sortedByCoord.out.bam"
echo ""
echo "Also: did you check with a subsample, that the reads don't need adapter removal?"
echo "Use TrimGalore!"
echo "trim_galore <fastq>"
echo ""


# checking software requirements
echo "$(date +'%T') Checking software requirements..."
echo "fasterq-dump..."
fasterq-dump --help 1> /dev/null
echo "fastqc..."
fastqc --help 1> /dev/null
echo "fastq_screen..."
fastq_screen --help 1> /dev/null
echo "trimmomatic..."
java -jar /opt/biosoftware/trimmomatic/0.39/trimmomatic-0.39.jar -version 1> /dev/null
echo "STAR..."
STAR --help 1> /dev/null
echo "samtools..."
samtools index 1> /dev/null
echo "featureCounts..."
featureCounts -v 1> /dev/null
echo "htseq-count..."
htseq-count -h 1> /dev/null
echo "$(date +'%T') All ok!"
echo ""
echo "create (touch) a file named 'stop' to stop the script before the start of next file's processing."
echo ""
echo "STRAND TO BE USED FOR COUNTING: ${STRAND}"
echo ""

mkdir -p qc
mkdir -p qc/raw
mkdir -p qc/trimmed
mkdir -p trimmed
mkdir -p aligned
mkdir -p counts


# Iterate over directories starting with SRR
for DIR in $SRADIR/SRR*
do
  # Check that DIR is actually a directory
  if [[ -d "${DIR}" ]]
  then

    # stopping point if the user 'touch'-ed a file named 'stop' in the current directory
    if [[ -f "stop" ]]
    then
      echo "$(date +'%T') Stop invoked by the presence of a file named 'stop'."
      exit
    fi
    
    echo "$(date +'%T') Processing ${DIR}"
    BASE=$(basename ${DIR})

    echo "$(date +'%T') Extracting..."
    fasterq-dump -e $THREADS ${SRADIR}/${BASE}
    
    echo "$(date +'%T') Running FastQC..."
    fastqc --outdir ./qc/raw -t 2 ${BASE}_1.fastq ${BASE}_2.fastq
    
    echo "$(date +'%T') Running FastQ-Screen on a subset of reads..."
    fastq_screen ${BASE}_1.fastq --top 100000,2000000 --outdir ./qc/raw --threads $THREADS
    fastq_screen ${BASE}_2.fastq --top 100000,2000000 --outdir ./qc/raw --threads $THREADS

    echo "$(date +'%T') Running Trimmomatic..."
    java -jar /opt/biosoftware/trimmomatic/0.39/trimmomatic-0.39.jar PE \
        -threads $THREADS \
        ${BASE}_1.fastq ${BASE}_2.fastq \
        ./trimmed/${BASE}_1_paired_trim.fq.gz ./trimmed/${BASE}_1_unpaired_trim.fq.gz ./trimmed/${BASE}_2_paired_trim.fq.gz ./trimmed/${BASE}_2_unpaired_trim.fq.gz \
        SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:20
    
    echo "$(date +'%T') Removing the fastq files..."
    rm ${BASE}_1.fastq
    rm ${BASE}_2.fastq
    
    echo "$(date +'%T') Running FastQC on the trimmed files..."
    fastqc --outdir ./qc/trimmed -t 2 ./trimmed/${BASE}_1_paired_trim.fq.gz ./trimmed/${BASE}_2_paired_trim.fq.gz
    
    echo "$(date +'%T') Aligning (and sorting)..."
    STAR --runThreadN $THREADS --genomeDir $GENOMEDIR \
    --readFilesIn ./trimmed/${BASE}_1_paired_trim.fq.gz ./trimmed/${BASE}_2_paired_trim.fq.gz \
    --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "./aligned/${BASE}."
    
    echo "$(date +'%T') Indexing..."
    samtools index -@ $THREADS ./aligned/${BASE}.Aligned.sortedByCoord.out.bam

    echo "$(date +'%T') Counting with featureCounts..."
    featureCounts -p -s $STRAND -T $THREADS -t exon -g gene_id -a $GTF -o ./counts/${BASE}.featureCounts.counts.txt ./aligned/${BASE}.Aligned.sortedByCoord.out.bam
    ## -p            paired
    ## -s $STRAND    stranded (feature strand in R1 -> use '1', feature strand in R2 -> use '2')
    ## -T $THREADS   number of threads
    ## -t exon       use "exon" features in the gtf to count
    ## -g gene_id    group features by "gene_id" (return counts per gene)
    ## -a $GTF       features file
    ## -o            output file

#    echo "$(date +'%T') Counting with htseq-count..."
#    htseq-count --format bam --stranded yes --order pos --max-reads-in-buffer 200000000 ./aligned/${BASE}.Aligned.sortedByCoord.out.bam $GTF > ./counts/${BASE}.htseq.counts.txt

  else
    echo "Skipping ${DIR}"
  fi
done
