#!/bin/bash 

# Need java and R installed

# sequence folder
SEQUENCE="/home/usr/sequence/"

# prefix of reference genome without .fa or .fasta
# bowtie2-build should be used to index the genome before running this pipeline
GENOME="/home/genome/161010_Chinese_Spring_v1.0_pseudomolecules"  

# key file
KEY="LF_100lines_key.txt"

# minor allele frequency for discovery step
MAF=0.01

DB_NAME_PREFIX="test"     # No .db
THREADS=8                 # for parallel
MEMORY=15                 # memory

# tassel software
TASSEL="/home/usr/bin/tassel-5.2.44-standalone/run_pipeline.pl"

# bowtie2 software
BOWTIE2="/home/usr/bin/bowtie2-2.3.4.2/bowtie2"

################################################## No change needed after here: ############################################### 
DB="${DB_NAME_PREFIX}.db"
TAGS="${DB_NAME_PREFIX}-tags.fastq"
M1="-Xms${MEMORY}G"
M2="-Xmx${MEMORY}G"

# Identifies tags from FASTQ files and store in the database:
/${TASSEL} ${M1} ${M2} -GBSSeqToTagDBPlugin -e PstI-MspI -i ${SEQUENCE} -db ${DB} -k ${KEY} -mnQS 10 -mxKmerNum 200000000  >> ${DB}-log

# Export unique tags (FASTQ)
${TASSEL}  ${M1} ${M2} -TagExportToFastqPlugin -db ${DB} -o ${TAGS} -c 20  >> ${DB}-log

# Alignment with bowtie2
${BOWTIE2} -p ${THREADS} --end-to-end -D 20 -R 3 -N 0 -L 10 -i S,1,0.25 -x ${GENOME} -U ${TAGS} -S ${TAGS}.sam

# Exporting alignment to the database
${TASSEL} ${M1} ${M2} -SAMToGBSdbPlugin -i ${TAGS}.sam -minMAPQ 20 -db ${DB} -aProp 0.64 -aLen 20 >> ${DB}-log

# SNP calling
${TASSEL} ${M1} ${M2} -DiscoverySNPCallerPluginV2 -db ${DB} -deleteOldData true -mnLCov 0.1 -mnMAF ${MAF} >> ${DB}-log

# Export all SNPs in the vcf format
${TASSEL} ${M1} ${M2} -ProductionSNPCallerPluginV2 -db ${DB} -eR 0.0001 -e PstI-MspI -i ${SEQUENCE} -k ${KEY} -minPosQS 0 -do true -mnQS 10 -o ${TAGS}-all.vcf >> ${DB}-log

# Count genotypes: Het, Homozygous allele1, Homozygous allele2, and missing
perl genotype_count_vcf.pl ${TAGS}-all.vcf

# SNP quality stat file
${TASSEL} ${M1} ${M2}  -SNPQualityProfilerPlugin -db ${DB} -statFile stat-${TAGS}.txt >> ${DB}-log

grep -i "ERROR" ${DB}-log >> ERROR.log

# Get SNPs passing at least one filter
bash get_SNP_1filter_passed.sh ${TAGS}



