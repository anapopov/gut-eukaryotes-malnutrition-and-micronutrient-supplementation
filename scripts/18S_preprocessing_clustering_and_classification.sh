#!/bin/bash

##########################################################################
#
#	Script for pre-processing of paired amplicon reads, clustering and classification.
#
#
#	Create file list.txt with sample names corresponding to file names.
#	e.g. Sample1_R1.fastq Sample1_R2.fastq Sample2_R1.fastq ... 
#	list.txt:	Sample1
#			Sample2
#
# 	Create primer.fasta file containing amplicon primer sequences
#
#	Environment: Linux
#	Locally installed tools: VSEARCH, USEARCH, TRIMMOMATIC, cutadapt
#	Downloaded databases: SILVA SSU reference fasta file
##########################################################################

ls *.fastq.gz | sed 's/_R.*$//g' | sort | uniq > list.txt

### Create directories
mkdir OUT
mkdir tmp

if [ -f "list.txt" ]; then
	echo "Checking list."
else
	echo 'Sample list file "list.txt" does not exist.'
	exit
fi

NUMFILES=$(wc -l < list.txt)


### Rename all reads to reflect sample ID and read number
for i in $(seq 1 $NUMFILES); do 

	FILE=$(sed -n ${i}p list.txt)
	echo $FILE

gunzip ${FILE}_R1.fastq.gz
gunzip ${FILE}_R2.fastq.gz

	awk -v var=$FILE '/@M/ {print "@" var ":" ++i; next} {print}' ${FILE}_R1.fastq > ${FILE}_R1mod.fastq
	awk -v var=$FILE '/@M/ {print "@" var ":" ++i; next} {print}' ${FILE}_R2.fastq > ${FILE}_R2mod.fastq

gzip ${FILE}_R1.fastq
gzip ${FILE}_R2.fastq

done


### Trim low quality sequences
for i in $(seq 1 $NUMFILES); do

        FILE=$(sed -n ${i}p list.txt)
        echo $FILE

	#paired
	java -jar ~/Work/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -trimlog OUT/${FILE}_PEtrim.log ${FILE}_R1mod.fastq ${FILE}_R2mod.fastq OUT/${FILE}_R1trim_paired200.fastq OUT/${FILE}_R1trim_unpaired200.fastq OUT/${FILE}_R2trim_paired200.fastq OUT/${FILE}_R2trim_unpaired200.fastq ILLUMINACLIP:primers.fasta:2:20:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:200 2>> OUT/trimmomatic.log

done


### Merge or join R1 and R2 reads from paired files
for i in $(seq 1 $NUMFILES); do 

	FILE=$(sed -n ${i}p list.txt)
	echo $FILE

	echo "Merging reads."
	vsearch --fastq_mergepairs OUT/${FILE}_R1trim_paired200.fastq --reverse OUT/${FILE}_R2trim_paired200.fastq --fastq_ascii 33 --fastqout OUT/${FILE}_merged.fastq --fastq_allowmergestagger --threads 39 --fastqout_notmerged_fwd OUT/${FILE}_notmerged_R1.fastq --fastqout_notmerged_rev OUT/${FILE}_notmerged_R2.fastq --log OUT/${FILE}_merged.log 2>> OUT/Merging.log

	echo "Joining unmerged reads with linker."
	usearch11.0.667 -fastq_join ${FILE}_notme
rged_R1.fastq -reverse ${FILE}_notmerged_R2.fastq -fastqout ${FILE}_joined.fastq
 -join_padgap NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN -join_padgapq I
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

done



# concatenate all sample files
cat OUT/Pa*_merged.fastq > OUT/Samples_merged.fastq
cat OUT/Pa*_joined.fastq > OUT/Samples_joined.fastq



# convert fastq to fasta
vsearch --fastq_filter OUT/Samples_merged.fastq --fastaout OUT/Samples_merged.fasta --fasta_width 0

vsearch --fastq_filter OUT/Samples_joined.fastq --fastaout OUT/Samples_joined.fasta --fasta_width 0



# dereplicate reads
vsearch --derep_fulllength OUT/Samples_merged.fasta --output OUT/Samples_merged_derep.fasta --sizeout --fasta_width 0 --relabel seq --log OUT/Samples_merged_derep.log 2>> OUT/Samples_merged_derep.log2

vsearch --derep_fulllength OUT/Samples_joined.fasta --output OUT/Samples_joined_derep.fasta --sizeout --fasta_width 0 --relabel seqjoin --log OUT/Samples_joined_derep.log 2>> OUT/Samples_joined_derep.log2



# pick OTUs with USEARCH UPARSE (97%ID fixed) (automatically removes chimera)
usearch11.0.667 -cluster_otus OUT/Samples_merged_derep_sorted.fasta -otus OUT/Samples_merged_uclust_min5.fasta -relabel otu -minsize 5

usearch11.0.667 -cluster_otus OUT/Samples_joined_derep_sorted.fasta -otus OUT/Samples_joined_uclust_min5.fasta -relabel otu -minsize 5



# assign reads to OTUs
/home/j/jparkin/apopovic/tools/usearch11.0.667 -otutab Samples_merged.fasta -otus Samples_merged_uclust_min5.fasta -id 0.97 -otutabout Samples_merged_uclust_min5_otu-table_97.txt -mapout Samples_merged_uclust_min5_otu-map_97.txt -dbmatched Samples_merged_uclust_min5_otu-matched_97.fasta -sizeout -strand both -threads 39

/home/j/jparkin/apopovic/tools/usearch11.0.667 -otutab Samples_joined.fasta -otus Samples_joined_uclust_min5.fasta -id 0.97 -otutabout Samples_joined_uclust_min5_otu-table_97.txt -mapout Samples_joined_uclust_min5_otu-map_97.txt -dbmatched Samples_joined_uclust_min5_otu-matched_97.fasta -sizeout -strand both -threads 39



# Classification I. SINA
# submit UCLUST otus (merged & joined) to SINA v1.2.11 web-terminal (min-neighbour 1, min % ID 90, similarity cutoff 70%, SILVA v138 Eukaryota)


# Classification II. USEARCH sintax with SILVA v138 database
# 1. download SILVA RefNR99 v138 dataset
# 2. extract eukaryotes only, unwrap sequences, convert from RNA to DNA
# 3. trim V4+V5 amplicon region using cutadapt (10% error)
	cutadapt -e 0.1 -g GCGGTAATTCCAGCTC SILVA_138_SSURef_NR99_tax_silva.fasta > temp.fasta 2> cutadaptV4V5.log
	cutadapt -e 0.1 -a GGAATTGACGGAAKGGC temp.fasta > SILVA_138_SSURef_NR99_V4V5_Eukaryota.fasta 2>> cutadaptV4V5.log
# 4. dereplicate
	vsearch --derep_fulllength SILVA_138_SSURef_NR99_V4V5_Eukaryota.fasta --output SILVA_138_SSURef_NR99_V4V5_Eukaryota_dereplicated.fasta --fasta_width 0
# 5. cluster to 97% using USEARCH UPARSE
# 6. manually added Enterocytozoon, Encephalitozoon and Entamoeba coli V4V5 regions from NCBI nr
# 7. create USEARCH sintax database
usearch11.0.667 -makeudb_usearch SILVA_138_SSURef_NR99_V4V5_Eukaryota_dereplicated_clust97.fasta -output SILVA_138_SSURef_NR99_V4V5_Eukaryota_clust97.udb


# run OTU classification
usearch11.0.667 -sintax Samples_merged_uclust_min5_otu-matched_97.fasta -db /media/sf_VirtualBox_Files/ReferenceDB/SilvaDB/v138/SILVA_138_SSURef_NR99_V4V5_Eukaryota_clust97.udb -tabbedout Samples_merged_derep_USEARCH_classified_SILVA138V4V5_clust97_MOD.txt -strand both -sintax_cutoff 0.8

usearch11.0.667 -sintax Samples_joined_uclust_min5_otu-matched_97.fasta -db /media/sf_VirtualBox_Files/ReferenceDB/SilvaDB/v138/SILVA_138_SSURef_NR99_V4V5_Eukaryota_clust97.udb -tabbedout Samples_merged_derep_USEARCH_classified_SILVA138V4V5_clust97_MOD.txt -strand both -sintax_cutoff 0.8


# Classification III. unassigned OTUs are submitted for BLAST against the NCBI nt database (downloaded Nov 28, 2017)

# all three sets of results merged, and classification method label added



echo "Complete."
