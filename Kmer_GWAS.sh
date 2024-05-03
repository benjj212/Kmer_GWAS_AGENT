### full pipeline 
### path to data
# contain the raw data as fastq file in this format accessions.fastq
path_raw_data="/media/benji/data_benji/001_data/002_GBS_data/unzip_fastq"
# folder containing all the program needed (bbmap/kmc_v3/ bin/kmer_add_strand_information etc..)
path_programs="/media/benji/data_benji/006_Programs"
# folder containing the reference genomes to be used as fasta including the .fai files 
path_to_genome="/media/benji/data_benji/001_data/001_genomes/Wheat/10genomes"
# fastq of all the reads merged from all accessions
path_full_fastq=/media/benji/data_benji/002_Analysis/014_GWAS_Ci1/raw_reads_merged_3/all_AGENT.fastq

path_name_to_pheno=/media/benji/data_benji/002_Analysis/011_Kmer_GWAS_final/acc_name_to_phenotype.txt

## base folder 
path_base="/media/benji/data_benji/005_manuscript/AGENT/base_data/"
# opening the folder and creating sub_folders 

cd $path_base
mkdir GBS_data_filtered
mkdir Mapping_to_10_genomes
mkdir GWAS_files
mkdir Phenotypes
mkdir Kmer_per_accessions
mkdir scripts

# folder where the filtered fastq gonna be generated 
path_nodup_fastq= $path_base/GBS_data_filtered
# folder where the mapped files (sam and bam files) gonna be created
path_mapped_files=$path_base/Mapping_to_10_genomes
# GWAS folder where the kmer table and other file gonna be written. Also where the results of the GWAS will be created as folder (for each phenotype)
path_gwas_folder=$path_base/GWAS_files
# folder containing the phenotypes as .txt files for each phenotype to run GWAS
path_phenotype=$path_base/Phenotypes
# Folder that will contain the kmers data for each accession in the collection (matching the sequencing file)
path_accessions=$path_base/Kmer_per_accessions
# Folder containing some specific R script to be used 
path_scripts=$path_base/scripts


# Filtering the fastq raw-data using Clumpify (removing duplicate reads)
for fastq_file in $path_raw_data/*.fastq
do 
fastq_only=$(basename "$fastq_file")
accessions="${fastq_only%.*}"
$path_programs/bbmap/clumpify.sh in=$fastq_file out=$path_nodup_fastq/filtered_fastq/$accessions.filtered.fastq zl=9 dedupe=t
done

# Mapping filtered reads to the 12 genomes
for fastq in $path_nodup_fastq/filtered_fastq/*fastq
do
for genomes in $path_to_genome/*.fasta
do
genos=$(basename "$genomes")
genomes_only="${genos%.*}"
accessions=$(basename "$fastq")
filename_only="${accessions%.*}"
### mapping
bwa mem $path_to_genome/$genos $path_nodup_fastq/filtered_fastq/$accessions -t 8 > $path_mapped_files/$filename_only.$genomes_only.sam
### sam to bam
samtools view -bS  $path_mapped_files/$filename_only.$genomes_only.sam >  $path_mapped_files/$filename_only.$genomes_only.bam
### sorting the bam files
samtools sort $path_mapped_files/$filename_only.$genomes_only.bam -o $path_mapped_files/$filename_only.$genomes_only.sorted.bam
done
done


### creating the folder and the txt file with the link to the fastq

Rscript $path_scripts/creating_link_path_folder.R $path_raw_data $path_accessions

# Generating the kmer table 
cd $path_gwas_folder
for files in $path_raw_data/*
do
accessions=$(basename "$files")
filename_only="${accessions%.*}"
filename_only="${filename_only%.*}"
cd $path_accessions/$filename_only
$path_programs/external_programs/kmc_v3 -t2 -k31 -ci2 @$filename_only.txt output_kmc_canon ./ 1> kmc_canon.1 2> kmc_canon.2
$path_programs/external_programs/kmc_v3 -t2 -k31 -ci0 -b @$filename_only.txt output_kmc_all ./ 1> kmc_all.1 2> kmc_all.2
echo @$filename_only.txt
$path_programs/bin/kmers_add_strand_information -c output_kmc_canon -n output_kmc_all -k 31 -o kmers_with_strand_3
done



cd $path_base
ll "$path_accessions" | tail -n +2 | awk -v path="$path_accessions" '{printf "%s/%skmers_with_strand_3\t%s\n", path, $NF, $NF, path, $NF}' > kmers_list_paths_3.txt

awk 'BEGIN{ FS=OFS="|" }{ sub(/.$/, "", $1) }1' kmers_list_paths_3.txt > kmers_list_paths_t.txt
awk 'NR > 2 { print }' < kmers_list_paths_t.txt > kmers_list_paths_3k.txt

$path_programs/bin/list_kmers_found_in_multiple_samples -l kmers_list_paths_3k.txt -k 31 --mac 5 -p 0 -o kmers_to_use_3
### generating kmer_table
$path_programs/bin/build_kmers_table -l kmers_list_paths_3k.txt -k 31 -a kmers_to_use_3 -o kmers_table_3
### generating kinship
$path_programs/bin/emma_kinship_kmers -t kmers_table_3 -k 31 --maf 0.05 > kmers_table_3.kinship
$path_programs/bin/filter_kmers -t kmers_table_3 -k kmers_list.txt -o output.txt

# Running the GWAS

for phenotype in $path_phenotype/*.txt
do
pheno=$(basename "$phenotype")
iso_only="${pheno%.*}"
#running GWAS
conda activate Benji_kmer_gwas
cd $path_gwas_folder/threshold_3_proper/
sudo python2.7 $path_programs/kmers_gwas.py --pheno $phenotype --kmers_table kmers_table_3 -l 31 -p 8 --outdir $iso_only --kmers_for_no_perm_phenotype 100001
done


# Extracting the reads containing the significant kmers
Rscript $path_scripts/Bot_listing_kmer_after_gwas_proper.R $path_gwas_folder/threshold_2/ $path_name_to_pheno

### extracting kmer data based on the list of significant kmers 
for phenotype in $path_phenotype/*.txt
do
pheno=$(basename "$phenotype")
iso_only="${pheno%.*}"
cd $path_gwas_folder/threshold_3_proper/$iso_only/kmers/
$path_programs/bin/filter_kmers -t $path_gwas_folder/threshold_3_proper/kmers_table_3 -k $iso_only.kmer_list.txt -o $iso_only.output.txt
done

# Running this part in R 
path_GWAS=$path_gwas_folder/threshold_3_proper/
Rscript $path_scripts/listing_kmer_fasta.R $path_GWAS

# Retrieve the coordinates of the reads 



for iso in CHE_96224 CHE_97251 IRN_GOR5 GRB_JIW2 KAZ_1b THUN12 ARG_4.2 JPN.Chikara
do 
# transforming the kmer list to proper fasta file 
awk -v OFS="\n" '$1=">" $1' $path_GWAS/$iso/kmers/$iso.kmer_list.fasta > $path_GWAS/$iso/kmers/$iso.kmer.fasta
# fetching reads that contain specific kmers 
$path_programs/fetch_reads_with_kmers-master/new/fetch_reads $path_full_fastq $path_full_fastq $path_GWAS/$iso/kmers/$iso.kmer.fasta 31 $path_GWAS/$iso/kmers/$iso
done

# extract the read name and the sequence of the kmer from the fastq
for iso in CHE_96224 CHE_97251 IRN_GOR5 GRB_JIW2 KAZ_1b THUN12 ARG_4.2 JPN.Chikara
do
cd $path_GWAS/$iso/kmers/
echo $iso.R1.fastq
awk 'NR%4==1 {split(substr($1,2), a, "_"); print a[1], a[2], $2}' $iso.R1.fastq > $iso.read.name.txt
# separate the read name and the kmer sequence into 2 differents files 
awk '{print $1}' $iso.read.name.txt > $iso.read.name.only.txt
awk '{print $2}' $iso.read.name.txt > $iso.read.seq.only.txt
done

#filter the bam file for only the reads in the file

for iso in CHE_96224 CHE_97251 IRN_GOR5 GRB_JIW2 KAZ_1b THUN12 ARG_4.2 JPN.Chikara
do 
for genomes in Genome_Triticum_aestivum Taes_ArinaLrFor_v1_genome Taes_Fielder_v1_genome Taes_Jagger_v1_1_genome Taes_Julius_v1_genome Taes_Lancer_v1_genome Taes_Landmark_v1_v1_genome Taes_Mace_v1_genome Taes_Renan_v1_genome Taes_Stanley_v1_2_genome Taes_SYMattis_v1_1_genome Taes_SynHex_v1_genome
do 
samtools merge -o $path_mapped_files/$genomes.full.bam $path_mapped_files/*$genomes.sorted.bam
samtools view -N $path_GWAS/$iso/kmers/${iso}.read.name.only.txt -b $path_mapped_files/$genomes.full.bam > $path_mapped_files/$genomes.$iso.selected.bam
# transform the bam to bed file
bedtools bamtobed -i $path_mapped_files/$genomes.$iso.selected.bam > $path_mapped_files/$genomes.$iso.selected.bed
done
done











