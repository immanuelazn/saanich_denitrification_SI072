# Create a new directory for UniProt-updated NirK processing

cd /root/data/kristi/
mkdir NirK_processed_revised

cd NirK_processed_revised



# Copy UniProt_NirK.fasta file to the new directory

cp /root/data/kristi/NirK_processed/UniProt_NirK.fasta /root/data/kristi/NirK_processed_revised/



# Classifying amino acid sequences 

cd /root/data/kristi/NirK_processed_revised/

treesapp assign \
  -n 4 \
  --trim_align \
  --refpkg_dir /root/data/imman/uniprot_data/NirK_UniProt_update/final_outputs/ \
  --fastx_input UniProt_NirK.fasta \
  --output /root/data/kristi/NirK_processed_revised/UniProt_NirK_assign/



# Download the data package SI072_sequence_data.tar.gz

wget https://zenodo.org/record/6323402/files/SI072_sequence_data.tar.gz && \
tar -xzvf SI072_sequence_data.tar.gz



# Classify ORFs predicted from genomes (SAGs).

treesapp assign \
  -n 4 \
  --trim_align \
  --refpkg_dir /root/data/imman/uniprot_data/NirK_UniProt_update/final_outputs/ \
  --fastx_input /mnt/datasets/2021w/saanich/SAGs/Concatenated_Med_Plus/Saanich_Med_Plus_SAGs.fasta \
  --output SI072_SAGs_assign/

# 
# # OUTPUT:
# #[2022-04-14 21:53:08] DEBUG: 	Initial alignments:	162
# 	Alignments discarded:	61
# 	Fragmented alignments:	70
# 	Inversions detected:	25
# 	Alignments scaffolded:	10
# 	Multi-alignments:	0
# 	Sequences identified:	66
# 
# 	Number of markers identified:
# 		NirK	66
# 
# # [2022-04-14 21:53:09] DEBUG: Number of query sequences in each marker's group:
# NirK	0	58
# NirK	1	4
# NirK	2	3
# NirK	3	1
# 

# Skip updating publicly available sequences 

# Update NirK refpkg with SAG sequences

treesapp update \
  --fast \
  --headless \
  --overwrite \
  --delete \
  --cluster \
  --trim_align \
  -n 4 \
  --output NirK_SAG_update/ \
  --skip_assign \
  --seqs2lineage /mnt/datasets/2021w/saanich/seq2lineage_Tables/Med_Plus_SAGs_GTDB_Taxonomies.tsv \
  --treesapp_output SI072_SAGs_assign/ \
  --refpkg_path /root/data/imman/uniprot_data/NirK_UniProt_update/final_outputs/NirK_build.pkl


# Output: 
# Sequence summary:
# 	Number of sequences: 942
# 	Longest sequence length: 568
# 	Shortest sequence length: 264
# 	Mean sequence length: 410.3
# 	Median sequence length: 383.0
# Extracting information from headers... done.
# Reading cached lineages in '/root/data/kristi/NirK_processed_revised/NirK_SAG_update/intermediates/accession_id_lineage_map.tsv'... done.
# Clustering sequences with MMSeqs' Linclust... done.
# Number of unique lineages:
# 	root       1
# 	domain     3
# 	phylum    22
# 	class     38
# 	order     80
# 	family   142
# 	genus    330
# 	species  634
# Unclassified and incomplete lineages account for 157/905 (17.3%) references.



# Classify ORFs predicted from genomes (MAGs)

treesapp assign \
  -n 4 \
  --trim_align \
  --refpkg_dir /root/data/kristi/NirK_processed_revised/NirK_SAG_update/final_outputs/ \
  --fastx_input /mnt/datasets/2021w/saanich/MAGs/Concatenated/All_SI072_Metawrap_MAGs.fa \
  --output SI072_MAGs_assign/



# Update the NirK reference package with MAG sequences

screen -S NirK_refpkg_MAG_update

treesapp update \
  --fast \
  --headless \
  --overwrite \
  --delete \
  --cluster \
  --trim_align \
  -n 4 \
  --output NirK_MAG_update/ \
  --skip_assign \
  --seqs2lineage /mnt/datasets/2021w/saanich/seq2lineage_Tables/SI072_MAGs_All_GTDB_taxonomies.tsv \
  --treesapp_output SI072_MAGs_assign/ \
  --refpkg_path /root/data/kristi/NirK_processed_revised/NirK_SAG_update/final_outputs/NirK_build.pkl


# Output: 
# [2022-04-16 22:40:21] DEBUG: 0 classified sequences did not meet minimum LWR of 
# 0.1 for updating
# 
# [2022-04-16 22:40:21] INFO:     Number of sequences: 75
#         Longest sequence length: 508
#         Shortest sequence length: 78
#         Mean sequence length: 343.1
#         Median sequence length: 366
# 
# [2022-04-16 22:40:22] DEBUG: Identified and replaced invalid ambiguity character
# s in 0 sequences.
# 
# [2022-04-16 22:40:22] INFO: Sequence summary:
#         Number of sequences: 965
#         Longest sequence length: 568
#         Shortest sequence length: 264
#         Mean sequence length: 410.6
#         Median sequence length: 383
# 
# [2022-04-16 22:40:29] INFO: Number of unique lineages:
#         root       1
#         domain     3
#         phylum    24
#         class     38
#         order     84
#         family   152
#         genus    342
#         species  636





# Assign taxonomic labels to Saanich Inlet metagenomic contigs using the NirK reference package

for f in /mnt/datasets/2021w/saanich/MetaG_Assemblies/SI072_*m_contig.fa; \
do sample=$( basename $f | sed 's/.fa//g')
treesapp assign \
-n 8 \
--trim_align \
--refpkg_dir /root/data/kristi/NirK_processed_revised/NirK_MAG_update/final_outputs/ \
-i $f \
--output SI072_MetaG_contigs_NirK_assign/${sample}_assign; done



# Overwrite previous abundance values with transcripts per million (TPM) values calculated from SI072 datasets ✅

screen -S NirK_abundance

for f in SI072_MetaG_contigs_NirK_assign/SI072_*assign; \
do sample=$( basename $f | sed 's/_contig_assign//g')
treesapp abundance \
-n 8 \
--treesapp_output $f \
--reads /mnt/datasets/2021w/saanich/MetaG_Trim_QC_Reads/${sample}_pe.1.fq.gz \
--reverse /mnt/datasets/2021w/saanich/MetaG_Trim_QC_Reads/${sample}_pe.2.fq.gz \
--report update \
--metric tpm; done


# Output files in SI072_MetaG_contigs_NirK_assign:

SI072_100m_contig_assign  SI072_135m_contig_assign  SI072_200m_contig_assign
SI072_10m_contig_assign   SI072_150m_contig_assign
SI072_120m_contig_assign  SI072_165m_contig_assign



# Calculate TPM values for the 7 metatranscriptomes

screen -S NirK_TPM_7

for f in SI072_MetaG_contigs_NirK_assign/SI072_*assign; \
do sample=$( basename $f | sed 's/_contig_assign//g')
treesapp abundance \
-n 8 \
--treesapp_output $f \
--reads /mnt/datasets/2021w/saanich/MetaT_Raw_Reads/${sample}_MetaT_QC_Filtered.fastq.gz \
--pairing pe \
--metric tpm \
--report append; done


# Annotate NirK query sequences with their respective paralog at all 7 depths

treesapp layer \
 -o /root/data/kristi/NirK_processed_revised/SI072_MetaG_contigs_NirK_assign/SI072_10m_contig_assign \
 --refpkg_dir /root/data/kristi/NirK_processed_revised/NirK_MAG_update/final_outputs/


treesapp layer \
 -o /root/data/kristi/NirK_processed_revised/SI072_MetaG_contigs_NirK_assign/SI072_100m_contig_assign \
 --refpkg_dir /root/data/kristi/NirK_processed_revised/NirK_MAG_update/final_outputs/


treesapp layer \
 -o /root/data/kristi/NirK_processed_revised/SI072_MetaG_contigs_NirK_assign/SI072_120m_contig_assign \
 --refpkg_dir /root/data/kristi/NirK_processed_revised/NirK_MAG_update/final_outputs/


treesapp layer \
 -o /root/data/kristi/NirK_processed_revised/SI072_MetaG_contigs_NirK_assign/SI072_135m_contig_assign \
 --refpkg_dir /root/data/kristi/NirK_processed_revised/NirK_MAG_update/final_outputs/


treesapp layer \
 -o /root/data/kristi/NirK_processed_revised/SI072_MetaG_contigs_NirK_assign/SI072_150m_contig_assign \
 --refpkg_dir /root/data/kristi/NirK_processed_revised/NirK_MAG_update/final_outputs/


treesapp layer \
 -o /root/data/kristi/NirK_processed_revised/SI072_MetaG_contigs_NirK_assign/SI072_165m_contig_assign \
 --refpkg_dir /root/data/kristi/NirK_processed_revised/NirK_MAG_update/final_outputs/


treesapp layer \
 -o /root/data/kristi/NirK_processed_revised/SI072_MetaG_contigs_NirK_assign/SI072_200m_contig_assign \
 --refpkg_dir /root/data/kristi/NirK_processed_revised/NirK_MAG_update/final_outputs/



# Concatenate the 7 classification files 

cat /data/kristi/NirK_processed_revised/SI072_MetaG_contigs_NirK_assign/SI072_*m_contig_assign/final_outputs/layered_classifications.tsv | head -n 1 >SI072_NirK_layered_classifications.tsv

tail -q -n +2 /data/kristi/NirK_processed_revised/SI072_MetaG_contigs_NirK_assign/SI072_*m_contig_assign/final_outputs/layered_classifications.tsv >>SI072_NirK_layered_classifications.tsv



# Copy SI072_NirK_layered_classifications.tsv to local computer

scp root@10.32.204.58:/root/data/kristi/NirK_processed_revised/SI072_NirK_layered_classifications.tsv <local file path>
