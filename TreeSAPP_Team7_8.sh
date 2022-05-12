######################################################################################################################################
##### NOSZ ##########################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################################

#Update of refpkgs w new sequences already peformed 
cd /data/mellaw_data/nosZ

#Testing purify of reference package

#STEP 1 SAGS
treesapp purity \
  -n 4 \
  -r /data/imman/uniprot_data/NosZ_UniProt_update/final_outputs/NosZ_build.pkl \
  --extra_info /data/mellaw_data/tutorial/TIGRFAM_info.tsv \
  -i /data/mellaw_data/tutorial/TIGRFAM_seed_named.faa \
  --output NosZ_purity

#Output
# Summarizing assignments for reference package NosZ
# Ortholog        Hits    Leaves  Tree-coverage   Description
# --------------------------------------------------------------------------------
# TIGR04246       12      12      0.9     nitrous-oxide reductase, Sec-dependent
# TIGR04244       15      15      1.2     nitrous-oxide reductase, TAT-dependent

#Classify amino acids using TreeSAPP Assign  
treesapp assign \
  -n 4 \
  --trim_align \
  --refpkg_dir /data/imman/uniprot_data/NosZ_UniProt_update/final_outputs \
  --fastx_input /data/imman/uniprot_data/NosZ_uniprot.fasta \
  --output UniProt_NosZ_assign/


#Classifying ORFs predicted from genomes
treesapp assign \
  -n 4 \
  --trim_align \
  --refpkg_dir /data/imman/uniprot_data/NosZ_UniProt_update/final_outputs \
  --fastx_input /mnt/datasets/2021w/saanich/SAGs/Concatenated_Med_Plus/Saanich_Med_Plus_SAGs.fasta \
  --output SI072_SAGs_NosZ_assign/

#Skip updating publically available sequences 

#Update ref package w SAGs #TODO RESTART HERE
treesapp update \
  --fast \
  --headless \
  --overwrite \
  --delete \
  --cluster \
  --trim_align \
  -n 4 \
  --output nosZ_SAG_update/ \
  --skip_assign \
  --seqs2lineage /mnt/datasets/2021w/saanich/seq2lineage_Tables/Med_Plus_SAGs_GTDB_Taxonomies.tsv \
  --treesapp_output SI072_SAGs_NosZ_assign/ \
  --refpkg_path /data/imman/uniprot_data/NosZ_UniProt_update/final_outputs/NosZ_build.pkl

#Output
# Sequence summary:
#         Number of sequences: 1300
#         Longest sequence length: 661
#         Shortest sequence length: 388
#         Mean sequence length: 587.1
#         Median sequence length: 590.0

# ReferencePackage instance of NosZ (D0601):
#         Molecule type:                                      'prot'
#         TreeSAPP version:                                   '0.11.3'
#         Profile HMM length:                                 '580'
#         Substitution model used for phylogenetic inference: 'LG+G4'
#         Number of reference sequences (leaf nodes):          1281
#         Software used to infer phylogeny:                   'FastTree'
#         Date of last update:                                '2022-04-15'
#         Description:                                        'Nitrous-oxide reductase (TAT-dependent|Sec-dependent)'

#Skip clade annotation 

#STEP 2 MAGS

#Classifying ORFs predicted from genomes (MAGS)
treesapp assign \
  -n 4 \
  --trim_align \
  --refpkg_dir /data/mellaw_data/nosZ/nosZ_SAG_update/final_outputs \
  --fastx_input /mnt/datasets/2021w/saanich/MAGs/Concatenated/All_SI072_Metawrap_MAGs.fa \
  --output SI072_MAGS_NosZ_assign/

#Update ref package w MAGs
treesapp update \
  --fast \
  --headless \
  --overwrite \
  --delete \
  --cluster \
  --trim_align \
  -n 4 \
  --output nosZ_MAG_update/ \
  --skip_assign \
  --seqs2lineage /mnt/datasets/2021w/saanich/seq2lineage_Tables/SI072_MAGs_All_GTDB_taxonomies.tsv \
  --treesapp_output SI072_MAGS_NosZ_assign/ \
  --refpkg_path /data/mellaw_data/nosZ/nosZ_SAG_update/final_outputs/NosZ_build.pkl

#Output
# Sequence summary:
#         Number of sequences: 1320
#         Longest sequence length: 661
#         Shortest sequence length: 393
#         Mean sequence length: 587.7
#         Median sequence length: 590.0

# Number of unique lineages:
#         root       1
#         domain     2
#         phylum    20
#         class     38
#         order     72
#         family   146
#         genus    431
#         species  958

# Summary of the updated reference package:
# ReferencePackage instance of NosZ (D0601):
#         Molecule type:                                      'prot'
#         TreeSAPP version:                                   '0.11.3'
#         Profile HMM length:                                 '580'
#         Substitution model used for phylogenetic inference: 'LG+G4'
#         Number of reference sequences (leaf nodes):          1290
#         Software used to infer phylogeny:                   'FastTree'
#         Date of last update:                                '2022-04-15'
#         Description:                                        'Nitrous-oxide reductase (TAT-dependent|Sec-dependent)'

#STEP 3 Calculating relative abundance

#Assign taxnomic labels to Sannich Inlet metagenomic contigs
mkdir SI072_MetaG_contigs_NosZ_assign_2

for f in /mnt/datasets/2021w/saanich/MetaG_Assemblies/SI072_*m_contig.fa; \
do sample=$( basename $f | sed 's/.fa//g')
treesapp assign \
  -n 8 \
  --trim_align \
  --refpkg_dir /data/mellaw_data/nosZ/nosZ_MAG_update/final_outputs/ \
  -i $f \
  --output SI072_MetaG_contigs_NosZ_assign_2/${sample}_assign; done
  
#Overwrite previous abundance values with TPM values calculated from SI072 datasets

#Run in screen
screen -S overwrite

for f in SI072_MetaG_contigs_NosZ_assign_2/SI072_*assign; \
do sample=$( basename $f | sed 's/_contig_assign//g') 
treesapp abundance \
  -n 8 \
  --treesapp_output $f \
  --reads /mnt/datasets/2021w/saanich/MetaG_Trim_QC_Reads/${sample}_pe.1.fq.gz \
  --reverse /mnt/datasets/2021w/saanich/MetaG_Trim_QC_Reads/${sample}_pe.2.fq.gz \
  --report update \
  --metric tpm; done

#detach from screen
Ctrl+a d
screen -ls 

screen -r

#Calculate TPM values from metatranscriptomes 
for f in SI072_MetaG_contigs_NosZ_assign_2/SI072_*assign; \
do sample=$( basename $f | sed 's/_contig_assign//g') 
treesapp abundance \
  -n 8 \
  --treesapp_output $f \
  --reads /mnt/datasets/2021w/saanich/MetaT_Raw_Reads/${sample}_MetaT_QC_Filtered.fastq.gz \
  --pairing pe \
  --metric tpm \
  --report append; done

# #annotate NosZ query sequences with respective paralod
#Peform for all 7 depths 
for f in SI072_MetaG_contigs_NosZ_assign_2/SI072_*assign; \
do treesapp layer \
  -o $f \
  --refpkg_dir /data/mellaw_data/nosZ/nosZ_MAG_update/final_outputs; done

#concatenate the 7 classication files
cat SI072_MetaG_contigs_NosZ_assign_2/SI072_*m_contig_assign/final_outputs/layered_classifications.tsv \
| head -n 1 > SI072_NosZ_layered_classifications.tsv

tail -q -n +2 SI072_MetaG_contigs_NosZ_assign_2/SI072_*m_contig_assign/final_outputs/layered_classifications.tsv >> SI072_NosZ_layered_classifications.tsv

#Attempt to resolve lack of data for 10m
treesapp assign \
  -n 8 \
  --trim_align \
  --refpkg_dir /data/mellaw_data/nosZ/nosZ_MAG_update/final_outputs/ \
  -i /mnt/datasets/2021w/saanich/MetaG_Assemblies/SI072_10m_contig.fa \
  --output SI072_MetaG_contigs_NosZ_assign_2/SI072_10m_contig_assign/



######################################################################################################################################
##### NORB ##########################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################################

#running on group 7 server:  10.32.204.57  Biome1874

#upload UniProt files to server
scp uniprot_NorB.fasta.gz 425:/data/tn/NorB
scp uniprot_NorB.tab.gz 425:/data/tn/NorB

#upload NorB refpkg to server
scp -r final_outputs/ 425:/data/tn/NorB

#unzip files in working directory (/data/tn/NorB)
gunzip uniprot_NorB.fasta.gz
gunzip uniprot_NorB.tab.gz

#Classifying amino acid sequences 
treesapp assign \
  -n 4 \
  --trim_align \
  --refpkg_dir final_outputs/ \
  --fastx_input uniprot_NorB.fasta \
  --output UniProt_NorB_assign/

#Create seqs_to_lineage files for update
echo -e "SeqID\tOrganism\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies" >UniProt_NorB.tsv
tail -q -n +2 uniprot_NorB.tab >>UniProt_NorB.tsv

#Update the refpkg with UniProt NorB sequences
treesapp update \
  --fast \
  --headless \
  --overwrite \
  --delete \
  --cluster \
  --trim_align \
  --min_taxonomic_rank o \
  -n 4 \
  --output NorB_UniProt_update/ \
  --skip_assign \
  --treesapp_output UniProt_NorB_assign/ \
  --refpkg_path final_outputs/NorB_build.pkl \
  --seqs2lineage UniProt_NorB.tsv



#Classifying ORFs predicted from genomes (SAGs) without --trim_align
treesapp assign \
  -n 4 \
  --refpkg_dir NorB_UniProt_update/final_outputs/ \
  --fastx_input /mnt/datasets/2021w/saanich/SAGs/Concatenated_Med_Plus/Saanich_Med_Plus_SAGs.fasta \
  --output SI072_SAGs_assign/

#Update refpkg with SAG sequences
treesapp update \
  --fast \
  --headless \
  --overwrite \
  --delete \
  --cluster \
  -n 4 \
  --output NorB_SAG_update/ \
  --skip_assign \
  --seqs2lineage /mnt/datasets/2021w/saanich/seq2lineage_Tables/Med_Plus_SAGs_GTDB_Taxonomies.tsv \
  --treesapp_output SI072_SAGs_assign/ \
  --refpkg_path NorB_UniProt_update/final_outputs/NorB_build.pkl

#Output after running treesapp update for SAGs with --trim_align:  
#       Number of sequences: 88
#       Longest sequence length: 480
#       Shortest sequence length: 44
#       Mean sequence length: 285.8
#       Median sequence length: 275.5
#ERROR - commands, line 770:
#No classified sequences exceed minimum length threshold of 486.
#None of the SAG sequences pass the minimum quality control threshold sequence length of 486. Therefore the reference package cannot be updated with SAGs from saanich inlet (for NorB)
#Tried again removing --trim_align and same result:
#       Number of sequences: 88
#       Longest sequence length: 478
#       Shortest sequence length: 45
#       Mean sequence length: 286.7
#       Median sequence length: 276.0
#ERROR - commands, line 770:
#No classified sequences exceed minimum length threshold of 489.


#Classifying ORFs predicted from genomes (MAGs) without --trim_align
treesapp assign \
  -n 4 \
  --refpkg_dir NorB_UniProt_update/final_outputs/ \
  --fastx_input /mnt/datasets/2021w/saanich/ MAGs/Concatenated/All_SI072_Metawrap_MAGs.fa \
  --output SI072_MAGs_assign/

#Update refpkg with MAG sequences without --trim_align
treesapp update \
  --fast \
  --headless \
  --overwrite \
  --delete \
  --cluster \
  -n 4 \
  --output NorB_MAG_update/ \
  --skip_assign \
  --seqs2lineage /mnt/datasets/2021w/saanich/seq2lineage_Tables/SI072_MAGs_All_GTDB_taxonomies.tsv \
  --treesapp_output SI072_MAGs_assign/ \
  --refpkg_path NorB_UniProt_update/final_outputs/NorB_build.pkl

#Same error output after running treesapp update for MAGs with and without --trim_align.
#         Number of sequences: 23
#         Longest sequence length: 471
#         Shortest sequence length: 45
#         Mean sequence length: 312.0
#         Median sequence length: 456
# ERROR - commands, line 770:
# No classified sequences exceed minimum length threshold of 489.


#Run command to look at MAGs classification.tsv
less SI072_MAGs_update_assign/final_outputs/classifications.tsv

#Shows that there are 11 sequences with length over 450 (455-472). These sequences are likely NorB and show up in c__Epsilonproteobacteria; o__Campylobacterales and c__Gammaproteobacteria
#With Connor, tried to change the treesapp command to lower the threshold to 450 length, however, there is a technical problem and treesapp was unable to process the command. The programmer of treesapp (Connor) was notified so that the issue can be resolved from the developer standpoint.



#Assign taxonomic labels to Saanich Inlet metagenomic contigs using the NorB reference #package
for f in /mnt/datasets/2021w/saanich/MetaG_Assemblies/SI072_*m_contig.fa; do sample=$( basename $f | sed 's/.fa//g'); treesapp assign -i $f --refpkg_dir NorB_UniProt_update/final_outputs/ --output SI072_MetaG_contigs_NorB_assign/${sample}_assign --trim_align -n 8; done


# overwrite the previous abundance values with the TPM values calculated from the seven #SI072 datasets
for f in SI072_MetaG_contigs_NorB_assign/SI072_*assign; do sample=$( basename $f | sed 's/_contig_assign//g'); treesapp abundance --treesapp_output $f --reads /mnt/datasets/2021w/saanich/MetaG_Trim_QC_Reads/${sample}_pe.1.fq.gz --reverse /mnt/datasets/2021w/saanich/MetaG_Trim_QC_Reads/${sample}_pe.2.fq.gz -n 8 --report update; done

#calculate TPM values from the seven metatranscriptomes
for f in SI072_MetaG_contigs_NorB_assign/SI072_*assign; do sample=$( basename $f | sed 's/_contig_assign//g'); treesapp abundance --treesapp_output $f --reads /mnt/datasets/2021w/saanich/MetaT_Raw_Reads/${sample}_MetaT_QC_Filtered.fastq.gz -n 8 --pairing pe --metric tpm --report append; done

#annotate the NorB query sequences with their respective paralog at each depth
treesapp layer \
 -o SI072_MetaG_contigs_NorB_assign/SI072_100m_contig_assign \
 --refpkg_dir NorB_UniProt_update/final_outputs/

treesapp layer \
 -o SI072_MetaG_contigs_NorB_assign/SI072_10m_contig_assign \
 --refpkg_dir NorB_UniProt_update/final_outputs/

treesapp layer \
 -o SI072_MetaG_contigs_NorB_assign/SI072_120m_contig_assign \
 --refpkg_dir NorB_UniProt_update/final_outputs/

treesapp layer \
 -o SI072_MetaG_contigs_NorB_assign/SI072_135m_contig_assign \
 --refpkg_dir NorB_UniProt_update/final_outputs/

treesapp layer \
 -o SI072_MetaG_contigs_NorB_assign/SI072_150m_contig_assign \
 --refpkg_dir NorB_UniProt_update/final_outputs/
treesapp layer \
 -o SI072_MetaG_contigs_NorB_assign/SI072_165m_contig_assign \
 --refpkg_dir NorB_UniProt_update/final_outputs/

treesapp layer \
 -o SI072_MetaG_contigs_NorB_assign/SI072_200m_contig_assign \
 --refpkg_dir NorB_UniProt_update/final_outputs/

#combine layered_classification tables into one file
cat /data/tn/NorB/SI072_MetaG_contigs_NorB_assign/SI072_*m_contig_assign/final_outputs/layered_classifications.tsv | head -n 1 >SI072_NorB_layered_classifications.tsv
tail -q -n +2 /data/tn/NorB/SI072_MetaG_contigs_NorB_assign/SI072_*m_contig_assign/final_outputs/layered_classifications.tsv >>SI072_NorB_layered_classifications.tsv


######################################################################################################################################
##### NARI ##########################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################################

#NarI reference package imported from Hallam lab github
#Refpkg file path
# /data/jw/NarI_seed/NarI_build.pkl

#Classify amino acids with TreeSAPP assign
cd/data/jw/

treesapp assign \
  -n 4 \
  --trim_align \
  --refpkg_dir /data/jw/NarI_seed/final_outputs/ \
  --fastx_input UniProt_NarI.fasta \
  --output /data/jw/UniProt_NarI_assign/


# Download the data package SI072_sequence_data.tar.gz from Zenodo and decompress it with tar

wget https://zenodo.org/record/6323402/files/SI072_sequence_data.tar.gz && \
tar -xzvf SI072_sequence_data.tar.gz

# Classify ORFs predicted from genomes (SAGs).
treesapp assign \
  -n 4 \
  --trim_align \
  --refpkg_dir  /data/jw/NarI_seed/final_outputs/ \
  --fastx_input /mnt/datasets/2021w/saanich/SAGs/Concatenated_Med_Plus/Saanich_Med_Plus_SAGs.fasta \
  --output SI072_SAGs_assign/

#Skip updating publicly available sequences

#Update NarI reference package with SAGs
treesapp update \
  --fast \
  --headless \
  --overwrite \
  --delete \
  --cluster \
  --trim_align \
  -n 4 \
  --output NarI_SAG_update/ \
  --skip_assign \
  --seqs2lineage /mnt/datasets/2021w/saanich/seq2lineage_Tables/Med_Plus_SAGs_GTDB_Taxonomies.tsv \
  --treesapp_output SI072_SAGs_assign/ \
  --refpkg_path/data/jw/NarI_seed/final_outputs/NarI_build.pkl

# Classify ORFs predicted from genomes (MAGs) 

#Updating with MAGs
treesapp assign \
  -n 4 \
  --trim_align \
  --refpkg_dir /data/jw/NarI_SAG_update/final_outputs \
  --fastx_input /mnt/datasets/2021w/saanich/MAGs/Concatenated/All_SI072_Metawrap_MAGs.fa \
  --output /data/jw/SI072_MAGs_NarI_assign_version2/

treesapp update \
  --fast \
  --headless \
  --overwrite \
  --delete \
  --cluster \
  --trim_align \
  -n 4 \
  --output NarI_MAG_update_version2/ \
  --skip_assign \
  --seqs2lineage /mnt/datasets/2021w/saanich/seq2lineage_Tables/SI072_MAGs_All_GTDB_taxonomies.tsv \
  --treesapp_output SI072_MAGs_NarI_assign_version2/ \
  --refpkg_path /data/jw/NarI_build.pkl


# Sequence summary:
#       Number of sequences: 242
#       Longest sequence length: 757
#       Shortest sequence length: 451
#       Mean sequence length: 681.9
#       Median sequence length: 694.0
# Extracting information from headers... done.
# Reading cached lineages in '/data/jw/NarI_MAG_update_version2/intermediates/accession_id_lineage_map.tsv'... done.
# Clustering sequences with MMSeqs' Linclust... done.
# Number of unique lineages:
#       root       1
#       domain     2
#       phylum    15
#       class     21
#       order     27
#       family    45
#       genus     96
#       species  165
# Unclassified and incomplete lineages account for 31/210 (14.8%) references.


#Calculating abundance
#Assign taxonomic labels to Sannich Inlet metagenomic contigs
mkdir SI072_MetaG_contigs_NarI_assign_3


for f in /mnt/datasets/2021w/saanich/MetaG_Assemblies/SI072_*m_contig.fa; \
do sample=$( basename $f | sed 's/.fa//g')
treesapp assign \
  -n 8 \
  --trim_align \
  --refpkg_dir /data/jw/NarI_MAG_update/final_outputs/ \
  -i $f \
  --output SI072_MetaG_contigs_NarI_assign_3/${sample}_assign;
done


#Overwrite previous abundance values with TPM values based on SI072 dataset calculations
for f in SI072_MetaG_contigs_NarI_assign_3/SI072_*assign; \
do sample=$( basename $f | sed 's/_contig_assign//g') 	
treesapp abundance \
  -n 8 \
  --treesapp_output $f \
  --reads /mnt/datasets/2021w/saanich/MetaG_Trim_QC_Reads/${sample}_pe.1.fq.gz \
  --reverse /mnt/datasets/2021w/saanich/MetaG_Trim_QC_Reads/${sample}_pe.2.fq.gz \
  --report update \
  --metric tpm; done

#Calculate TPM values from metatranscriptomes 
for f in SI072_MetaG_contigs_NarI_assign_3/SI072_*assign; \
do sample=$( basename $f | sed 's/_contig_assign//g') 
treesapp abundance \
  -n 8 \
  --treesapp_output $f \
  --reads /mnt/datasets/2021w/saanich/MetaT_Raw_Reads/${sample}_MetaT_QC_Filtered.fastq.gz \
  --pairing pe \
  --metric tpm \
  --report append; done


# Annotate NarI query sequences with their respective paralog at all 7 depths
for FILE in /data/jw/SI072_MetaG_contigs_NarI_assign_3/*
do 
	treesapp layer \
		-o $FILE \
		--refpkg_dir /data/jw/NarI_MAG_update_version2/final_outputs/
done



# Concatenate 7 classification files 
cat /data/jw/SI072_MetaG_contigs_NarI_assign_3/SI072_*m_contig_assign/final_outputs/layered_classifications.tsv | head -n 1 >SI072_NarI_layered_classifications.tsv

tail -q -n +2 /data/jw/SI072_MetaG_contigs_NarI_assign_3/SI072_*m_contig_assign/final_outputs/layered_classifications.tsv >>SI072_NarI_layered_classifications.tsv


######################################################################################################################################
##### NIRK ##########################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################################


# Absolute path of NirK reference package (Immanuel updated from UniProt):
# /root/data/imman/uniprot_data/NirK_UniProt_update/final_outputs/NirK_build.pkl
 

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



# Download the data package SI072_sequence_data.tar.gz from Zenodo and decompress it with tar
wget https://zenodo.org/record/6323402/files/SI072_sequence_data.tar.gz && \
tar -xzvf SI072_sequence_data.tar.gz

# Classify ORFs predicted from genomes (SAGs). 
treesapp assign \
  -n 4 \
  --trim_align \
  --refpkg_dir /root/data/imman/uniprot_data/NirK_UniProt_update/final_outputs/ \
  --fastx_input /mnt/datasets/2021w/saanich/SAGs/Concatenated_Med_Plus/Saanich_Med_Plus_SAGs.fasta \
  --output SI072_SAGs_assign/


# OUTPUT:
#[2022-04-14 21:53:08] DEBUG: 	Initial alignments:	162
	# Alignments discarded:	61
	# Fragmented alignments:	70
	# Inversions detected:	25
	# Alignments scaffolded:	10
	# Multi-alignments:	0
	# Sequences identified:	66

	# Number of markers identified:
	# 	NirK	66

# # [2022-04-14 21:53:09] DEBUG: Number of query sequences in each marker's group:
# NirK	0	58
# NirK	1	4
# NirK	2	3
# NirK	3	1


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


# # Output: 
# [2022-04-16 22:40:21] DEBUG: 0 classified sequences did not meet minimum LWR of 
# 0.1 for updating

# [2022-04-16 22:40:21] INFO:     Number of sequences: 75
#         Longest sequence length: 508
#         Shortest sequence length: 78
#         Mean sequence length: 343.1
#         Median sequence length: 366

# [2022-04-16 22:40:22] DEBUG: Identified and replaced invalid ambiguity character
# s in 0 sequences.

# [2022-04-16 22:40:22] INFO: Sequence summary:
#         Number of sequences: 965
#         Longest sequence length: 568
#         Shortest sequence length: 264
#         Mean sequence length: 410.6
#         Median sequence length: 383

# [2022-04-16 22:40:29] INFO: Number of unique lineages:
#         root       1
#         domain     3
#         phylum    24
#         class     38
#         order     84
#         family   152
#         genus    342
#         species  636

 

# Assign taxonomic labels to Saanich Inlet metagenomic contigs using the NirK reference package ✅

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
# SI072_100m_contig_assign  SI072_135m_contig_assign  SI072_200m_contig_assign
# SI072_10m_contig_assign   SI072_150m_contig_assign
# SI072_120m_contig_assign  SI072_165m_contig_assign



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
# scp root@10.32.204.58:/root/data/kristi/NirK_processed_revised/SI072_NirK_layered_classifications.tsv <local file path>



######################################################################################################################################
##### NAPA ##########################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################################

cd /data/jhe/napA
# All commands were run in a script after defining the shell #!/bin/bash
# command was run as follows:
	nohup bash <script> &
 
# Classify SAGs with TreeSAPP assign
treesapp assign \
	-n 4  \
	--trim_align \
	--refpkg_dir /root/Refpkgs/Nitrogen_metabolism/Denitrification/NapA/seed_refpkg/final_outputs \
	--fastx_input /mnt/datasets/2021w/saanich/SAGs/Concatenated_Med_Plus/Saanich_Med_Plus_SAGs.fasta \
	--output SI072_SAGs_assign/

# Classify MAGs with TreeSAPP assign
treesapp assign \
	-n 4  \
	--trim_align \
	--refpkg_dir /root/Refpkgs/Nitrogen_metabolism/Denitrification/NapA/seed_refpkg/final_outputs \
	--fastx_input /mnt/datasets/2021w/saanich/MAGs/Concatenated/All_SI072_Metawrap_MAGs.fa \
	--output SI072_MAGs_assign/

# Updating reference packages with SAGs
treesapp update  \
	--fast \
	--headless \
	--overwrite \
	--delete \
	--cluster \
	--trim_align \
	-n 4 \
	--output NapA_SAG_update/ \
	--skip_assign \
	--seqs2lineage /mnt/datasets/2021w/saanich/seq2lineage_Tables/Med_Plus_SAGs_GTDB_Taxonomies.tsv \
	--treesapp_output SI072_SAGs_assign \
	--refpkg_path /root/Refpkgs/Nitrogen_metabolism/Denitrification/NapA/seed_refpkg/final_outputs/NapA_build.pkl \

# Updating reference packages with MAGs (after SAG update)
treesapp update  \
	--fast \
	--headless \
	--overwrite \
	--delete \
	--cluster \
	--trim_align \
	-n 4 \
	--output NapA_MAG_update/ \
	--skip_assign \
	--seqs2lineage /mnt/datasets/2021w/saanich/seq2lineage_Tables/SI072_MAGs_All_GTDB_taxonomies.tsv \
	--treesapp_output SI072_MAGs_assign \
	--refpkg_path /data/jhe/napA/NapA_SAG_update/final_outputs/NapA_build.pkl \

# Check purity of updated reference packages (most recent update was with MAGs)
treesapp purity \
	-n 4 \
	-r NapA_MAG_update/final_outputs/NapA_build.pkl \
	--extra_info /data/jhe/tutorial/TIGRFAM_info.tsv \
	-i /data/jhe/tutorial/TIGRFAM_seed_named.faa \
	--output NapA_purity

for FILE in /mnt/datasets/2021w/saanich/MetaG_Assemblies/SI072_*m_contig.fa
do 
	sample=$( basename $FILE | sed 's/.fa//g');
 
	treesapp assign \
		--fastx_input $FILE \
		--refpkg_dir /data/jhe/napA/NapA_MAG_update/final_outputs/ \
		--output /data/jhe/napA/SI072_MetaG_contigs_NapA_assign/${sample}_assign \
		--trim_align \
		-n 8
done

# Calculate TPM metrics for metgenome reads
for FILE in /mnt/datasets/2021w/saanich/MetaG_Trim_QC_Reads/SI072_*m_pe.1.fq.gz
do 
	sample=$( basename $FILE | sed 's/.1.fq.gz//g');
	depth=$( basename $FILE | sed 's/_pe.1.fq.gz//g');
 
	treesapp abundance \
		-n 4 \
		--treesapp_output SI072_MetaG_contigs_NapA_assign/${depth}_contig_assign/ \
		--reads $FILE \
		--reverse /mnt/datasets/2021w/saanich/MetaG_Trim_QC_Reads/${sample}.2.fq.gz \
		--report update \
		--metric tpm
done
 
 
 
# Calculate TPM metrics for metatranscriptome reads
for FILE in /mnt/datasets/2021w/saanich/MetaT_Raw_Reads/*
do 
	depth=$( basename $FILE | sed 's/_MetaT_QC_Filtered.fastq.gz//g');
 
	treesapp abundance \
		-n 4 \
		--treesapp_output SI072_MetaG_contigs_NapA_assign/${depth}_contig_assign/ \
		--reads $FILE \
		--pairing pe \
		--report append \
		--metric tpm
done
 

# Annotate NapA query sequences with their paralogs
for FILE in /data/jhe/napA/SI072_MetaG_contigs_NapA_assign/*
do 
	treesapp layer \
		-o $FILE \
		--refpkg_dir NapA_MAG_update/final_outputs/
done


######################################################################################################################################
##### CHECK PURITY OF PKL FILES ##########################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################################################

cd /root/Refpkgs/Nitrogen_metabolism/Denitrification/pkl_only/
for file in /root/Refpkgs/Nitrogen_metabolism/Denitrification/pkl_only/*.pkl; do \
treesapp purity \
  -n 4 \
  -r $file \
  --extra_info /root/data/imman/TIGRFAM_info.tsv \
  -i /root/data/imman/TIGRFAM_seed_named.faa \
  --output "${file%*_build.pkl}"_purity; done


# Retrieved the uniprot databasefiles manually from uniprot using:
# "nitrate reductase" napa NOT taxonomy:"uncultured bacterium" (137 reviewed, 5965 unreviewed)
# "nitrate reductase" gene:nari NOT taxonomy:"uncultured bacterium" (2 reviewed, 5777 unreviewed)
# "copper-containing nitrite reductase" gene:nirk NOT taxonomy:"uncultured bacterium"  (7 reviewed, 3046 unreviewed)
# "nitrous-oxide reductase" gene:nosz NOT taxonomy:"uncultured bacterium" (12 reviewed, 5707 unreviewed)
# "nitric oxide reductase subunit b" gene:norb NOT taxonomy:"uncultured bacterium" (2 reviewed, 515 unreviewed)
# Named each as as "$gene_name"_uniprot.fasta
# placed in /data/imman/uniprot_data


# run treeSAPP assign on all of the retrieved UniProt files in Denitrification pathway
# Note: we had to manually move the pkl files to the front of each directory for this to work.
#/root/data/imman/uniprot_data
for file in *.fasta; do \
treesapp assign \
  -n 4 \
  --trim_align \x
  --refpkg_dir /root/Refpkgs/Nitrogen_metabolism/Denitrification/"${file%*_uniprot.fasta}" \
  --fastx_input $file \
  --output UniProt_"${file%*_uniprot.fasta}"_assign/ ;
done | tee uniprot_assign_cmd_output.txt
\


######################################################################################################################################
##### Update pkl with Uniprot data ##########################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################################################

#Create seqs_to_lineage files for update
cd /root/data/imman/uniprot_data
for uniprot_out in */; do 
uniprot_out_no_ext=${uniprot_out%*_assign/} 
uniprot_clean=${uniprot_out_no_ext#*UniProt_}
echo -e "SeqID\tOrganism\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies" > "$uniprot_out"UniProt_"$uniprot_clean".tsv
tail -q -n +2 "$uniprot_clean"_uniprot.tab >> "$uniprot_out"UniProt_"$uniprot_clean".tsv; done

# Run update for all refpkgs
for uniprot_out in *_assign/; do
uniprot_out_no_ext=${uniprot_out%*_assign/} 
uniprot_clean=${uniprot_out_no_ext#*UniProt_}
treesapp update \
  --fast \
  --headless \
  --overwrite \
  --delete \
  --cluster \
  --trim_align \
  --min_taxonomic_rank o \
  -n 4 \
  --output "$uniprot_clean"_UniProt_update/ \
  --skip_assign \
  --treesapp_output UniProt_"$uniprot_clean"_assign/ \
  --refpkg_path /root/Refpkgs/Nitrogen_metabolism/Denitrification/pkl_only/"$uniprot_clean"_build.pkl \
  --seqs2lineage UniProt_"$uniprot_clean"_assign/UniProt_"$uniprot_clean".tsv; done

######################################################################################################################################
##### Cloning Reference Packages ##########################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################################
#####################################################################################################################################################################################################################################################

cd ~
mkdir Refpkgs
cd Refpkgs

#  Pull only 1 depth to make sure it actually works -- crashing when pulling with all depths
git init
git remote add origin https://github.com/hallamlab/RefPkgs.git
git config core.sparseCheckout true
echo "Nitrogen_metabolism/Denitrification" >> .git/info/sparse-checkout
echo "TreeSAPP_reference_package_creation.tsv" >> .git/info/sparse-checkout
echo "refpkg_manager.py" >> .git/info/sparse-checkout
git pull --depth=1 origin master


#
cd /data/imman
for fold in /root/Refpkgs/Nitrogen_metabolism/Denitrification/*; do 
treesapp package \
  view tree \
  --refpkg_path "$fold"/seed_refpkg/final_outputs/"${fold#*/root/Refpkgs/Nitrogen_metabolism/Denitrification/}"_build.pkl > \
  "${fold#*/root/Refpkgs/Nitrogen_metabolism/Denitrification/}"_labelled_tree.txt;
Done

# Created a folder including only .pkl files of NorB, NosZ, NarI, NirK, and NapA
