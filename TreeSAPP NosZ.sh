# TreeSAPP NosZ

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
# treesapp layer \
#  -o SI072_MetaG_contigs_NosZ_assign/ \
#  --refpkg_dir /data/mellaw_data/nosZ/nosZ_MAG_update/final_outputs
