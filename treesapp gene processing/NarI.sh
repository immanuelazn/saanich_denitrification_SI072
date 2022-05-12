
#Classify amino acids with TreeSAPP assign
cd/data/jw/

treesapp assign \
  -n 4 \
  --trim_align \
  --refpkg_dir /data/jw/NarI_seed/final_outputs/ \
  --fastx_input UniProt_NarI.fasta \
  --output /data/jw/UniProt_NarI_assign/


# Download the data package SI072_sequence_data.tar.gz 

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

Updating with MAGs

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

# 
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
# 
# 



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

