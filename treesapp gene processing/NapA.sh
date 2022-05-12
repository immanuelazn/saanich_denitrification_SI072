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

# Classify MAGs with TreeSAPP assign
treesapp assign \
	-n 4  \
	--trim_align \
	--refpkg_dir /root/Refpkgs/Nitrogen_metabolism/Denitrification/NapA/seed_refpkg/final_outputs \
	--fastx_input /mnt/datasets/2021w/saanich/MAGs/Concatenated/All_SI072_Metawrap_MAGs.fa \
	--output SI072_MAGs_assign/


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

