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

# #Same error output after running treesapp update for MAGs with and without --trim_align.
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
