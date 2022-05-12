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
