# edit tip names ##
sed 's/|.*//g' input.fasta > input_unaligned.fasta ;
grep ">" input_unaligned.fasta ; 
grep ">" input_unaligned.fasta | wc -l ;
ls -lh ; 

# run nextclade #
nextclade dataset get --name hMPXV --output-dir DB
nextclade run --input-dataset DB/ --output-all out input_unaligned.fasta
mv out/nextclade.aligned.fasta out/nextclade.tsv .
mv nextclade.aligned.fasta sequences_total.fasta

## extract selected sequences ##
seqtk subseq sequences_total.fasta ext.txt > sequences.fasta
seqtk subseq sequences_total.fasta list2.txt > sequences_1.fasta
seqtk subseq snps2.fasta list2.txt > snps3.fasta

## run snp-sites ## 
snp-sites -o snps.fasta sequences.fasta ; 
aliview snps.fasta ; 

# run raxml-nextsrain #

raxmlHPC-PTHREADS -p 123568 -m GTRCAT -s snps.phy -T 31 -# 10 -n nwk ; 
mv RAxML_bestTree.nwk raw_tree.nwk ;
mkdir raxml ; 
mv RAx* raxml/ ; 
conda activate nextstrain ; 
augur refine --alignment sequences.fasta --tree raw_tree.nwk --metadata metadata.tsv --output-tree refine_tree.nwk --output-node-data node_Data.json --timetree --coalescent opt --gen-per-year 52 --root least-squares --covariance --date-confidence --date-inference marginal --branch-length-inference marginal --year-bounds 2017 2022 --divergence-units mutations-per-site --seed 12548746 ; 
augur ancestral --tree refine_tree.nwk --alignment sequences.fasta --output-node-data ancestral.json --inference marginal --keep-overhangs ; 
augur translate --tree refine_tree.nwk --ancestral-sequences ancestral.json --reference-sequence ref_monkey.gb --output-node-data translate.json ; 
augur traits --tree refine_tree.nwk --metadata metadata.tsv --columns continent country province lineage age sex clade --confidence --output-node-data traits.json ; 
AUGUR_RECURSION_LIMIT=10000 augur export v2 --tree refine_tree.nwk --node-data ancestral.json node_Data.json traits.json translate.json --output final.json --auspice-config auspice_config.json --geo-resolutions country --color-by-metadata continent country province lineage age sex clade --panels tree map entropy frequencies --metadata metadata.tsv --lat-longs lat_longs.tsv ; 
ls -lh ; 

# only phyplogeny # 
raxmlHPC-PTHREADS -p 123568 -m GTRCAT -s snps.phy -T 31 -# 10 -n nwk ; 
mv RAxML_bestTree.nwk raw_tree.nwk ;
mkdir raxml ; 
mv RAx* raxml/ ; 
ls -lh ; 


