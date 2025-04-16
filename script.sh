#!/bin/bash

#load("Figures/DA_analysis/TB_TBCOVID_family_type/TB_TBCOVID_family_saved_plot.RData")
#write.csv(top_taxa, "picrust_analysis/TBTBCOVID/TBTBCOVID.csv", row.names = FALSE)
while read line; do  

    if grep -qw "$line" ./picrust_analysis/pathways_out_final/path_abun_unstrat_descrip_perseq.tsv; then
        pwy=$(grep -w "$line" ./picrust_analysis/pathways_out_final/path_abun_unstrat_descrip_perseq.tsv| cut -f 1)

        grep $pwy ./picrust_analysis/pathways_out_final/path_abun_contrib.tsv > aa.txt

        grep -Fw -f  ./picrust_analysis/TBTBCOVID/TBTBCOVIDids.csv aa.txt | cut -f 3 | sort > extracted_seqids.txt
        
        while read pattern; do   
        grep -wF "$pattern" ./picrust_analysis/export_taxa_picrust/taxonomy_picrust.tsv;
        done < extracted_seqids.txt | cut -f 2 |sort |  uniq -c > temp.txt
        sed 's/^[ \t]*//' temp.txt > "./picrust_analysis/TBTBCOVID/TBTBCOVID_${line}.txt"
    fi
    rm temp.txt
    rm aa.txt
    rm extracted_seqids.txt

done < ./picrust_analysis/TBTBCOVID/TBTBCOVID.csv