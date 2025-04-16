library(data.table)
load(paste0("data/species_phseq.RData"))
mat_otu=fread("picrust_analysis/pathways_out_final/path_abun_unstrat_descrip_perseq.tsv",header=TRUE,sep="\t")
mat_otu<-as.data.frame(mat_otu)

rownames(mat_otu)<-paste0(mat_otu$description," (",mat_otu$pathway,")") #for creating Fig S7 so that y labels include pathway name as well description

sample_ids<-intersect(sample_data(physeq)$X.OTU, colnames(mat_otu))
mat_otu<-mat_otu[,sample_ids]
mat_tt<-data.frame("family"=rownames(mat_otu))

tt=tax_table(mat_tt)
rownames(tt)=mat_tt$family
colnames(tt)="family"
physeq_sam=sample_data(physeq)
physeq_taxtable=phyloseq::tax_table(tt)
mat_otu_final<-as.matrix(mat_otu)
colnames(mat_otu_final) <-sample_names(physeq_sam)
rownames(mat_otu_final)<-rownames(tt)
physeq_otu=otu_table(round(mat_otu_final),TRUE)
physeq = phyloseq(physeq_otu, physeq_taxtable,physeq_sam)

save(physeq,file="data/family_phseq.RData")    



