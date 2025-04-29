library(data.table)
library(mixOmics)
library(circlize)
library(pheatmap)
library(seriation)
library(ComplexHeatmap)
library(stringr)
library(RColorBrewer)
cl_cb <- function(hcl, mat){
  # Recalculate manhattan distances for reorder method
  dists <- dist(mat, method = "manhattan")
  
  # Perform reordering according to OLO method
  hclust_olo <- reorder(hcl, dists)
  return(hclust_olo)
}


############# Fig 8 ######################################
load("../data/family_phseq.RData")
load("../Figures/DA_analysis/TB_TBCOVID_family_type/top_taxa.RData")

mat_otu=fread("../picrust_analysis/pathways_out_final/path_abun_unstrat_descrip_perseq.tsv",header=TRUE,sep="\t")

mat_otu$description <- gsub("pyridoxal 5'", "pyridoxal 5", mat_otu$description)
#superpathway of pyridoxal 5'-phosphate biosynthesis and salvage (PWY0-845) to be replaced without ' in the above file


top_level<-read.csv("../data/metacyc_pathways_info_prokaryotes_top_level.tsv",header = FALSE,sep = "\t")

physeq_filtered <- subset_samples(physeq, X.type %in% c("TB", "TBCOVID"))
mat_otu_final<-otu_table(physeq_filtered)
#mat_otu_final<-log(mat_otu_final[top_taxa,]+1)
mat_otu_final=logratio.transfo(t(mat_otu_final+1),logratio = 'CLR')
mat_otu_final<-t(mat_otu_final)
mat_otu_final<-mat_otu_final[top_taxa,]

mat_otu$description<-paste0(mat_otu$description," (",mat_otu$pathway,")")

indices<-match(rownames(mat_otu_final),mat_otu$description)

ind_top<-match(mat_otu$pathway[indices],top_level$V1)

colnames(mat_otu_final)<-sample_data(physeq_filtered)$X.OTU
descp<-top_level$V2[ind_top]
rownames(mat_otu_final)<-top_level$V1[ind_top]


class(mat_otu_final)<-"matrix"

status=read.csv("../data/metadata_status.csv",check.names = FALSE)



status$Challenging[status$Challenging==0]="No"
status$Challenging[status$Challenging==1]="Yes"

dists <- dist(t(mat_otu_final), method = "manhattan")
o = seriate(dists, method = "OLO")
row_h<-rowAnnotation(`Pathway Type`= descp,show_legend = TRUE,
                     col=list(`Pathway Type` =setNames(brewer.pal(5, "Set3"),unique(descp))))

ann_colors <- list(
  `Single/Multiple episodes` = setNames(c("#C04BDB","#FFF2CC"), unique(status$`Single/Multiple episodes`)),
  `Treatment status` = setNames(c("#6495ED","#FFD700","#D5006D","#32CD32"), unique(status$Treatment)),
  `Drug resistance status`=setNames(c("red","#585858"),unique(status$`Drug Resistance`)),
  `Group`=setNames(c("#4B0082","#ADFF2F"),unique(sample_data(physeq_filtered)$X.type)),
  "Adverse outcome"=setNames(c("#F55F74","#9E7F2A"),unique(as.factor(status$Challenging))))

top_h<-HeatmapAnnotation( "Single/Multiple episodes"=status$`Single/Multiple episodes`,`Treatment status`=status$Treatment,`Drug resistance status` = status$`Drug Resistance`,"Group" = sample_data(physeq_filtered)$X.type,col=ann_colors,
                          show_legend = TRUE)
#bottom_h<-HeatmapAnnotation("Adverse outcome"=as.factor(status$Challenging),col=ann_colors)
bottom_h<-HeatmapAnnotation("Adverse outcome"=as.factor(status$Challenging),"Gender"=sample_data(physeq_filtered)$Gender,col=ann_colors)

status$`Drug Resistance`<-factor(trimws(status$`Drug Resistance`))
#col=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
h<-Heatmap(mat_otu_final,bottom_annotation=bottom_h,cluster_columns=hclust(dist(t(mat_otu_final)), method = "ward.D2"),left_annotation=row_h,top_annotation=top_h,col=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
           heatmap_legend_param = list(title = "clr(Pathway abundance)",direction="horizontal"),show_heatmap_legend = TRUE,
           heatmap_width = unit(24, "cm"), 
           heatmap_height = unit(18, "cm"),cluster_rows=TRUE,column_split = 2,
           show_column_dend = TRUE)

                          
                          
png("../Figures/manuscript_results/pathogen_heatmap_1_gender.png", width=5000,height=2900,res = 300)
draw(h, annotation_legend_side = "bottom",heatmap_legend_side="top")
dev.off()


######################################################

load("../data/species_phseq.RData")
status=read.csv("../data/metadata_status.csv",check.names = FALSE)

sp<-c("Acinetobacter baumannii", "Klebsiella pneumoniae", "Haemophilus influenzae", "Stenotrophomonas maltophilia", "Alloscardovia omnicolens")
physeq_filtered <- subset_samples(physeq, X.type %in% c("TB", "TBCOVID"))
mat_otu_final_species<-otu_table(physeq_filtered)
mat_otu_final_species=logratio.transfo(t(mat_otu_final_species+1),logratio = 'CLR')
mat_otu_final_species<-t(mat_otu_final_species)
mat_otu_final_species<-(mat_otu_final_species[sp,])
colnames(mat_otu_final_species)<-sample_data(physeq_filtered)$X.type



h<-Heatmap(mat_otu_final_species,cluster_columns=hclust(dist(t(mat_otu_final)),method = "ward.D2"),
           heatmap_legend_param = list(title = "clr(Species abundance)",direction="horizontal"),show_heatmap_legend = TRUE,
           heatmap_width = unit(24, "cm"), 
           heatmap_height = unit(8, "cm"),show_column_dend = FALSE,column_split = 2)
h

png("../Figures/manuscript_results/heatmap_2_gender.png", width=5000,height=2700,res = 300)
draw(h, annotation_legend_side = "bottom",heatmap_legend_side="top")
dev.off()


sp<-c("Prevotella melaninogenica","Capnocytophaga gingivalis","Veillonella parvula","Escherichia coli")

mat_otu_final_species<-(mat_otu_final_species[sp,])

class(mat_otu_final_species) <- "matrix"
top_h<-HeatmapAnnotation( "Gender"=sample_data(physeq_filtered)$Gender, show_legend = TRUE)
h<-Heatmap(mat_otu_final_species,cluster_columns=hclust(dist(t(mat_otu_final_species)),method = "ward.D2"),top_annotation=top_h,col=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
           heatmap_legend_param = list(title = "clr(Species abundance)",direction="horizontal"),show_heatmap_legend = TRUE,
           heatmap_width = unit(24, "cm"), 
           heatmap_height = unit(8, "cm"),show_column_dend = FALSE,column_split = 2, show_row_names = FALSE)




ordered_indices <- order(sample_data(physeq_filtered)$X.type)
matrix_data_ordered <- mat_otu_final_species[, ordered_indices]
group<-sample_data(physeq_filtered)$X.type
group_ordered <- group[ordered_indices]

# Then plot

top_h<-HeatmapAnnotation( "Gender"=sample_data(physeq_filtered)$Gender[ordered_indices], show_legend = TRUE)
Heatmap(matrix_data_ordered,top_annotation=top_h,
        column_split = group_ordered,  # Split columns by group
        cluster_columns = FALSE,        # Do not cluster columns within groups
        show_column_names = TRUE,heatmap_width = unit(24, "cm"), 
        heatmap_height = unit(8, "cm"),)


h<-Heatmap(mat_otu_final_species,cluster_columns=hclust(dist(t(mat_otu_final_species)),method = "ward.D2"),top_annotation=top_h,col=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
           heatmap_legend_param = list(title = "clr(Species abundance)",direction="horizontal"),show_heatmap_legend = TRUE,
           heatmap_width = unit(24, "cm"), 
           heatmap_height = unit(8, "cm"),show_column_dend = FALSE,column_split = 2, show_row_names = TRUE)

h


png("../Figures/manuscript_results/gender_heatmap2.png", width=5000,height=2700,res = 300)
draw(h, annotation_legend_side = "bottom",heatmap_legend_side="top")
dev.off()


col_clusters<-column_order(h)
cluster_identity<-rep("Cluster 2",48)

cluster_identity[c(col_clusters[[1]])]<-"Cluster 1"


p_values <- apply(mat_otu_final_species, 1, function(row) {
  # Perform t-test between the two groups
  t_test_result <- t.test(row[cluster_identity == "Cluster 1"], row[cluster_identity == "Cluster 2"])
  return(t_test_result$p.value)  # Return the p-value from the t-test
})

fdr_corrected_p_values <- p.adjust(p_values, method = "BH")

############# Fig S5 ######################################

load("../data/species_phseq.RData")
sp<-c("Acinetobacter baumannii", "Klebsiella pneumoniae", "Haemophilus influenzae", "Stenotrophomonas maltophilia", "Alloscardovia omnicolens")
pathogen_list<-c()
count_list<-c()
group_list<-c()
dtype<-list("TB","TBCOVID","Control","COVID")
for(pathogen in sp)
{
  for(typ in 1:4)
  {
    physeq_filtered <- subset_samples(physeq, X.type %in% dtype[typ])
    mat_otu_final<-otu_table(physeq_filtered)
    mat_otu_final<-(mat_otu_final[pathogen,])
    mat_otu_final[mat_otu_final<=5]=0
    count_list<-c(count_list,sum( mat_otu_final!=0))
    pathogen_list<-c(pathogen_list,pathogen)
    group_list<-c(group_list,unlist(dtype[typ]))
  }
}
pathogen_data<-data.frame(Pathogen=pathogen_list,count=count_list,Group=group_list)
pgplot<-ggplot(pathogen_data, aes(x = Pathogen, y=count,fill = Group)) +
  geom_col(position = "stack") +
  scale_fill_viridis_d(option = "plasma") +
  labs(x = "Pathogen", y = "Number of Individuals") +
  ggplot2::theme_bw()+
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(size = 10, vjust = 1, hjust = 1, angle = 45),
    axis.text.y = ggplot2::element_text(size = 10, hjust = 1),
    axis.title = ggplot2::element_text(size = 10),
    legend.title = ggplot2::element_text(size = 10, face = 'bold'),
    legend.text = ggplot2::element_text(size = 10),
    axis.line = ggplot2::element_line(colour = 'black', size = .25),
    axis.line.x = ggplot2::element_line(colour = 'black', size = .25),
    axis.line.y = ggplot2::element_line(colour = 'black', size = .25),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    legend.position = "top"
  ) 

png_file<-paste0("Figures/manuscript_results/pathogen.png")
png(png_file, res = 300, width = 2000, height = 2000)
stdout <- capture.output(print(pgplot))
dev.off()

