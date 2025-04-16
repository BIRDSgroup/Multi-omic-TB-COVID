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
load("data/family_phseq.RData")
load("Figures/DA_analysis/TB_TBCOVID_family_type/top_taxa.RData")

mat_otu=fread("picrust_analysis/pathways_out_final/path_abun_unstrat_descrip_perseq.tsv",header=TRUE,sep="\t")
#superpathway of pyridoxal 5'-phosphate biosynthesis and salvage (PWY0-845) to be replaced without ' in the above file


top_level<-read.csv("data/metacyc_pathways_info_prokaryotes_top_level.tsv",header = FALSE,sep = "\t")

physeq_filtered <- subset_samples(physeq, X.type %in% c("TB", "TBCOVID"))
mat_otu_final<-otu_table(physeq_filtered)
mat_otu_final<-log(mat_otu_final[top_taxa,]+1)

mat_otu$description<-paste0(mat_otu$description," (",mat_otu$pathway,")")

indices<-match(rownames(mat_otu_final),mat_otu$description)

ind_top<-match(mat_otu$pathway[indices],top_level$V1)

colnames(mat_otu_final)<-sample_data(physeq_filtered)$X.OTU
rownames(mat_otu_final)<-paste0(top_level$V2[ind_top]," (",top_level$V1[ind_top],")")

colors <- c("#ADFF2F", "#4B0082")
unique_ids <- unique(sample_data(physeq_filtered)$X.type)


color_map <- setNames(colors, unique_ids)
row_colors <- color_map[sample_data(physeq_filtered)$X.type]

status<-rep("cured/lost_to_followup",48)
status[c(10,11,15,19,21,28,29,30,34,42,44)]<-"relapse/DR/died"
anno_col<-data.table(Group = sample_data(physeq_filtered)$X.type,Status=status)
anno_colors<-list(Group = color_map)
rownames(anno_col)<-colnames(mat_otu_final)


class(mat_otu_final)<-"matrix"

status=read.csv("data/metadata_status.csv",check.names = FALSE)

generate_colors <- function(labels, palette = "Set1") {
  colorRampPalette(brewer.pal(8, palette))(length(unique(labels)))
}


ann_colors <- list(Group = c("TB" = "#ADFF2F", "TBCOVID" = "#4B0082")  )
dists <- dist(t(mat_otu_final), method = "manhattan")
o = seriate(dists, method = "OLO")
column_h<-HeatmapAnnotation(Group = sample_data(physeq_filtered)$X.type,col=ann_colors,show_legend = TRUE)

ann_colors <- list(
  `Single/Multiple episodes` = setNames(generate_colors(status$`Single/Multiple episodes`), unique(status$`Single/Multiple episodes`)),
  `Treatment` = setNames(generate_colors(status$Treatment, "Set2"), unique(status$Treatment)),
  `Drug Resistance`=setNames(generate_colors(status$`Drug Resistance`,"Set1"),unique(status$`Drug Resistance`)))

top_h<-HeatmapAnnotation( "Single/Multiple episodes"=status$`Single/Multiple episodes`,"Treatment status"=status$Treatment,"Drug resistance status" = status$`Drug Resistance`,col=ann_colors,
                          show_legend = TRUE)
status$`Drug Resistance`<-factor(trimws(status$`Drug Resistance`))
h<-Heatmap(mat_otu_final,top_annotation=top_h,bottom_annotation = column_h,column_order = get_order(o),col=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
           heatmap_legend_param = list(title = "log(Pathway abundance)",direction="horizontal"),show_heatmap_legend = TRUE,
           heatmap_width = unit(24, "cm"), 
           heatmap_height = unit(18, "cm"))
h
png("Figures/manuscript_results/pathway_heatmap.png", width=5000,height=2700,res = 300)
draw(h, annotation_legend_side = "bottom",heatmap_legend_side="top")
dev.off()

############# Fig S5 ######################################

load("data/species_phseq.RData")
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


