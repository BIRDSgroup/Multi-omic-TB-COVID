source("utils.R")

#generate phyloseq objects at genus and species levels
gen_physeq("genus")
gen_physeq("species")

load('../data/repeatsamples_genus_phseq.RData')

counts<-mat_otu_final
mat_otu_final<-t(mat_otu_final) + 1
mat_otu_final<-logratio.transfo(mat_otu_final,logratio = 'CLR')
mat_otu_final<-t(mat_otu_final)
class(mat_otu_final)<-"matrix"

##########   Fig 2a ########################################

plot_scatter<-function(id,id1,id2){
  mat_df<-counts[,c(id,id1,id2)]
  mat_df<-mat_df[!apply(mat_df, 1, function(x) all(x == 0)), ]
  mat_df<-t(logratio.transfo(t(mat_df+1),logratio = 'CLR'))
  class(mat_df)<-"matrix"
  mat_df<-as.data.frame(mat_df)
  df.gathered <- mat_df %>%
    as_data_frame() %>%
    gather(key = "Sample id", value = "value",
           id1,id2)
  return(df.gathered)
}

df<-plot_scatter("SM007","SM039","SM088")
df<-ggplot(df, aes(x = SM007, y = value)) +
  ggplot2::geom_point(
    ggplot2::aes(fill = `Sample id`),
    alpha = 1 ,
    size = 3,
    shape = 21,
    stroke = 0.35,
    position = ggplot2::position_jitterdodge()
  ) +
  facet_wrap(~ `Sample id`, scales = "free") +
  ggplot2::scale_fill_brewer(palette = "Spectral")+
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top")

temp_plot1<- df + 
  nature_theme("SM007", "clr(relative abundance)") +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(colour = "black")
  ) +
  ggplot2::xlab("SM007") +
  ggplot2::ylab("Genus read counts (CLR-transformed)") +
  ggplot2::theme(legend.position = "none") 
temp_plot1
ggsave("../Figures/manuscript_results/Sample_scatter.png",temp_plot1,width=5,height=3,units="in")

##########   Fig 2b  ########################################

subset<-counts[,repeatsamples$X.OTU]

meta_sub<-metadata[metadata$X.OTU %in% repeatsamples$X.OTU,]
combined_names=paste(meta_sub$X.OTU, meta_sub$X.id, sep = "_")


unique_ids <- unique(meta_sub$X.id)
library(RColorBrewer)

colors <- c("#FF6347", "#FFD700", "#FF8C00", "#ADFF2F", "#8B0000", 
            "#4B0082", "#20B2AA", "#FFFF00")

color_map <- setNames(colors, unique_ids)
row_colors <- color_map[meta_sub$X.id]

colnames(subset)<-combined_names
cor_matrix <- cor(subset, method = "spearman")
anno_row<-data.frame(Group = meta_sub$X.id)
anno_colors<-list(Group = color_map)
rownames(anno_row)<-rownames(cor_matrix)

p<-pheatmap(cor_matrix,show_colnames = F,cluster_cols =F,clustering_distance_rows="correlation",fontsize = 10,annotation_row = anno_row,annotation_colors = anno_colors, annotation_legend = FALSE)

ggsave("../Figures/manuscript_results/corrplot.png",p,width=7,height=5,units="in")

##########   Fig 2c ########################################

load("../data/species_phseq.RData")
#Loading taxonomy data
mat_otu_final<-as.data.frame(otu_table(physeq))
cal_clr<-function(mat_otu_final){
  total_counts_per_sample <- rowSums(mat_otu_final)
  mat_otu_final<-sweep(mat_otu_final, 1, total_counts_per_sample, "/")
  mat_otu_final<-mat_otu_final*100
  
  mat_otu_final<-logratio.transfo(mat_otu_final,logratio = 'CLR')
  return(mat_otu_final)
}

mat_otu1<-t(mat_otu_final+1)
mat_otu1<-t(cal_clr(mat_otu1))

metadata=read.csv("../data/metadata_R.csv",header=TRUE,sep=",")

colnames(mat_otu1)<-sample_data(physeq)$X.OTU

#find the row indices of the otu table corresponding to Tuberculosis species
row_ix1 <- which(grepl("tuberculosis", rownames(otu_table(physeq))))

tuberculosis_df<- data.frame(matrix(ncol = 82, nrow = 3))
tuberculosis_df[1,]<-as.numeric(mat_otu1[row_ix1,])
tuberculosis_df[2,]<-sample_data(physeq)$X.type
tuberculosis_df[3,]<-sample_data(physeq)$X.OTU
tuberculosis_df<-as.data.frame(t(tuberculosis_df))
colnames(tuberculosis_df)<-c("CLR-transformed counts","Type","Sample Number")
tuberculosis_df$`CLR-transformed counts`<-as.numeric(tuberculosis_df$`CLR-transformed counts`)
p<-ggplot(tuberculosis_df, aes(x = `Sample Number`, y = Type)) +
  geom_tile()+
  facet_wrap(Type ~., ncol=1,scales="free",labeller = label_wrap_gen(multi_line = TRUE)) +
  geom_raster(aes(fill = `CLR-transformed counts`), interpolate=FALSE) +
  ylab("")+
  scale_fill_gradient2(low = "blue", mid="white", high = "red",midpoint=0) +
  theme_minimal()+
  theme(aspect.ratio = 0.2,axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text = element_text(color = "black",size=20), axis.title = element_text(color = "black",size=20),
        legend.title = element_text(color = "black",size=20),strip.text = element_text(color = "black",size=20),
        legend.text=element_text(color = "black",size=20),legend.position = "bottom",legend.direction = "horizontal",
        axis.text.y = element_blank())
p 
ggsave("../Figures/manuscript_results/heatmap.png",plot = p, dpi = 300 ,width=16,height=16)


##########   Fig 2d ########################################


load('../data/species_phseq.RData')
physeq<-subset_samples(physeq, X.type %in% c("TB","TBCOVID"))
mat_otu_final<-as.data.frame(otu_table(physeq))
class(mat_otu_final)<-"data.frame"

input<-t(mat_otu_final)
input=input+0.01
input=logratio.transfo(input,logratio = 'CLR')

metadata<-as.data.frame(sample_data(physeq))
class(metadata)<-"data.frame"

for(j in 32:32)
{
  input_df<-data.frame(input[,j])
  colnames(input_df)<-colnames(input)[j]
  input_df$y<-input[,j]
  x_label="Smear"
  input_df$x<-as.character(metadata[[x_label]])
  y_label=colnames(input)[j]
  x_axis_label_names <- unique(input_df[['x']])
  renamed_levels <- as.character(levels(metadata[,x_label]))
  if (length(renamed_levels) == 0) {
    renamed_levels <- x_axis_label_names
  }
  
  for (name in x_axis_label_names) {
    total <- length(which(input_df[['x']] == name))
    new_n <- paste(name, " (n=", total, ")", sep="")
    input_df[c(which(input_df$x == name)),'x'] <- new_n
    renamed_levels<- replace(renamed_levels, renamed_levels == name, new_n)
  }
  input_df$xnames <- factor(input_df[['x']], levels=renamed_levels)
  
  temp_plot <-
    ggplot2::ggplot(
      data = input_df, ggplot2::aes(xnames, y)) +
    ggplot2::geom_boxplot(
      ggplot2::aes(fill = x),
      outlier.alpha = 0.0,
      na.rm = TRUE,
      alpha = .5,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = x),
      alpha = 0.75 ,
      size = 1,
      shape = 21,
      stroke = 0.15,
      color = 'black',
      position = ggplot2::position_jitterdodge()
    ) +
    ggplot2::scale_fill_brewer(palette = "Spectral")
  
  temp_plot1<- temp_plot + 
    nature_theme(input_df[, 'x'], y_label) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black")
    ) +
    ggplot2::xlab("Severity") +
    ggplot2::ylab('Mycobacterium tuberculosis (CLR-transformed)') +
    ggplot2::theme(legend.position = "none") 
  
  #temp_plot1
  png_file<-paste0("../Figures/manuscript_results/TB Severity.png")
  png(png_file, res = 300, width = 600, height = 900)
  stdout <- capture.output(print(temp_plot1))
  dev.off()
}

input_df["Group"]<-as.character(metadata[["X.type"]])

temp_plot <-
  ggplot2::ggplot(
    data = input_df, ggplot2::aes(xnames, y)) +
  ggplot2::geom_boxplot(
    ggplot2::aes(fill = Group),
    outlier.alpha = 0.0,
    na.rm = TRUE,
    alpha = .5,
    show.legend = TRUE
  ) +
  ggplot2::geom_point(
    ggplot2::aes(fill = Group),
    alpha = 0.75 ,
    size = 1,
    shape = 21,
    stroke = 0.15,
    color = 'black',
    position = ggplot2::position_jitterdodge()
  ) + scale_fill_manual(values=c("#E69F00", "#00BFC4")) 

temp_plot1<- temp_plot + 
  nature_theme(input_df[, 'x'], y_label) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(colour = "black")
  ) +
  ggplot2::xlab("Severity") +
  ggplot2::ylab('Mycobacterium tuberculosis (CLR-transformed)') +
  ggplot2::theme(legend.position = "bottom") 

png_file<-paste0("../Figures/manuscript_results/TB Severity_group.png")
png(png_file, res = 300, width = 600, height = 1100)
stdout <- capture.output(print(temp_plot1))
dev.off()

##########   Fig S5 ########################################
batch_effect_analysis("genus")
batch_effect_analysis("species")
