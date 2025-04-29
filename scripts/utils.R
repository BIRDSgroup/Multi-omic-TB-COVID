library('dplyr')
library('tidyr')
library('ggplot2')
library(ggpubr)
library(mixOmics) #log ratio transformation
library("vegan")
library(pheatmap)
library(viridis)
library(phyloseq) 
library(vegan) 
library(reshape2) 
library(data.table)
library(dplyr)
library(ggnewscale) 
library(ggsignif)

######## Alpha diversity #############
plot_alphadiversity<-function()
{
  richness_values <- estimate_richness(physeq,measures = c("Observed", "Shannon","Simpson"))
  richness_values["SampleType"]<-sample_data(physeq)$X.type
  p_val<-c(0,0,0)
  anova_result <- aov(Observed ~ SampleType, data = richness_values)
  p_val[1] <- summary(anova_result)[[1]]$`Pr(>F)`[1]
  anova_result <- aov(Shannon ~ SampleType, data = richness_values)
  p_val[2]<- summary(anova_result)[[1]]$`Pr(>F)`[1]
  anova_result <- aov(Simpson ~ SampleType, data = richness_values)
  p_val[3] <- summary(anova_result)[[1]]$`Pr(>F)`[1]
  
  a_d<-physeq %>%                                                              #phyloseq object
    plot_richness(x = "X.type", measures = c("Observed", "Shannon","Simpson"), nrow = 4)+                            
    geom_boxplot(aes(fill = X.type), show.legend = FALSE,alpha=0.75,outliers = FALSE,outlier.shape = NA)+  ggplot2::geom_point(
      aes(fill = X.type),
      alpha = 08 ,
      size = 1,
      shape = 21,
      stroke = 0.15,
      color = 'black',
      position = ggplot2::position_jitterdodge()
    ) +
    stat_anova_test(method="one_way",step.increase = 2,vjust=1,label = "ANOVA, p-value = {p}")+
    theme_linedraw()+                                                     #change theme to classic
    xlab(NULL)+                                                           #no label on x-axis
    theme(axis.text.y.left = element_text(size = 12),                     #adjust y-axis text
          axis.text.x = element_text(size = 12, hjust = 1, vjust = 1, angle = 45),           #adjust x-axis label position
          axis.title.y = element_text(size = 12))+                        #adjust y-axis title
    theme(strip.text = element_text(face = "bold", size = 10))+           #adjust headings
    scale_fill_manual(values=c("Control"="#A1C400","COVID"="#B03060","TB"="#E69F00","TBCOVID"= "#00BFC4")) +                          #set fill colors
    theme(plot.title=element_text(size = 11, face = "bold", hjust = 0.5),legend.position = "bottom") #change title size, face and position
  a_d$layers <- a_d$layers[-1]
  return(a_d)
}

######### PCoA plot ###################
plot_pcoa <- function() {
  
  ps_subset <- subset_samples(physeq, sample_data(physeq)$X.type %in% pair)
  metadata<-sample_data(ps_subset)
  class(metadata)="data.frame"
  dist_matrix <- vegdist(t(otu_table(ps_subset)), method = "bray")
  pcoa_result <- cmdscale(dist_matrix, k = 3, eig = TRUE)
  
  pca_com<-as.data.frame(pcoa_result$points)
  
  density.lwd = 0.2
  title.cex = 1.5
  legend.cex = 0.7
  legend.title.cex = 0.75
  color.set <- color.mixo(seq_len(10))
  
  Type<-as.factor(metadata$X.type)
  minx<-min(pca_com$V1)
  maxx<-max(pca_com$V1)
  miny<-min(pca_com$V2)
  maxy<-max(pca_com$V2)
  
  ######Permanova
  if(pair[1]=="Control")
  {
    formula=as.formula("dist_matrix~X.type+X.batch_no")
  }
  else
  {
    formula=as.formula("dist_matrix~X.type+X.batch_no+Diabetes+TOBACCO.Smoking+Alcohol+Smear+Age+Height+Weight+Gender")
  }
  set.seed(31125) 
  permanova_results_bray <- adonis2(formula, data=metadata,permutations = 999,by="margin")
  bs<-betadisper(dist_matrix,metadata$X.type)
  
  bs_res<-anova(bs)
  
  print(bs_res)
  
  #print(permanova_results_bray)
  
  print(paste('R2=',round(permanova_results_bray$R2[1],3),', P=',permanova_results_bray$`Pr(>F)`[1],sep=''))
  
  
  return(ggplot(data = pca_com, aes(x = V1, y = V2,
                                    color = Type)) +
           geom_point(size=2,alpha=0.8) + xlab(paste0('PC1: ', round(as.numeric(pcoa_result$eig[1]/sum(pcoa_result$eig)*100)),'% expl.var')) +
           ylab(paste0('PC2: ', round(as.numeric(pcoa_result$eig[2]/sum(pcoa_result$eig)*100)),'% expl.var')) +
           scale_shape_manual(values = c(1,19,2,17,4)) +
           scale_color_manual(values = color.set) + theme_bw() +
           scale_x_continuous(limits = c(minx-0.5,maxx+0.5)) + scale_y_continuous(limits = c(miny-0.5,maxy+0.5)) +
           scale_color_manual(labels=pair, values=c("#E69F00", "#00BFC4")) +
           stat_ellipse()+
           theme(legend.position = c(0.8,0.2), legend.box = 'horizontal',
                 legend.direction = 'vertical',
                 legend.key.height = unit(0.2, 'cm'),
                 legend.key.width = unit(0.1, 'cm'),
                 legend.title = element_text(size =12),
                 legend.spacing.x = unit(0.1, 'cm'),
                 legend.spacing.y = unit(0.1, 'cm'),
                 legend.text = element_text(size = 12),
                 axis.text.x = ggplot2::element_text(size = 10, vjust = 1, hjust = 1),
                 axis.text.y = ggplot2::element_text(size = 10, hjust = 1),))
}


############# top taxa  ##########################

plot_toptaxa<-function(t_level)
{
  physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))
  
  # Aggregate by species and calculate the mean relative abundance for each group
  species_abundance <- physeq_rel %>%
    tax_glom(taxrank = t_level) %>%                      
    psmelt() %>%                                           
    group_by(across(all_of(t_level)), X.type) %>%                           
    summarize(MeanAbundance = mean(Abundance, na.rm = TRUE)) %>%
    ungroup()
  
  # Get the top 10 species based on overall mean abundance
  top_species <- species_abundance %>%
    group_by(.data[[t_level]]) %>%
    summarize(TotalMeanAbundance = mean(MeanAbundance)) %>%
    arrange(desc(TotalMeanAbundance)) %>%
    slice_head(n = 10) %>%
    pull(.data[[t_level]])
  
  # Filter the original data for the top 10 species
  filtered_data <- species_abundance %>%
    filter(.data[[t_level]] %in% top_species)
  
  
  filtered_data 
  sorted_data<-arrange(filtered_data,desc(MeanAbundance))
  
  # Create the bar plot
  temp<-ggplot(sorted_data, aes(x = reorder(.data[[t_level]], MeanAbundance), y = MeanAbundance, fill = X.type)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(x = NULL, y = "Mean Relative Abundance", fill = "Type") +
    scale_fill_manual(values=c("Control"="#A1C400","COVID"="#B03060","TB"="#E69F00","TBCOVID"= "#00BFC4")) +
    coord_flip() + # Flip coordinates for a horizontal bar plot
    theme_bw() +
    theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10),axis.text.y = element_text(angle = 0, hjust = 1,size = 10),legend.position = "bottom")
  return(temp)
}


############ plots theme #########################
nature_theme <- function(x_axis_labels, y_label) {
  # set default text format based on categorical and length
  angle = NULL
  hjust = NULL
  size = 8
  if (max(nchar(x_axis_labels), na.rm=TRUE) > 5) {
    angle = 45
    hjust = 1
    size = 6
  }
  axis_title_size = 10
  if (nchar(y_label) > 15) {
    axis_title_size = 8
  }
  if (nchar(y_label) > 25) {
    axis_title_size = 6
  }
  return ( ggplot2::theme_bw() + ggplot2::theme(
    axis.text.x = ggplot2::element_text(size = size, vjust = 1, hjust = hjust, angle = angle),
    axis.text.y = ggplot2::element_text(size = 8, hjust = 1),
    axis.title = ggplot2::element_text(size = axis_title_size),
    plot.title = ggplot2::element_text(size = 7, face = 'bold'),
    legend.title = ggplot2::element_text(size = 6, face = 'bold'),
    legend.text = ggplot2::element_text(size = 6),
    axis.line = ggplot2::element_line(colour = 'black', size = .25),
    axis.line.x = ggplot2::element_line(colour = 'black', size = .25),
    axis.line.y = ggplot2::element_line(colour = 'black', size = .25),
  )
  )
}

############ generate phyloseq object  #########################
gen_physeq<-function(f_level)
{
    mat_s<-read.table("data/merged_homd_tax.txt",header=TRUE,sep=",")  
    mat_otu=read.table("data/merged_feature_table.tsv",header=FALSE,sep="\t")
    mat_otu_df<-as.data.frame(mat_otu)

    # Split taxonomy string into components
    tax_components <- strsplit(mat_otu_df$V1, ";")

    # Extract taxonomic levels into separate columns
    df_split <- as.data.frame(matrix(unlist(sapply(tax_components, function(x) sub("^[a-z]+__", "", x))), ncol = 7, byrow = TRUE))

    # Rename columns
    colnames(df_split) <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

    mat_otu_df<-cbind(df_split,mat_otu_df)
    mat_otu_df<-mat_otu_df[,-which(names(mat_otu_df) == "V1")]
    metadata=read.csv("data/metadata_R.csv",header=TRUE,sep=",")



    colnames(mat_otu_df)[8:dim(mat_otu_df)[2]]<-c(as.vector(((metadata[c(1:dim(metadata)[1]),c(1:1)]))))

    mat_otu_df_summarize<-function(f_level,mat_otu_df_copy){
      mat_otu_df_reduce<-mat_otu_df_copy[(mat_otu_df_copy[f_level]!="__")[1:dim(mat_otu)[1],1],]
      if(f_level!="species")
      {
      summary_df <- mat_otu_df_reduce %>%
        group_by(!!sym(f_level)) %>%
        summarise(across(.cols = SM001:SM093, .fns = sum, na.rm = TRUE))
      }
      else
        {summary_df<-mat_otu_df_reduce}
      return(summary_df)
    }

    summary_df<-mat_otu_df_summarize(f_level,mat_otu_df)
    
    summary_df_species<-mat_otu_df_summarize("species",mat_otu_df)

    if(f_level!="species"){
      mat_otu_final=summary_df[c(1:dim(summary_df)[1]),c(2:dim(summary_df)[2])]
      mat1_filtered<-data.frame(summary_df[[f_level]])
    }else
      {
      mat_otu_final=summary_df[c(1:dim(summary_df)[1]),c(8:dim(summary_df)[2])]

      mat1_filtered<-mat_s
      }

    colnames(mat1_filtered)<-c(f_level)

    repeatsamples <- metadata %>%
      group_by(X.id) %>%
      filter(n() > 1) %>%
      ungroup()
    colnames(mat_otu_final)<-metadata$X.OTU
    #save.image(paste0("data/repeatsamples_",f_level,"_phseq.RData"))
    save(mat_otu_final, repeatsamples,metadata,file = paste0("../data/repeatsamples_",f_level,"_phseq.RData"))
    non_repeatsamples <- metadata %>%
      group_by(X.id) %>%
      filter(n() == 1) %>%
      ungroup()
    colnames(mat_otu_final)<-metadata$X.OTU
    non_repeatsamples<-as.data.frame(non_repeatsamples)
    repeatsamples<-as.data.frame(repeatsamples)
    s_id<-c()
    species_df<-summary_df_species[c(1:dim(summary_df_species)[1]),c(2:dim(summary_df_species)[2])]
    for(i in unique(repeatsamples$X.id))
    {
      read_count<-c()
      ids<-repeatsamples[repeatsamples$X.id==i,]
      for(j in ids$X.OTU)
      {
          s<-sum(species_df[,j])
          read_count<-append(read_count,s)
      }
      s_id<-append(s_id,ids$X.OTU[which(read_count==max(read_count))])
      non_repeatsamples<-bind_rows(non_repeatsamples,ids[which(read_count==max(read_count)),])
    }
    non_repeatsamples$X.OTU
    mat_otu_final<-mat_otu_final[,non_repeatsamples$X.OTU]
    metadata<-non_repeatsamples

    #create phyloseq object
    tt=tax_table(mat1_filtered)
    colnames(tt)<-colnames(mat1_filtered)
    rownames(tt)<-c(mat1_filtered[[f_level]])
    physeq_sam=sample_data(metadata)
    physeq_taxtable=phyloseq::tax_table(tt)
    mat_otu_final<-as.matrix((mat_otu_final))
    colnames(mat_otu_final) <-sample_names(physeq_sam)
    rownames(mat_otu_final)<-rownames(tt)
    physeq_otu=otu_table(mat_otu_final,TRUE)
    physeq = phyloseq(physeq_otu, physeq_taxtable,physeq_sam)
    sample_data(physeq)$X.batch_no<-as.factor(sample_data(physeq)$X.batch_no)
    sample_data(physeq)$X.type<-as.factor(sample_data(physeq)$X.type)
    sample_data(physeq)$Diabetes<-as.factor(sample_data(physeq)$Diabetes)
    sample_data(physeq)$Alcohol<-as.factor(sample_data(physeq)$Alcohol)
    sample_data(physeq)$Smear<-as.factor(sample_data(physeq)$Smear)
    sample_data(physeq)$Gender<-as.factor(sample_data(physeq)$Gender)
    sample_data(physeq)$TOBACCO.Smoking<-as.factor(sample_data(physeq)$TOBACCO.Smoking)
    sample_data(physeq)$Age<-as.numeric(sample_data(physeq)$Age)
    sample_data(physeq)$Height<-as.numeric(sample_data(physeq)$Height)
    sample_data(physeq)$Weight<-as.numeric(sample_data(physeq)$Weight)

    save(physeq,file=paste0("../data/",f_level,"_phseq.RData"))
}

batch_effect_analysis<-function(t_level)
{
  
  load(paste0("../data/",t_level,"_phseq.RData"))
  
  metadata<-sample_data(physeq)
  class(metadata)="data.frame"
  mat_otu_final<-(logratio.transfo(t(otu_table(physeq)+1),logratio = 'CLR'))
  data<-mixOmics::pca(as.matrix(mat_otu_final),ncomp=3)
  
  pca_com<-as.data.frame(do.call(rbind, data$variates))
  
  variance<-do.call(rbind, data$prop_expl_var)
  density.lwd = 0.2
  title.cex = 1.5
  legend.cex = 0.7
  legend.title.cex = 0.75
  color.set <- color.mixo(seq_len(10))
  
  batch <- as.factor(metadata$X.batch_no)
  Type<-as.factor(metadata$X.type)
  pMain <- ggplot(data = pca_com, aes(x = PC1, y = PC2,
                                      color = batch)) +
    geom_point(size=2,alpha=0.8) + xlab(paste0('PC1: ', round(as.numeric(variance[1]*100)),
                                               '% expl.var')) +
    ylab(paste0('PC2: ', round(as.numeric(variance[2]*100)),
                '% expl.var')) +
    scale_shape_manual(values = c(1,19,2,17,4)) +
    scale_color_manual(values = color.set) + 
    theme_bw() +
    
    scale_x_continuous(limits = c(-18,18)) + scale_y_continuous(limits = c(-15,25)) +
    #scale_color_manual(labels=c('TB','TBCOVID','COVID','Control'), values=c("#E69F00", "#00BFC4","#","")) +
    theme(legend.position = 'right', legend.box = 'horizontal',
          legend.direction = 'vertical',
          legend.key.height = unit(0.2, 'cm'),
          legend.key.width = unit(0.1, 'cm'),
          legend.title = element_text(size = rel(legend.title.cex)),
          legend.spacing.x = unit(0.1, 'cm'),
          legend.spacing.y = unit(0.1, 'cm'),
          legend.text = element_text(size = rel(legend.cex)))
  pMain
  
  xlim.update <- layer_scales(pMain)$x$get_limits()
  
  ylim.update <- layer_scales(pMain)$y$get_limits()
  
  pTop <- ggplot(data = pca_com, aes(x = PC1, fill = batch, linetype = NULL)) +
    geom_density(size = density.lwd, alpha = 0.5) + ylab('Density') +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = rel(0.8)),
          plot.title = element_text(hjust = 0.5, size = rel(title.cex)),
          axis.line = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_blank(), legend.position = 'none') +
    scale_fill_manual(values = color.set) +
    scale_x_continuous(limits = xlim.update) 
  
  
  pTop
  pRight <- ggplot(data = pca_com, aes(x = PC2,
                                       fill = batch, linetype = NULL)) +
    geom_density(size = density.lwd, alpha = 0.5) +  coord_flip() +
    ylab('Density') +
    theme(axis.title.x = element_text(size = rel(0.8)),
          axis.title.y = element_blank(), axis.line = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(),
          panel.background = element_blank(), legend.position = 'none') +
    scale_fill_manual(values = color.set) +
    scale_x_continuous(limits = ylim.update)
  
  library(gridExtra)
  legend <- get_legend(pMain)
  PCAplot<-grid.arrange(pTop, legend, pMain + theme(legend.position = 'none'),
                        pRight, ncol = 2, nrow = 2,
                        widths = c(3, 1), heights = c(1, 3))
  
  
  ggsave(paste0("../Figures/manuscript_results/pca_batch_",t_level,".png"),plot = PCAplot, dpi = 300,width=3,height=3)

  
  dist_matrix <- vegdist(t(otu_table(physeq)), method = "bray")
  dissimilarity_values <- as.matrix(dist_matrix)
  calculate_avg_dissimilarity <- function(group, batch, dist_matrix) {
    indices <- which(metadata$X.type == group & metadata$X.batch_no == batch)
    print(indices)
    values <- c()
    if(length(indices)>1){
      for (i in 1:(length(indices) - 1)) {
        for (j in (i + 1):length(indices)) {
          if (i < j) {
            values <- c(values,dissimilarity_values[indices[i],indices[j]])
          }
        }
      }}
    if(length(indices)==1){
      values <- c(values,dissimilarity_values[indices,indices])
    }
    
    return(values)
  }
  
  result_df <- data.frame(Group = character(),
                          Batch = character(),
                          Dissimilarity = numeric(),
                          stringsAsFactors = FALSE)
  for (group in unique(metadata$X.type)) {
    for (batch in unique(metadata$X.batch_no)) {
      
      indices <- which(metadata$X.type == group & metadata$X.batch_no == batch)
      avg_dissimilarity <- calculate_avg_dissimilarity(group, batch, dist_matrix)
      result_df <- rbind(result_df, data.frame(Group = rep(group, length(avg_dissimilarity )),
                                               Batch = rep(batch, length(avg_dissimilarity )),
                                               Dissimilarity = avg_dissimilarity))
    }
  }
  
  print(result_df )
  
  subtype_palette <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854")
  
  p <- ggplot(result_df, aes(x = Group, y = Dissimilarity, fill = Group)) +
    geom_boxplot() +
    #scale_fill_manual(values = subtype_palette) +  # Apply custom colors
    facet_wrap(~ Batch, nrow = 1, scales = "free") +  # Facet in a single row
    theme_bw() +  # White background theme
    theme(legend.position = "bottom")+
    scale_fill_manual(values=c("Control"="#A1C400","COVID"="#B03060","TB"="#E69F00","TBCOVID"= "#00BFC4")) +     
    ylab("BrayCurtis dissimilarity")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels
  p
  ggsave(paste0("../Figures/manuscript_results/bray_",t_level,".png"),plot = p, dpi = 300,width=10,height=5)
  
  
}
