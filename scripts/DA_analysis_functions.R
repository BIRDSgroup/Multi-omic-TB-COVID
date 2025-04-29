library(ANCOMBC)
library(corncob)
library(LinDA)
library(pheatmap)
library(reshape2)
library(dplyr)
library(forcats)
library(ggplot2)
library(ggh4x)
library(phyloseq)

DA_groups<-function(physeq,type1,type2)
{
  new_directory <- paste0("../Figures/DA_analysis/",type1,"_",type2,"_",t_level,"_",grouptype)
  
  if (!dir.exists(new_directory)) {
    dir.create(new_directory, recursive = TRUE)  
  }
  
  # Set the working directory to the newly created directory
  setwd(new_directory)
  
  physeq_filtered <- subset_samples(physeq, X.type %in% c(type1, type2))
  
  num_samples <- nsamples(physeq_filtered)
  threshold_samples <- ceiling(0.05 * num_samples)
  
  # Step 3: Calculate the number of samples with at least 5 reads for each taxon
  taxa_with_5_reads <- rowSums(otu_table(physeq_filtered) > 5) > threshold_samples
  
  # Step 4: Subset the phyloseq object to include only the taxa that meet the threshold
  physeq_sub <- prune_taxa(taxa_with_5_reads, physeq_filtered)
  
  
  if(type1=="TB" & type2=="TBCOVID")
  {
    bothTB=TRUE
  }
  
  sig_taxa<-sig_taxa(physeq_sub,type1,type2,bothTB)
  setwd("../../../scripts")
  return(sig_taxa)
}




DA_adjusted<-function(physeq,bothTB,t_level,adjusted=TRUE)
{
  metadata<-data.frame(sample_data(physeq))
  if(adjusted)
  {
  if(bothTB)
  {
    c_formula="~X.type+X.batch_no+Diabetes+TOBACCO.Smoking+Alcohol+Smear+Age+Height+Weight+Gender"
    L <- matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0,0,0,0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0,0,0,0), nrow = 2, byrow = TRUE)
  }
  else
  {
    c_formula="~X.type+X.batch_no"
    L <- matrix(c(1, 0, 0,0,0,0, 0,1,0,0,0,0), nrow = 2, byrow = TRUE)
  }
    descp="Adjusted"
    }
  else
  {
    c_formula="~X.type"
    L <- matrix(c(1, 0, 0,1), nrow = 2, byrow = TRUE)
    descp="Undjusted"
  }
 print(c_formula)    
 corncob_out<-differentialTest(formula = as.formula(c_formula),
                               formula_null = ~ 1,
                               phi.formula =  as.formula(c_formula), 
                               phi.formula_null =  as.formula(c_formula),
                               test = "LRT",fdr="BH",data = physeq,boot=F,fdr_cutoff = conf)
 corncob_taxa<-corncob_out$significant_taxa
 corncob_out$p_fdr[corncob_out$significant_taxa]
 
 ancom_out<- ancombc( phyloseq = physeq, 
             formula=substring(c_formula,2),
             tax_level=t_level,
             p_adj_method = "BH", 
             prv_cut=0,
             lib_cut = 0, 
             group = "X.type", 
             struc_zero = TRUE, 
             neg_lb = FALSE, 
             tol = 1e-5, 
             max_iter = 100, 
             conserve = TRUE, 
             alpha = 0.05, 
             global = TRUE
             )
 ancom_taxa<-ancom_out$res$q_val$taxon[ancom_out$res$q_val[[3]]<=conf]
 ancom_dtaxa<-ancom_out$res$q_val$taxon[ancom_out$res$q_val[[3]]==0]
 all_taxa<-ancom_out$res$q_val$taxon
 
 linda_out<- linda(otu_table(physeq),
       metadata,
       formula=c_formula,
       prev.cut = 0,
       p.adj.method = "BH",
       alpha = conf)
 
 lout<-linda.wald.test(linda_out, L, 'LM', alpha = conf)
 
 linda_taxa<-rownames(lout)[lout$reject]

 taxa<-list(c=corncob_taxa,a=ancom_taxa,l=linda_taxa,a_d=ancom_dtaxa)
 
 df<-data.frame("Corncob p-val"=corncob_out$p[all_taxa],"Corncob q-val"=corncob_out$p_fdr[all_taxa],
                "ANCOMBC p-val"=c(ancom_out$res$p_val[ancom_out$res$p_val$taxon==all_taxa,][[3]]),"ANCOMBC q-val"=c(ancom_out$res$q_val[ancom_out$res$q_val$taxon==all_taxa,][[3]]),"ANCOMBC LFC"=c(ancom_out$res$lfc[ancom_out$res$lfc$taxon==all_taxa,][[3]]),
                "LinDA p-val"=lout[all_taxa,"pvalue"],"LinDA q-val"=lout[all_taxa,"padj"],"LinDA LFC"=linda_out$output[[1]][all_taxa,"log2FoldChange"])
 
 write.csv(df, file = paste0(type1,"_",type2,"_",t_level,"_",descp,".csv"), row.names = TRUE)
 
 
 
 return(list(corncob_out = corncob_out,ancom_out=ancom_out,linda_out=linda_out,taxa=taxa))
 
}

sig_taxa<-function(physeq,type1,type2,bothTB)
{
  output1<-DA_adjusted(physeq,bothTB,t_level)
  
  
  output2<-DA_adjusted(physeq,bothTB,t_level,FALSE)
  
  taxa<-sort(unique(c(unlist(output1$taxa),unlist(output2$taxa))))
  
  if(length(taxa)!=0)
  {
    df<-data.frame("Taxa"=taxa,"Unadjusted"=rep("",length(taxa)),"Adjusted"=rep("",length(taxa)))
    
    method<-c("Unadjusted","Adjusted")
    for(m in method)
    {
      if(m=="Unadjusted")
      {
        t<-"output2"}
      else
      {
        t<-"output1"}
      output<-get(t)
      for(i in taxa)
      {
        if(i %in% output$taxa$c)
        {
          df[df["Taxa"]==i,m]<-paste0(df[df["Taxa"]==i,m],"$")
        }
        if(i %in% output$taxa$a)
        {
          df[df["Taxa"]==i,m]<-paste0(df[df["Taxa"]==i,m],"+")
        }
        if(i %in% output$taxa$l)
        {
          df[df["Taxa"]==i,m]<-paste0(df[df["Taxa"]==i,m],"*")
        }
        if(i %in% output$taxa$a_d)
        {
          df[df["Taxa"]==i,m]<-paste0(df[df["Taxa"]==i,m],"^") #exclusive taxa present in a group
        }
          
      }
    }
    
    write.csv(df, file = paste0(type1,"_",type2,"_",t_level,"_","significant.csv"), row.names = FALSE)
    plot_taxa<-df[nchar(df$Adjusted) >= 2, ]$Taxa
    
    #Plotting
    samp_frac <- output1$ancom_out$samp_frac
    samp_frac[is.na(samp_frac)] <- 0
    ps <- log(otu_table(physeq) + 1) -samp_frac

    lfc<-output1$ancom_out$res$lfc[output1$ancom_out$res$lfc$"taxon" %in% plot_taxa,][[3]]
    upper<-c(lfc + 1.96 * output1$ancom_out$res$se[output1$ancom_out$res$se$"taxon" %in% plot_taxa,][[3]])
    lower<- c(lfc - 1.96 * output1$ancom_out$res$se[output1$ancom_out$res$se$"taxon" %in% plot_taxa,][[3]])
    
    ancom_df<-data.frame("LFC"=lfc,"CI.upper"=upper,"CI.lower"=lower)
    ancom_df["Taxa"]<-output1$ancom_out$res$lfc[output1$ancom_out$res$lfc$"taxon" %in% plot_taxa,]$taxon
    ancom_df["Method"]<-rep("ANCOM-BC",length(plot_taxa))
    
    
    lfc<-output1$linda_out$output[[1]][plot_taxa,"log2FoldChange"]
    upper<-c(lfc + 1.96 * output1$linda_out$output[[1]][plot_taxa,"lfcSE"])
    lower<- c(lfc - 1.96 * output1$linda_out$output[[1]][plot_taxa,"lfcSE"])
    linda_df<-data.frame("LFC"=lfc,"CI.upper"=upper,"CI.lower"=lower)
    linda_df["Taxa"]<-rownames(output1$linda_out$output[[1]][plot_taxa,])
    linda_df["Method"]<-rep("LinDA",length(plot_taxa))
    LFC_df<-rbind(ancom_df,linda_df)
    
    #####only for manuscript plotting get the top pathways <=1.5 linda LFC
    abs_lfc<-abs(linda_df$LFC)
    sorted_linda<- linda_df[order(-abs(linda_df$LFC)),]
    top_taxa <- sorted_linda[abs(sorted_linda$LFC)>=1.5,"Taxa"]
    #l_taxa<-load(top_taxa.RData)
    #top_taxa <- sorted_linda[,l_taxa]
    #######################################################################
    
    ps<-t(ps)
    if(length(plot_taxa)!=0){
      ps<-ps[,plot_taxa]
      ps<-as.data.frame(ps)
      ps["Subject Group"]<-sample_data(physeq)$X.type
      
      melted_df<-melt(ps,id.vars="Subject Group")
      
     
      LFC_df["Subject Group"]<-rep(NA,dim(LFC_df)[1])
      LFC_df["variable"]<-rep(NA,dim(LFC_df)[1])
      LFC_df["value"]<-rep(NA,dim(LFC_df)[1])
      LFC_df["Plot type"]<-rep("Log Fold Change",dim(LFC_df)[1])
      
      melted_df["LFC"]<-rep(NA,dim(melted_df)[1])
      melted_df["CI.upper"]<-rep(NA,dim(melted_df)[1])
      melted_df["CI.lower"]<-rep(NA,dim(melted_df)[1])
      melted_df["Taxa"]<-rep(NA,dim(melted_df)[1])
      melted_df["Method"]<-rep(NA,dim(melted_df)[1])
      melted_df["Plot type"]<-rep("Log (Bias corrected Abundances)",dim(melted_df)[1])
      
      #######################Manuscript Plotting###################
      melted_df<-melted_df[melted_df$variable %in% top_taxa,]
      LFC_df<-LFC_df[LFC_df$Taxa %in% top_taxa,]               
      #############################################################
      
      plot_df<-rbind(melted_df,LFC_df)
      LFC_df<-LFC_df[order(LFC_df$LFC),]
      
      plot_df
            
      odd<-order(plot_df["Taxa"])
      
      
      p<-ggplot()+
        geom_boxplot(data=plot_df %>% filter(!is.na(variable)),aes(x=variable,y=value,fill=`Subject Group`),na.rm = TRUE)+
        geom_errorbar(data=plot_df %>% filter(!is.na(Taxa)),aes(x=Taxa, ymin=CI.lower, ymax=CI.upper,color=Method),position = position_dodge(0.75),width=0)+
        geom_point(data=plot_df %>% filter(!is.na(Taxa)),aes(x=Taxa, y=LFC, color=Method), position=position_dodge(0.75), size=1.75)+
        geom_hline(data=plot_df %>% filter(!is.na(Taxa)),aes(yintercept=0),size=0.5, linetype='dashed', alpha=0.5)+
        facet_nested(. ~ `Plot type`,drop = TRUE, scales='free', space='free_y', switch='y',
                     strip=strip_nested(text_y=list(element_text(angle=0))),
                     labeller=labeller(group=label_wrap_gen(width=10),
                     sub_group=label_wrap_gen(width=10))) +
        scale_fill_manual(values=c("Control"="#A1C400","COVID"="#B03060","TB"="#E69F00","TBCOVID"= "#00BFC4","cured/lost_to_followup"="#E69F00","relapse/DR/died"= "#00BFC4")) +
        coord_flip() +
        scale_linetype_manual(values=c("11", "solid"))+
        scale_color_manual(values=c("blue", "red"))+
        scale_linetype_manual(values=c("11", "solid")) +
        
          theme(legend.position="top", legend.key=element_blank(),
                legend.title=element_text(size=16), legend.text=element_text(size=16),
                axis.title.x=element_blank(), axis.text.x=element_text(size=16),
                axis.title.y=element_blank(), axis.text.y=element_text(size=16),
                strip.text=element_text(size=16),
                strip.background=element_rect(fill='gray90', color='gray'),
                strip.placement="outside", panel.spacing.y=unit(0.5, "lines"))
      save(p, file = paste0(type1,"_",type2,"_",t_level,"_saved_plot.RData"))
      save(top_taxa,file=("top_taxa.RData"))
      ggsave(paste0(type1,"_",type2,"_",t_level,"_.png"), plot = p, ,height=8,width=f_width, units = "in", dpi = 300)
      return(melted_df[,c("Subject Group","variable","value")])
    }
  }
  return(data.frame())
}

