source("DA_analysis_functions.R")

grouptype="type"
bothTB=FALSE
conf<-0.20
f_width<-21       #genus,species 12 pathway 21
m_width=14        #species 14 
level<-c("family" )

for(t_level in level)
{
      
      load(paste0("../data/",t_level,"_phseq.RData"))
      metadata<-sample_data(physeq)
      print(t_level)
      if(t_level=="family")
      {
        print(t_level)
        conf<-0.01
      }
      
      type1="TB"
      type2="TBCOVID"
      TB_TBCOVID<-DA_groups(physeq,type1,type2)
      
      type1="Control"
      type2="TB"
      Control_TB<-DA_groups(physeq,type1,type2)
      
      type1="Control"
      type2="TBCOVID"
      Control_TBCOVID<-DA_groups(physeq,type1,type2)
      
      type1="Control"
      type2="COVID"
      Control_COVID<-DA_groups(physeq,type1,type2)
      
      if(length(Control_TB)!=0){
        Control_TB["Type"]=rep("Control TB",dim(Control_TB)[1])
      }
      if(length(Control_TBCOVID)!=0){
        Control_TBCOVID["Type"]=rep("Control TBCOVID",dim(Control_TBCOVID)[1])
      }
      if(length(Control_COVID)!=0){
        Control_COVID["Type"]=rep("Control COVID",dim(Control_COVID)[1])
      }
      
      df<-rbind(Control_TB,Control_TBCOVID,Control_COVID)
      p<-ggplot()+
        geom_boxplot(data=df %>% filter(!is.na(variable)),aes(x=variable,y=value,fill=`Subject Group`),na.rm = TRUE)+
        facet_nested(. ~ `Type`,drop = TRUE, scales='free', space='free_y', switch='y',
                     strip=strip_nested(text_y=list(element_text(angle=0))),
                     labeller=labeller(group=label_wrap_gen(width=10),
                                       sub_group=label_wrap_gen(width=10)))+
        coord_flip() +
        scale_fill_manual(values=c("Control"="#A1C400","COVID"="#B03060","TB"="#E69F00","TBCOVID"= "#00BFC4")) +
        theme(legend.position="top", legend.key=element_blank(),
              legend.title=element_text(size=16), legend.text=element_text(size=16),
              axis.title.x=element_blank(), axis.text.x=element_text(size=16),
              axis.title.y=element_blank(), axis.text.y=element_text(size=16),
              strip.text=element_text(size=16),
              strip.background=element_rect(fill='gray90', color='gray'),
              strip.placement="outside", panel.spacing.y=unit(0.5, "lines"))
      save(p, file = paste0("../Figures/manuscript_results/",t_level,"_",grouptype,"_saved_plot.RData"))
      ggsave(paste0("../Figures/manuscript_results/",t_level,"_",grouptype,".png"), plot = p, dpi = 300,height=14,width=m_width)
    
    
}

