
source("utils.R")

##########ICC calculation########################
library("psych")
load('../data/repeatsamples_genus_phseq.RData')
counts<-mat_otu_final

shannon_matrix_na <- matrix(NA, nrow = length(unique(repeatsamples$X.id)), ncol = 3)
rownames(shannon_matrix_na)<-unique(repeatsamples$X.id)
for(i in unique(repeatsamples$X.id))
{
  rep_ids<-repeatsamples[repeatsamples$X.id==i,]$X.OTU
  mat_df<-counts[,rep_ids]
  mat_df<-mat_df[!apply(mat_df, 1, function(x) all(x == 0)), ]
  d<-diversity(t(mat_df))
  for(j in 1:length(d))
  {
    shannon_matrix_na[i,j]<-c(d[j])
  }
}
ICC(shannon_matrix_na,missing=TRUE,lmer=TRUE)


##########   Fig S1 ########################################

##########CV calculation########################
library("psych")
ids<-c("PT331","PT258") #the ones with atleast 3 replicates

type<-c("TB","TBCOVID")
load("../data/genus_phseq.RData")


std<-c()
mean<-c()
cv<-c()
type_vector<-c()
rep_type<-c()

mat_otu_final<-t(counts)
total_counts_per_sample <- rowSums(mat_otu_final)
mat_otu_final<-sweep(mat_otu_final, 1, total_counts_per_sample, "/")
mat_otu_final<-mat_otu_final*100
k=1
for(i in ids)
{
  OTU_ids<-repeatsamples[repeatsamples$X.id==i,]$X.OTU
  subset<-t(mat_otu_final[OTU_ids,])
  subset<-subset[!apply(subset, 1, function(x) all(x == 0)), ]
  std<-c(std,apply(subset, 1, sd))
  mean<-c(mean,apply(subset, 1, mean))
  type_vector<-c(type_vector,rep(type[k],dim(subset)[1])  )
  rep_type<-c(rep_type,rep("technical replicates",dim(subset)[1])  )
  print(length(std))
  print(dim(subset))
              
  subset<-sample_data(physeq)[sample_data(physeq)$X.type==type[k],]
  subset<-t(mat_otu_final[subset[!subset$X.id %in% repeatsamples$X.id,]$X.OTU,])
  subset<-subset[!apply(subset, 1, function(x) all(x == 0)), ]
  std<-c(std,apply(subset, 1, sd))
  mean<-c(mean,apply(subset, 1, mean))
  type_vector<-c(type_vector,rep(type[k],dim(subset)[1]))  
  rep_type<-c(rep_type,rep("biological samples",dim(subset)[1]))  

  k=k+1
}

cv<-log10(std/mean)
mean<-log10(mean)

df<-data.frame("Group"=type_vector,"Type"=rep_type,"CV"=cv,"Mean"=mean)

p<-ggplot(df,aes(x=Mean,y=CV,color=Type))+
  ggplot2::geom_point(
    ggplot2::aes(fill = Type),
    alpha = 1 ,
    size = 2,
    shape = 21,
    stroke = 0.35,
    position = ggplot2::position_jitterdodge()
  )+
  scale_fill_manual(values=c("red","blue"))+
  scale_color_manual(values=c("red","blue"))+
  facet_wrap(~Group)+
  labs( x = expression("log"[10]*"mean relative abundance"),
        y =  expression("log"[10]*"CV"))+
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom") +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "bottom")
png_file<-paste0("manuscript_results/CV_group_new.png")
png(png_file, res = 300, width = 1900, height = 1500)
stdout <- capture.output(print(p))
dev.off()

