library(patchwork)
source("utils.R")

##########   Fig 3a and S2a ########################################

t_level="genus"
load(paste0("data/",t_level,"_phseq.RData"))
temp<-plot_alphadiversity()
ggsave("Figures/manuscript_results/alphadiversity_genus1.png",plot=temp,dpi=300,width=4,height=6)

t_level="species"
load(paste0("data/",t_level,"_phseq.RData"))
temp<-plot_alphadiversity()
ggsave("Figures/manuscript_results/alphadiversity_species1.png",plot=temp,dpi=300,width=5,height=6)


##########   Fig 3b and S2b ########################################

load("data/species_phseq.RData")
pair<-c("Control","TB")
plot1 <- plot_pcoa()
pair<-c("Control","COVID")
plot2 <- plot_pcoa()
pair<-c("Control","TBCOVID")
plot3 <- plot_pcoa()
pair<-c("TB","TBCOVID")
plot4 <- plot_pcoa()

final_plot <- wrap_plots(list(plot1,plot2,plot3,plot4), ncol = 2)
final_plot
ggsave("Figures/manuscript_results/pca_species.png",plot=final_plot,dpi=300)

load("data/genus_phseq.RData")
pair<-c("Control","TB")
plot1 <- plot_pcoa()
pair<-c("Control","COVID")
plot2 <- plot_pcoa()
pair<-c("Control","TBCOVID")
plot3 <- plot_pcoa()
pair<-c("TB","TBCOVID")
plot4 <- plot_pcoa()

final_plot <- wrap_plots(list(plot1,plot2,plot3,plot4), ncol = 2)
final_plot
ggsave("Figures/manuscript_results/pca_genus.png",plot=final_plot,dpi=300)


##########   Fig 3c and S2c ########################################

t_level="genus"
load(paste0("data/",t_level,"_phseq.RData"))
temp<-plot_toptaxa(t_level)
ggsave("Figures/manuscript_results/top_genus.png",plot=temp,dpi=300,width=5,height=6)

t_level="species"
load(paste0("data/",t_level,"_phseq.RData"))
temp<-plot_toptaxa(t_level)
ggsave("Figures/manuscript_results/top_species.png",plot=temp,dpi=300,width=5,height=6)
