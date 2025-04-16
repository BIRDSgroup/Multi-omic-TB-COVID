library(phyloseq)
library(NetCoMi)

load(paste0("data/species_phseq.RData"))

physeq_filtered <- subset_samples(physeq, X.type %in% c("TB"))
#otu_table(physeq_filtered)[otu_table(physeq_filtered)<5]=0
taxa_to_keep <- rowSums(otu_table(physeq_filtered) > 0) / ncol(otu_table(physeq_filtered))>=0.1
otu_data_filtered <- otu_table(physeq_filtered)[taxa_to_keep, ]

physeq1 <- phyloseq(otu_table(otu_data_filtered, taxa_are_rows = taxa_are_rows(physeq_filtered)),
                       sample_data(physeq_filtered),
                       tax_table(physeq_filtered))

physeq_filtered <- subset_samples(physeq, X.type %in% c("TBCOVID"))
#otu_table(physeq_filtered)[otu_table(physeq_filtered)<5]=0
taxa_to_keep <- rowSums(otu_table(physeq_filtered) > 0) / ncol(otu_table(physeq_filtered))>=0.1
otu_data_filtered <- otu_table(physeq_filtered)[taxa_to_keep, ]

physeq2 <- phyloseq(otu_table(otu_data_filtered, taxa_are_rows = taxa_are_rows(physeq_filtered)),
                    sample_data(physeq_filtered),
                    tax_table(physeq_filtered))


net_genus <- netConstruct(data=physeq1,data2=physeq2, 
                          taxRank = "species",
                          measure = "spieceasi",
                          measurePar = list(method = "mb",pulsar.params = list(rep.num = 10),symBetaMode = "ave"),
                          sparsMethod = "threshold",
                          thresh=0.3,
                          normMethod = "none",
                          dissFunc = "signed",
                          verbose = 3,seed = 123456)


diff_net <- diffnet(net_genus,adjust = "adaptBH",alpha=0.05)


plot(diff_net, 
     cexNodes = 3,
     cexlabels=4,
     cexTitle = 4,
     mar = c(2,2,8,5),
     legendGroupnames = c("group 'no'", "group 'yes'"),
     legendPos = c(0.7,1.6))

net_genus_analyse <- netAnalyze(net_genus, 
                           centrLCC = FALSE,
                           avDissIgnoreInf = TRUE,
                           sPathNorm = FALSE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = c("degree", "eigenvector"),
                           hubQuant = 0.9,
                           lnormFit = TRUE,
                           normDeg = FALSE,
                           normBetw = FALSE,
                           normClose = FALSE,
                           normEigen = FALSE)


plot(net_genus_analyse, 
     sameLayout = TRUE, 
     nodeColor = "cluster",
     nodeSize = "fix",
     labelScale = FALSE,
     rmSingles="inboth",
     cexNodes = 1, 
     cexLabels = 0.3,
     cexHubLabels = 0.2,
     cexTitle = 2,
     groupNames = c("TB", "TBCOVID"),
     hubBorderCol  = "gray40")

png_file<-paste0("Figures/manuscript_results/network.png")
png(png_file, res = 300, width = 2500, height = 1600)
stdout <- capture.output(print(plot(net_genus_analyse,
                                    sameLayout = TRUE, 
                                    nodeColor = "cluster",
                                    nodeSize = "fix",
                                    labelScale = FALSE,
                                    rmSingles="inboth",
                                    cexNodes = 1, 
                                    cexLabels = 0.4,
                                    cexHubLabels = 0.5,
                                    cexTitle = 1,
                                    groupNames = c("TB", "TBCOVID"),
                                    hubBorderCol  = "gray40")))
dev.off()

#legend("bottom", title = "estimated association:", legend = c("+","-"), 
#       col = c("#009900","red"), inset = 0.02, cex = 4, lty = 1, lwd = 4, 
#       bty = "n", horiz = TRUE)
