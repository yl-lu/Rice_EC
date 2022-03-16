# clustering of EC in rice
library(ClustVarLV)
library(ClusterR)
library(pheatmap)
library(grid)
library(amap)
library(RColorBrewer)

setwd("C:/project-rice_EC")

oridf <- read.table("TPMmeans.csv",
                 header = T,
                 row.names = 1,
                 sep = ",")

df <- oridf[,c("WT_ZT06_SD","WT_ZT06_LD",
               "WT_ZT14_SD","WT_ZT14_LD",
               "WT_ZT18_SD","WT_ZT18_LD",
               "Elf3_ZT06_SD","Elf3_ZT06_LD",
               "Elf3_ZT14_SD","Elf3_ZT14_LD",
               "Elf3_ZT18_SD","Elf3_ZT18_LD",
               "Phyb_ZT06_SD","Phyb_ZT06_LD",
               "Phyb_ZT14_SD","Phyb_ZT14_LD",
               "Phyb_ZT18_SD","Phyb_ZT18_LD")]

df <- log2(df+1)
df <- df[apply(df,1,max)>2, ]
WT_SDminusLD <- data.frame(WT_ZT06_SDminusLD = df$WT_ZT06_SD - df$WT_ZT06_LD,
                           WT_ZT14_SDminusLD = df$WT_ZT14_SD - df$WT_ZT14_LD,
                           WT_ZT18_SDminusLD = df$WT_ZT18_SD - df$WT_ZT18_LD)
rownames(WT_SDminusLD) <- rownames(df)
data <- scale(WT_SDminusLD)
opt_gmm <- Optimal_Clusters_GMM(data,
                               max_clusters = 20,
                               criterion = "BIC", 
                               dist_mode = "eucl_dist",
                               seed_mode = "static_spread",
                               km_iter = 30,
                               em_iter = 30,
                               var_floor = 1e-10, 
                               plot_data = T)
gmm <- GMM(data,
           gaussian_comps = 8,
           dist_mode = "eucl_dist",
           seed_mode = "static_spread",
           km_iter = 30,
           em_iter = 30,
           verbose = F)
pr <- predict_GMM(data,
                  gmm$centroids,
                  gmm$covariance_matrices,
                  gmm$weights)
names(pr$cluster_labels) <- rownames(data)
clustersTree <- sort(pr$cluster_labels)
gaps <- which((clustersTree[-1] - clustersTree[-length(clustersTree)]) != 0)
gene.cluster <- as.data.frame(clustersTree)
df.cluster <- cbind(gene.cluster,data[rownames(gene.cluster),])
df.cluster[,-1][df.cluster[,-1] > 2] <- 2
df.cluster[,-1][df.cluster[,-1] < -2] <- -2
my_gene_clusters <- data.frame(cluster = df.cluster[,1])
rownames(my_gene_clusters) <- rownames(df.cluster)
colorsVec <- brewer.pal(8, "Set2")
names(colorsVec) <- c(1:8)
my_colors <- list(cluster = colorsVec)
phm <- pheatmap(df.cluster[,-1],
                color = colorRampPalette(c("#c1207e","snow","#669900"))(200),
                border_color = "snow",
                cluster_rows = F,
                cluster_cols = F,                
                treeheight_row = 30, 
                treeheight_col = 30,
                fontsize_col = 8,
                fontsize_row = 10,
                cutree_rows = 12,
                angle_col = 90,
                gaps_row = gaps,
                annotation_colors = my_colors,
                annotation_row = my_gene_clusters,
                annotation_legend = T,
                labels_col = c("ZT6","ZT14","ZT18"),
                annotation_names_row = T,                  
                annotation_names_col = T,
                show_rownames = F,
                display_numbers = F,
                drop_levels = T
)
# ----
# line plot
#-----------------------------
library(ggplot2)
library(fields)
setwd("C:/project-rice_EC")
df <- read.table("clustering.originalTPM_correspondingWithRatio.8cluster_labels.xls", header = T, row.names = 1)
mycolnames <- colnames(df)

df <- data.frame(log2(df[,1:18]+1),
                 df$cluster)
colnames(df) <- mycolnames
head(df)

df.select <- read.table("select.geneID.clusters.geneSymbol.TPM.xls",
                        header = T, row.names = 1, sep = "\t")
df.lux_ZT14_SD <- read.table("TPM_lux_ZT14_SD.txt",
                             header = T, row.names = 1, sep = "\t")
df.select$lux_ZT14_SD <- df.lux_ZT14_SD[df.select$gene_id,]

df.select <- df.select[,-1]

df.select <- data.frame(log2(df.select[,c(1:18,20)]+1),
                        df.select$cluster)
colnames(df.select) <- gsub("df.select.cluster","cluster",colnames(df.select))
head(df.select)

cluster <- NULL
for (k in 1:8){
  cluster[[k]] <- df[df$cluster==(k-1),]
}

gene <- NULL
p <- NULL

for (i in 1:length(df.select[,1])){
  gene[[i]] <- df.select[i,]
  reshapeList <- strsplit(colnames(gene[[i]]), split = '_')
  myline <- NULL
  for (j in 1:(length(reshapeList)-1)){
    myline[[j]] <- as.vector(unlist(c(reshapeList[[j]], gene[[i]][j])))
  }
  df.gene <- data.frame(matrix(unlist(myline),
                               nrow = length(reshapeList)-1,
                               byrow = TRUE),
                        stringsAsFactors = FALSE)
  colnames(df.gene) <- c("Genotype","ZT","Day_length","TPM")
  df.gene$Genotype <- gsub("elf3", "elf3-1 elf3-2", df.gene$Genotype)
  df.gene <- subset(df.gene,df.gene$Genotype == "WT" |
                            df.gene$Genotype == "elf3-1 elf3-2" |
                            df.gene$Genotype == "lux")
  df.gene$TPM <- as.numeric(df.gene$TPM)
  df.gene$cluster <- rep(gene[[i]]$cluster,length(df.gene$TPM))
  df.gene$ZT <- as.integer(gsub("ZT","",df.gene$ZT))
  df.gene$Genotype <- factor(df.gene$Genotype, levels = c("WT","phyb","elf3-1 elf3-2","lux"))
  df.gene$Day_length <- as.factor(df.gene$Day_length)
  df.gene$symbol <- rep(rownames(df.select[i,]),length(df.gene$TPM))
  
  mysymbol <- unique(df.gene$symbol)
  myclusterName <- unique(df.gene$cluster)
  
  mycluster <- cluster[[unique(df.gene$cluster)+1]]
  mysummary <- summary(mycluster)
  
  df.gene$sd <- apply(mycluster,2,sd)[1:length(rownames(df.gene))]
  df.gene$mean <- apply(mycluster,2,mean)[1:length(rownames(df.gene))]
  df.gene$median<- apply(mycluster,2,median)[1:length(rownames(df.gene))]
  df.gene$upperquartile <- as.numeric(gsub("1st Qu.:","",mysummary[2,][1:length(rownames(df.gene))]))
  df.gene$lowerquartile <- as.numeric(gsub("3rd Qu.:","",mysummary[5,][1:length(rownames(df.gene))]))

  p[[i]] <- ggplot(df.gene) +
    geom_point(aes(x = ZT, y = TPM,
                   color = Day_length,
                   group = Day_length),
               size = 1.5) +
    geom_ribbon(aes(x = ZT, y = TPM,
                    ymin = median-sd, ymax = median+sd,
                    alpha = 0.02,
                    group = Day_length,
                    fill = Day_length),
                inherit.aes = F) +
    geom_line(aes(x = ZT, y = TPM,
                  color = Day_length,
                  group = Day_length),
              size = 1.5) +
    scale_color_manual(values = c("red2","deepskyblue3")) +
    scale_fill_manual(values = c("lightcoral","lightskyblue")) +
    facet_wrap(~Genotype, scales = "fixed") +
    scale_x_continuous(breaks = c(6,14,18),
                       labels = c("6","14","18"),
                       position = 'bottom') +
    theme(axis.text.x = element_text(size = 21,
                                     color = "black",
                                     face = "plain",
                                     vjust = 0,
                                     hjust = 0.5,
                                     angle = 0)) +
    theme(axis.text.y = element_text(size = 21,
                                     color = "black",
                                     face = "plain",
                                     vjust = 0.5,
                                     hjust = 1,
                                     angle = 0)) +
    theme(axis.title.x = element_text(size = 21,
                                      color = "black",
                                      face = "plain",
                                      vjust = -0.3,
                                      hjust = 0.5,
                                      angle = 0)) +
    theme(legend.title = element_text(size = 21,
                                      color = "black",
                                      face = "plain",
                                      vjust = 5,
                                      hjust = 0.5,
                                      angle = 0)) +
    theme(legend.text = element_text(size = 18,
                                     color = "black",
                                     face = "plain",
                                     angle = 0)) +
    theme_bw() +
    theme(panel.border = element_blank(),
          panel.grid = element_blank(),
          strip.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black")) +
    labs(x = paste(mysymbol," Cluster ",myclusterName, sep = "")) +
    labs(y = "log2(TPM+1)") 
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
pdf("lineplots.WT_elf3-1elf3-2_lux.selectedGenes.pdf", width = 5, height = length(rownames(df.select))*(20/11))
multiplot(p[[1]],p[[2]],p[[3]],
          p[[4]],p[[5]],p[[6]],
          p[[7]],p[[8]],p[[9]],
          p[[10]],p[[11]],p[[12]],
          p[[13]],p[[14]],p[[15]],
          p[[16]],p[[17]],p[[18]],
          p[[19]],p[[20]],p[[21]],
          p[[22]],p[[23]],p[[24]],
          p[[25]],p[[26]],p[[27]],
          p[[28]],
          cols = 1)
dev.off()

#---
# extract gene sets after clustering WT_SD/WT_LD,
# such as extract log2(TPM ratio) for each cluster in phyb_LD/WT_LD
setwd("C:/project-rice_EC")
oridf <- read.table("TPMmeans.csv",
                    header = T,
                    row.names = 1,
                    sep = ",")
df <- oridf[,c("WT_ZT06_SD","WT_ZT06_LD",
               "WT_ZT14_SD","WT_ZT14_LD",
               "WT_ZT18_SD","WT_ZT18_LD",
               "Elf3_ZT06_SD","Elf3_ZT06_LD",
               "Elf3_ZT14_SD","Elf3_ZT14_LD",
               "Elf3_ZT18_SD","Elf3_ZT18_LD",
               "Phyb_ZT06_SD","Phyb_ZT06_LD",
               "Phyb_ZT14_SD","Phyb_ZT14_LD",
               "Phyb_ZT18_SD","Phyb_ZT18_LD")]
df <- log2(df+1)
df.cluster <- read.table("WT_log_SD-LD.plotOrder.xls", header = T,
                         row.names = 1, sep = "\t")
head(df.cluster)
head(df)
# WT_SDminusLD <- data.frame(WT_ZT06_SDminusLD = df$WT_ZT06_SD-df$WT_ZT06_LD,
#                            WT_ZT14_SDminusLD = df$WT_ZT14_SD-df$WT_ZT14_LD,
#                            WT_ZT18_SDminusLD = df$WT_ZT18_SD-df$WT_ZT18_LD)
# rownames(WT_SDminusLD) <- rownames(df)

phyb_LDminusWT_LD <- data.frame(ZT06_phyb_LDminusWT_LD = df$Phyb_ZT06_LD - df$WT_ZT06_LD,
                                ZT14_phyb_LDminusWT_LD = df$Phyb_ZT14_LD - df$WT_ZT14_LD,
                                ZT18_phyb_LDminusWT_LD = df$Phyb_ZT18_LD - df$WT_ZT18_LD)
rownames(phyb_LDminusWT_LD) <- rownames(df)
head(phyb_LDminusWT_LD)
logTPMRatio_phybLDminusWTLD <- cbind(df.cluster$clustersTree, phyb_LDminusWT_LD[rownames(df.cluster),])
colnames(logTPMRatio_phybLDminusWTLD) <- c("Cluster","ZT06_phyb_LD_vs_WT_LD","ZT14_phyb_LD_vs_WT_LD","ZT18_phyb_LD_vs_WT_LD")
head(logTPMRatio_phybLDminusWTLD)
write.table(logTPMRatio_phybLDminusWTLD,
            file = "phyb_LD_vs_WT_LD_clusteringOrder.xls",
            row.names = T, col.names = T, quote = FALSE, sep = '\t')

# mark out the ELF3-bound genes in clustering heatmap of all ZTs in WT_SD/WT_LD
allWTdf <- oridf[rownames(df.cluster),]
allWTdf <- allWTdf[,c("WT_ZT00_SD","WT_ZT00_LD",
                      "WT_ZT02_SD","WT_ZT02_LD",
                      "WT_ZT06_SD","WT_ZT06_LD",
                      "WT_ZT10_SD","WT_ZT10_LD",
                      "WT_ZT14_SD","WT_ZT14_LD",
                      "WT_ZT18_SD","WT_ZT18_LD")]
allWTdf <- log2(allWTdf+1)
str(allWTdf)
allWT_SDminusLD <- data.frame(WT_ZT00_SDminusLD = allWTdf$WT_ZT00_SD - allWTdf$WT_ZT00_LD,
                              WT_ZT02_SDminusLD = allWTdf$WT_ZT02_SD - allWTdf$WT_ZT02_LD,
                              WT_ZT06_SDminusLD = allWTdf$WT_ZT06_SD - allWTdf$WT_ZT06_LD,
                              WT_ZT10_SDminusLD = allWTdf$WT_ZT10_SD - allWTdf$WT_ZT10_LD,
                              WT_ZT14_SDminusLD = allWTdf$WT_ZT14_SD - allWTdf$WT_ZT14_LD,
                              WT_ZT18_SDminusLD = allWTdf$WT_ZT18_SD - allWTdf$WT_ZT18_LD)
rownames(allWT_SDminusLD) <- rownames(allWTdf)
my_gene_clusters <- data.frame(cluster = df.cluster[,1])
rownames(my_gene_clusters) <- rownames(df.cluster)
library(RColorBrewer)
colorsVec <- brewer.pal(8, "Set2")
names(colorsVec) <- c(1:8)
my_colors <- list(cluster = colorsVec)
allWT_SDminusLD[allWT_SDminusLD > 2] <- 2
allWT_SDminusLD[allWT_SDminusLD < -2] <- -2

df.reorder <- read.table("ELF3-bound_genes_clusters_sort.xls", row.names = 1)
colnames(df.reorder) <- c("Cluster","if_bound")
head(df.reorder)
head(allWT_SDminusLD)
phm <- pheatmap(allWT_SDminusLD[rownames(df.reorder),],
                color = colorRampPalette(c("#c1207e","snow","#669900"))(200),
                border_color = "snow",
                cluster_rows = F,
                cluster_cols = F,                
                clustering_distance_rows = "correlation",
                clustering_method = "ward.D",                      
                treeheight_row = 30, 
                treeheight_col = 30,
                fontsize_col = 8,
                fontsize_row = 10,
                cutree_rows = 12,
                angle_col = 90,
                gaps_row = gaps,
                # gaps_col = my_gaps_col,
                annotation_colors = my_colors,
                annotation_row = my_gene_clusters,
                annotation_legend = T,
                labels_col = c("ZT0","ZT2","ZT6",
                               "ZT10","ZT14","ZT18"),
                annotation_names_row = T,                  
                annotation_names_col = T,
                show_rownames = F,
                display_numbers = F,
                drop_levels = T
)


elf3BoundGenes <- read.table("C:/project-rice_EC/ChIP-seq/ELF3.targetGenes.bed", header = F, row.names = 4, sep = "\t")
df.cluster <- read.table("WT_log_SD-LD.plotOrder.xls", header = T,
                         row.names = 1, sep = "\t")
df.cluster$clustersTree <- df.cluster$clustersTree + 1 
head(elf3BoundGenes)
head(df.cluster)
my_gene_clusters <- data.frame(cluster = df.cluster[,1])
rownames(my_gene_clusters) <- rownames(df.cluster)
colorsVec <- brewer.pal(8, "Set2")
names(colorsVec) <- c(1:8)
my_colors <- list(cluster = colorsVec)

elf3BoundGenesClusters <- read.table("ELF3-bound_genes_clusters_sort.xls", header = F, row.names = 1)
colnames(elf3BoundGenesClusters) <- c("Cluster","if_bound")
elf3BoundGenesClusters$Cluster <- 0
head(elf3BoundGenesClusters[,-1])
phm <- pheatmap(elf3BoundGenesClusters,
                color = colorRampPalette(c("white","black"))(2),
                border_color = "snow",
                cluster_rows = F,
                cluster_cols = F,                
                treeheight_row = 30, 
                treeheight_col = 30,
                fontsize_col = 8,
                fontsize_row = 10,
                cutree_rows = 12,
                angle_col = 90,
                gaps_row = gaps,
                annotation_colors = my_colors,
                annotation_row = my_gene_clusters,
                annotation_legend = T,
                annotation_names_row = T,                  
                annotation_names_col = T,
                show_rownames = F,
                display_numbers = F,
                drop_levels = T
)
