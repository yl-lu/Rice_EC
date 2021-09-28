library(motifStack)
library(ade4)
library(RColorBrewer)
library(magrittr)

wd <- "C:/project-rice_elf3/motif/ELF3_summits.known/"
#setwd("C:/project-rice_elf3/motif/ELF3_summits.denovo/")
setwd(wd)
pfms <- importMatrix(dir(wd, "pfm$", full.names = TRUE))
hc <- clusterMotifs(pfms)
phylog <- ade4::hclust2phylog(hc)
leaves <- names(phylog$leaves)
pfms <- pfms[leaves]
groupDistance <- 20
motifSig <- motifSignature(pfms,
                           phylog,
                           cutoffPval = 0.0001,
                           groupDistance = groupDistance,
                           rcpostfix = "(RC)",
                           min.freq = 1,
                           trim = 0.3,
                           #families = 12,
                           sort = T)

sig <- signatures(motifSig)
gpCol <- sigColor(motifSig)
color <- brewer.pal(9, "Set1")

pdf(file = paste("groupDistance_",
                 groupDistance,
                 ".ELF3.summits.known.pdf",
                 sep = ""),
    width = 9,
    height = 7)
motifPiles(phylog = phylog,
           pfms = pfms,
           pfms2 = sig,
           r.tree = 0.45,
           col.tree = c(rep(color[1],9), # myb
                        rep(color[2],6), # G2like
                        rep(color[3],2), # C2C2gata
                        rep(color[4],4), # C2C2dof
                        rep(color[5],1), # bZIP
                        rep(color[6],1), # LBD
                        rep(color[7],1), # VRN1
                        rep(color[3],1), # C2C2gata
                        rep(color[5],14), # bZIP
                        rep(color[8],1), # bHLH
                        rep(color[5],1), # bZIP
                        rep(color[8],3), # bHLH
                        rep(color[9],4)), # NAC
           col.leaves = c(rep(color[1],9), # myb
                          rep(color[2],6), # G2like
                          rep(color[3],2), # C2C2gata
                          rep(color[4],4), # C2C2dof
                          rep(color[5],1), # bZIP
                          rep(color[6],1), # LBD
                          rep(color[7],1), # VRN1
                          rep(color[3],1), # C2C2gata
                          rep(color[5],14), # bZIP
                          rep(color[8],1), # bHLH
                          rep(color[5],1), # bZIP
                          rep(color[8],3), # bHLH
                          rep(color[9],4)), # NAC
           col.pfms2 = gpCol,
           #r.anno = c(0.02, 0.03, 0.04),
           # col.anno = list(colors()[1:length(pfms)],
           #               colors()[(length(pfms)+1):(2*length(pfms))],
           #               colors()[(2*length(pfms)+1):(3*length(pfms))]),
           motifScale = "logarithmic",
           plotIndex = TRUE,
           groupDistance = groupDistance,
           groupDistanceLineCol = "grey")
dev.off()
