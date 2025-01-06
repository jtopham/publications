# ------------------------------------------------------------------------------
#                           _      _      _
#                        __(.)< __(.)> __(.)=
#                        \___)  \___)  \___)
# ------------------------------------------------------------------------------
# FILENAME : NET_prj.R
#
# AUTHOR : James Topham
#
# DATE : Jan 2024
#
# DESCRIPTION :
#
#   -  NET analysis
# ------------------------------------------------------------------------------
#                           _      _      _
#                        __(.)< __(.)> __(.)=
#                        \___)  \___)  \___)
# ------------------------------------------------------------------------------
library(ggplot2)
library(reshape2)
library(pheatmap)
library(DEqMS)
library(ggbio)
library(GenomicRanges)

# color palettes
source("/projects/pangen/analysis/jtopham/scripts/custom_pal.R")

# gene sets and GSEA script
load("/projects/pangen/analysis/jtopham/data/msigdb_apr2021.RData")
source("/projects/pangen/analysis/jtopham/scripts/my_gsea.R")

# ------------------------------------------------------------------------------
# load processed data
# ------------------------------------------------------------------------------
load("/projects/pangen/analysis/jtopham/data/all_clin_feb2024.RData")
load("/projects/pangen/analysis/jtopham/data/subtypes_hg38_feb2024.RData")
load("/projects/pangen/analysis/jtopham/data/rna_hg38_feb2024.RData")

# POG/PanGen RNA data
rownames(all.rna) <- all.rna$geneid
all.rna <- all.rna[, -c(1:6)]

# CPTAC-3
load("/projects/pangen/analysis/jtopham/data/cptac_org_rnaseq.RData")

# ICGC and TCGA
load("/projects/pangen/analysis/jtopham/data/ext_rnaseq.RData")

# combine all RNAseq data
all.rna <- all.rna[rownames(all.rna) %in% 
                     intersect(rownames(all.rna), rownames(ext)), ]
all.rna <- all.rna[rownames(all.rna) %in% 
                     intersect(rownames(all.rna), cptac$gene), ]

ext <- ext[rownames(ext) %in% intersect(rownames(all.rna), rownames(ext)), ]
ext <- ext[match(rownames(all.rna), rownames(ext)), ]
cptac <- cptac[cptac$gene %in% intersect(rownames(all.rna), cptac$gene), ]
cptac <- cptac[match(rownames(all.rna), cptac$gene), -1]

pan <- cbind(all.rna, ext, cptac)
load("/projects/pangen/analysis/jtopham/data/all_clin_feb2024.RData")
pan <- pan[, colnames(pan) %in% c(grep("POG|PAN", colnames(pan), value = T),
                                  clin$sample[clin$cohort %in% 
                                                c("POG/PanGen", "TCGA", "CPTAC")])]
batches <- data.frame(read.delim(paste0("/projects/pangen/analysis/jtopham/",
                                        "data/pog_batches_v2.tsv"),
                                 header = F, sep = '\t'), stringsAsFactors = F)
colnames(batches) <- c("id", "p2")
batches <- batches[batches$id %in% colnames(pan), ]

# mapper for external IDs
dmap <- read.delim(paste0("/projects/pangen/analysis/jtopham/data/",
                          "external_paad/specimen.tsv"), header = T, sep = '\t',
                   stringsAsFactors = F)
dmap <- dmap[!duplicated(dmap$icgc_donor_id), ]

# manually add batches for non-POG samples
batches <- rbind(batches, cbind(id = colnames(ext)[colnames(ext) %in% 
                                                     dmap$icgc_donor_id[dmap$project_code == "PAAD-US"]],
                                p2 = "TCGA"))
batches <- rbind(batches, cbind(id = colnames(cptac)[grep("C3", colnames(cptac))],
                                p2 = "CPTAC3"))

batches$p2 <- as.character(batches$p2)
batches$id <- as.character(batches$id)
batches <- batches[batches$id %in% colnames(pan), ]

# KRAS mutation status
load("/projects/pangen/analysis/jtopham/data/kras_mutType_hg38.RData")

# protein data
load("/projects/pangen/analysis/jtopham/data/protein_data_sep2021_norm.RData")
prot <- prot[, colnames(prot) %in% c("NumberPSM", batches$id)]
prot2 <- equalMedianNormalization(prot[, -1]) # median centre samples

# ------------------------------------------------------------------------------
# normalize RNAseq data
# ------------------------------------------------------------------------------
# remove genes with zero variance in any batch
for(b in unique(batches$p2)){
  if(sum(colnames(pan) %in% batches$id[batches$p2 == b]) == 0){next}
  vars <- apply(pan[, colnames(pan) %in%
                      batches$id[batches$p2 == b]], 1, var)
  vars <- vars[!is.na(vars)]
  pan <- pan[rownames(pan) %in% names(vars)[vars > 0], ]
}

tmp <- batches[batches$id %in% colnames(pan), ]
pan2 <- ComBat(dat = as.matrix(pan[, match(tmp$id, colnames(pan))]), 
               batch = tmp$p2, par.prior = TRUE)
pan2 <- data.frame(pan2, stringsAsFactors = F)
colnames(pan2) <- gsub("[.]", "-", colnames(pan2))
pan3 <- data.frame(t(pan2), stringsAsFactors = F)

# ------------------------------------------------------------------------------
# identify NET signature
# ------------------------------------------------------------------------------
netg <- read.delim("/projects/pangen/analysis/jtopham/data/net_genes.txt",
                   sep = "\t", stringsAsFactors = F, header = F)
netsig <- pan3[, colnames(pan3) %in% c("CD177", "OLFM4", "CCDC25", "ILK", netg$V1)]
netsig$sample <- rownames(netsig)
netsig$type <- "met"
netsig$type[netsig$sample %in% batches$id[batches$p2 %in% c("TCGA", "CPTAC3")]] <- "resect"
netsig$CD177_2 <- as.character(netsig$CD177 >= 0.75)
netsig$OLFM4_2 <- as.character(netsig$OLFM4 >= 2.5)

anno <- netsig[, colnames(netsig) %in% c("CD177", "OLFM4", "type", "CCDC25", "ILK", "CD177_2",
                                         "OLFM4_2")]

netsig2 <- as.matrix(t(netsig[, !colnames(netsig) %in% c("sample", "CD177", "OLFM4", "type",
                                                "CCDC25", "ILK", "CD177_2", "OLFM4_2")]))
netsig$CD177_2 <- NULL
netsig$OLFM4_2 <- NULL

# correlations b/w markers and each gene to overlay
tmp <- data.frame(t(netsig2), stringsAsFactors = F)
cors <- NULL
for(i in 1:4){
  for(j in 1:ncol(tmp)){
    cors <- rbind(cors, cbind(gene = colnames(tmp)[j], marker = colnames(anno)[i],
                              cor = cor(tmp[, j], anno[, i]), method = "spearman"))
  }
}
cors <- data.frame(cors, stringsAsFactors = F)
cors$cor <- as.numeric(cors$cor)
cors <- dcast(gene ~ marker, data = cors, value.var = c("cor"))

# z score
netsig2 <- data.frame(t(apply(netsig2, 1, FUN = function(x){
  return((x - mean(unlist(x))) / sd(unlist(x)))})),
  stringsAsFactors = F)

# PDAC-specific coexpresssed NET sig
ca <- function(x, k){
  return(cutree(hclust(x, method = "ward.D2"), k = k))
}

cc <- ConsensusClusterPlus(as.matrix(t(netsig2)), maxK = 8, reps = 100,
                           pItem = 0.8, pFeature = 1, writeTable = FALSE,
                           clusterAlg = "ca", distance = "euclidean",
                           seed = 123, plot = 'pdf')

# extract cc results for k
k=6
tmp <- cc[[k]]$consensusMatrix
colnames(tmp) <- rownames(netsig2)
rownames(tmp) <- rownames(netsig2)
tmp <- tmp[rev(cc[[k]]$consensusTree$order), cc[[k]]$consensusTree$order]

# gene cluster assignments
ganno <- data.frame(cc[[k]]$consensusClass, stringsAsFactors = F)
colnames(ganno)[1] <- "cluster"
ganno <- cbind(ganno, cors[match(rownames(ganno), cors$gene), -1])

# 10 x 7
pheatmap(tmp, cluster_rows = F, cluster_cols = F, annotation_col = ganno,
         color = purples_pal)

# 23 and 8 genes
netgcc <- rownames(ganno)[ganno$cluster == 1] # CD177 + ILK
netgcc2 <- rownames(ganno)[ganno$cluster == 4] # OLFM4 + ILK
gclust <- colnames(tmp)[colnames(tmp) %in% c(netgcc, netgcc2)]

tmp <- ganno; tmp$gene <- rownames(tmp); tmp <- melt(tmp, id.vars = c("gene", "cluster"))
tmp$variable <- factor(tmp$variable, levels = c("CD177", "OLFM4", "ILK", "CCDC25"))
tmp$cluster <- as.character(tmp$cluster)
tmp$anno <- tmp$cluster %in% c("1", "4")
tmp$cluster <- factor(tmp$cluster, levels = c("3", "2", "6", "5", "4", "1"))

# 7 x 5
p <- ggplot(tmp, aes(x = cluster, y = value, fill = anno))
p + geom_boxplot(outlier.shape = NA, aes(group = cluster)) + geom_jitter(width = 0.1) +
  facet_wrap(~ variable, scales = "free") +
  scale_fill_manual(values = c("#b3bad8", "#f28729")) +
  theme(panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(linetype = "dashed", color = "gray89"))

# anova
summary(aov(value ~ cluster, data = tmp[tmp$variable == "CD177", ])) # 0.71
summary(aov(value ~ cluster, data = tmp[tmp$variable == "OLFM4", ])) # 0.89
summary(aov(value ~ cluster, data = tmp[tmp$variable == "CCDC25", ])) # 0.000203
summary(aov(value ~ cluster, data = tmp[tmp$variable == "ILK", ])) # 0.017
pairwise.wilcox.test(tmp$value[tmp$variable=="ILK"], tmp$cluster[tmp$variable=="ILK"],
                     p.adjust.method = "BH")

wilcox.test(tmp$value[tmp$variable == "ILK" & tmp$cluster == "4"],
            tmp$value[tmp$variable == "ILK" & tmp$cluster %in% c("2", "3", "5", "6")])$p.value

colnames(netsig2) <- gsub("[.]", "-", colnames(netsig2))

# ------------------------------------------------------------------------------
# cluster samples on PDAC NET signature
# ------------------------------------------------------------------------------
netsig3 <- netsig2[rownames(netsig2) %in% c(netgcc, netgcc2), ]
cc <- ConsensusClusterPlus(as.matrix(netsig3), maxK = 8, reps = 100,
                           pItem = 0.8, pFeature = 1, writeTable = FALSE,
                           clusterAlg = "ca", distance = "euclidean",
                           seed = 123, plot = 'pdf')
# extract cc results for k
k=6
tmp <- cc[[k]]$consensusMatrix
colnames(tmp) <- colnames(netsig3)
rownames(tmp) <- colnames(netsig3)
tmp <- tmp[rev(cc[[k]]$consensusTree$order), cc[[k]]$consensusTree$order]

tt <- data.frame(cc[[k]]$consensusClass, stringsAsFactors = F)
colnames(tt)[1] <- "cluster"
anno$clust <- tt$cluster[match(rownames(anno), rownames(tt))]

# 10 x 7
pheatmap(tmp,
         cluster_rows = F, cluster_cols = F, show_colnames = F,
         annotation_col = anno, color = greys_pal)

ganno2 <- data.frame(cbind(gene = c(netgcc, netgcc2),
                           group = c(rep("netgcc", length(netgcc)),
                                     rep("netgcc2", length(netgcc2)))),
                     stringsAsFactors = F)
rownames(ganno2) <- ganno2$gene; ganno2$gene <- NULL

netsig4 <- netsig3[match(gclust, rownames(netsig3)), match(colnames(tmp), colnames(netsig3))]

# 10 x 7
pheatmap(netsig4, annotation_row = ganno2,
         cluster_rows = F, cluster_cols = F, show_colnames = F,
         annotation_col = anno, color = blue_pal)

tmp <- anno
tmp$sample <- rownames(tmp); tmp <- melt(tmp, id.vars = c("sample", "clust", "type"))
tmp <- tmp[!grepl("_2", tmp$variable), ]
tmp$value <- as.numeric(tmp$value); tmp$clust <- as.character(tmp$clust)
tmp$clust2 <- tmp$clust; tmp$clust2[tmp$clust %in% c("2", "5", "6")] <- "other"
tmp$clust2 <- factor(tmp$clust2, levels = c("other", "3", "1", "4"))
tmp$variable <- factor(tmp$variable, levels = c("CD177", "OLFM4", "ILK", "CCDC25"))
tmp$anno <- tmp$clust %in% c("1", "3", "4")

# 7 x 5 
p <- ggplot(tmp, aes(x = clust2, y = value, fill = anno))
p + geom_boxplot(outlier.shape = NA, aes(group = clust2)) + geom_jitter(width = 0.1) +
  facet_wrap(~ variable, scales = "free")+
  scale_fill_manual(values = c("#b3bad8", "#f28729")) +
  theme(panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(linetype = "dashed", color = "gray89"))

# anova
summary(aov(value ~ clust2, data = tmp[tmp$variable == "CD177", ])) # 0.00335
summary(aov(value ~ clust2, data = tmp[tmp$variable == "OLFM4", ])) # 0.00286
summary(aov(value ~ clust2, data = tmp[tmp$variable == "CCDC25", ])) # 0.177
summary(aov(value ~ clust2, data = tmp[tmp$variable == "ILK", ])) # 2.17e-11
pairwise.wilcox.test(tmp$value[tmp$variable=="CD177"], tmp$clust2[tmp$variable=="CD177"],
                     p.adjust.method = "BH")

# ------------------------------------------------------------------------------
# clinical correlates
# ------------------------------------------------------------------------------
tmpclin <- clin[clin$sample %in% rownames(anno), ]
tmpclin$group <- anno$clust[match(tmpclin$sample, rownames(anno))]
tmpclin <- tmpclin[tmpclin$group %in% c("1", "3", "4"), ]

p1 <- ggsurv(survfit(Surv(os, censor) ~ group, 
               data = tmpclin[tmpclin$cohort=="POG/PanGen", ]),
       surv.col = c("#39c6d2", "#545454", "#7289da"),
       size.est = 1.5, cens.size = 3.5, back.white = TRUE,
       order.legend = FALSE, xlab = "time (months)")
p2 <- ggsurv(survfit(Surv(os, censor) ~ group, 
               data = tmpclin[tmpclin$cohort=="TCGA", ]),
       surv.col = c("#39c6d2", "#545454", "#7289da"),
       size.est = 1.5, cens.size = 3.5, back.white = TRUE,
       order.legend = FALSE, xlab = "time (months)")
p3 <- ggsurv(survfit(Surv(os, censor) ~ group, 
               data = tmpclin[tmpclin$cohort=="CPTAC", ]),
       surv.col = c("#39c6d2", "#545454", "#7289da"),
       size.est = 1.5, cens.size = 3.5, back.white = TRUE,
       order.legend = FALSE, xlab = "time (months)")

# 15 x 4
plot_grid(p1, p3, p2, ncol = 3)

# log-rank p values
1 - pchisq(survdiff(Surv(os, censor) ~ group, 
                    data = tmpclin[tmpclin$cohort=="POG/PanGen" &
                                     tmpclin$group %in% c(1, 3, 4), ])$chisq,
           length(survdiff(Surv(os, censor) ~ group,
                           data = tmpclin[tmpclin$cohort=="POG/PanGen" &
                                            tmpclin$group %in% c(1, 3, 4), ])$n) - 1)

# HR and CI
res.cox <- coxph(Surv(os, censor) ~ group, 
                 data = tmpclin[tmpclin$cohort=="POG/PanGen" &
                                  tmpclin$group %in% c(1, 3, 4), ])
summary(res.cox)

# other vs NETsig high
tmpclin2 <- clin[clin$sample %in% rownames(anno), ]
tmpclin2$group <- anno$clust[match(tmpclin2$sample, rownames(anno))]
tmpclin2$group <- tmpclin2$group %in% c("1", "3", "4")

p1 <- ggsurv(survfit(Surv(os, censor) ~ group, 
                     data = tmpclin2[tmpclin2$cohort=="POG/PanGen", ]),
             #surv.col = c("#39c6d2", "#545454", "#7289da"),
             size.est = 1.5, cens.size = 3.5, back.white = TRUE,
             order.legend = FALSE, xlab = "time (months)")
p2 <- ggsurv(survfit(Surv(os, censor) ~ group, 
                     data = tmpclin2[tmpclin2$cohort=="TCGA", ]),
             #surv.col = c("#39c6d2", "#545454", "#7289da"),
             size.est = 1.5, cens.size = 3.5, back.white = TRUE,
             order.legend = FALSE, xlab = "time (months)")
p3 <- ggsurv(survfit(Surv(os, censor) ~ group, 
                     data = tmpclin2[tmpclin2$cohort=="CPTAC", ]),
             #surv.col = c("#39c6d2", "#545454", "#7289da"),
             size.est = 1.5, cens.size = 3.5, back.white = TRUE,
             order.legend = FALSE, xlab = "time (months)")
# 15 x 4
plot_grid(p1, p3, p2, ncol = 3)

# log-rank p values
1 - pchisq(survdiff(Surv(os, censor) ~ group, 
                    data = tmpclin2[tmpclin2$cohort=="POG/PanGen", ])$chisq,
           length(survdiff(Surv(os, censor) ~ group,
                           data = tmpclin2[tmpclin2$cohort=="POG/PanGen", ])$n) - 1)

# HR and CI
res.cox <- coxph(Surv(os, censor) ~ group, 
                 data = tmpclin2[tmpclin2$cohort=="POG/PanGen", ])
summary(res.cox)

# show all groups
tmpclin <- clin[clin$sample %in% rownames(anno), ]
tmpclin$group <- anno$clust[match(tmpclin$sample, rownames(anno))]
tmpclin$group[!tmpclin$group %in% c("1", "3", "4")] <- "other"

p1 <- ggsurv(survfit(Surv(os, censor) ~ group, 
                     data = tmpclin[tmpclin$cohort=="POG/PanGen", ]),
             #surv.col = c("#39c6d2", "#545454", "#7289da"),
             size.est = 1.5, cens.size = 3.5, back.white = TRUE,
             order.legend = FALSE, xlab = "time (months)")
p2 <- ggsurv(survfit(Surv(os, censor) ~ group, 
                     data = tmpclin[tmpclin$cohort=="TCGA", ]),
             #surv.col = c("#39c6d2", "#545454", "#7289da"),
             size.est = 1.5, cens.size = 3.5, back.white = TRUE,
             order.legend = FALSE, xlab = "time (months)")
p3 <- ggsurv(survfit(Surv(os, censor) ~ group, 
                     data = tmpclin[tmpclin$cohort=="CPTAC", ]),
             #surv.col = c("#39c6d2", "#545454", "#7289da"),
             size.est = 1.5, cens.size = 3.5, back.white = TRUE,
             order.legend = FALSE, xlab = "time (months)")
# 15 x 4
plot_grid(p1, p3, p2, ncol = 3)

# ------------------------------------------------------------------------------
# Aug 2024 update : add seven additional genes
# ------------------------------------------------------------------------------
tmp <- pan3[, colnames(pan3) %in% c("ACTN1", "ACTN4", "ITGB1", "MTOR",
                                    "AKT2", "CD44")]
tmp$sample <- rownames(tmp); tmp <- melt(tmp, id.vars = c("sample"))
tmp$variable <- factor(tmp$variable, levels = c("ITGB1", "ACTN1", "ACTN4", 
                                                "CD44", "MTOR", "AKT2"))
tmp$cluster <- anno$clust[match(tmp$sample, rownames(anno))]
tmp$clust2 <- tmp$cluster; tmp$clust2[tmp$cluster %in% c("2", "5", "6")] <- "other"
tmp$clust2 <- factor(tmp$clust2, levels = c("other", "3", "1", "4"))
tmp$anno <- tmp$cluster %in% c("1", "3", "4")

# 8 x 6
p <- ggplot(tmp, aes(x = clust2, y = value, fill = clust2))
p + geom_boxplot(outlier.shape = NA, aes(group = clust2)) + geom_jitter(width = 0.1) +
  facet_wrap(~ variable, scales = "free")+
  scale_fill_manual(values = c("#dfdfde", "#86ad73", "#b15827", "#ccbea3")) +
  theme(panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(linetype = "dashed", color = "gray89"))

# anova
summary(aov(value ~ clust2, data = tmp[tmp$variable == "ACTN1", ]))
pairwise.wilcox.test(tmp$value[tmp$variable=="ACTN1"], tmp$clust2[tmp$variable=="ACTN1"],
                     p.adjust.method = "BH")

# additional genes
prot3 <- prot2[rownames(prot2) %in% c("ACTN1", "ACTN4", "ACTB", "MTOR",
                                      "AKT1", "AKT2", "CD44"), ]
prot3$gene <- rownames(prot3); prot3 <- melt(prot3, id.vars = c("gene"))
prot3$clust <- anno$clust[match(prot3$variable, rownames(anno))]
prot3$clust2 <- prot3$clust; prot3$clust2[prot3$clust %in% c("2", "5", "6")] <- "other"
prot3$clust2 <- factor(prot3$clust2, levels = c("other", "3", "1", "4"))
prot3$gene <- factor(prot3$gene, levels = c("ACTN1", "ACTN4", "ACTB", "MTOR",
                                            "AKT1", "AKT2", "CD44"))
prot3$anno <- prot3$clust %in% c("1", "3", "4")

# 8 x 6
p <- ggplot(prot3, aes(x = clust2, y = value, fill = anno))
p + geom_boxplot(outlier.shape = NA, aes(group = clust2)) + geom_jitter(width = 0.1) +
  facet_wrap(~ gene, scales = "free")+
  scale_fill_manual(values = c("#b3bad8", "#f28729")) +
  theme(panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(linetype = "dashed", color = "gray89"))

# anova
summary(aov(value ~ clust2, data = prot3[prot3$gene  == "ACTN1", ]))

pairwise.wilcox.test(prot3$value[prot3$gene=="CD44"], prot3$clust2[prot3$gene=="CD44"],
                     p.adjust.method = "BH")

# ------------------------------------------------------------------------------
# Aug 2024 update : table showing prop. of clusters 1/2 that 
# have 'upreg' of the 7 genes
# ------------------------------------------------------------------------------
tmp <- pan3[, colnames(pan3) %in% c("ACTN1", "ACTN4", "ACTB", "MTOR",
                                    "AKT1", "AKT2", "CD44", "ITGB1")]
tmp$sample <- rownames(tmp); tmp <- melt(tmp, id.vars = c("sample"))
tmp$cluster <- anno$clust[match(tmp$sample, rownames(anno))]
tmp$cohort <- batches$p2[match(tmp$sample, batches$id)]
tmp$cohort[tmp$cohort %in% c("2.x", "3.x", "4.x")] <- "met"
tmp$cohort[tmp$cohort != "met"] <- "resected"
tmp$str <- paste0(tmp$variable, "__", tmp$cohort)
tmp <- split(tmp, f = tmp$str)
tmp <- lapply(tmp, FUN=function(x){
  x$up <- x$value >= quantile(x$value, 0.75)
  return(x)
})
tmp <- do.call(rbind, tmp)
tmp$cluster <- as.character(tmp$cluster)
tmp <- tmp[tmp$cluster %in% c("3", "1"), ]
tmp$str <- paste0(tmp$str, "_", tmp$cluster)
tmp <- split(tmp, f = tmp$str)
tmp2 <- lapply(tmp, FUN=function(x){
  return(cbind(group = unique(x$str),
               num = sum(x$up),
               perc = 100 * (sum(x$up) / nrow(x))
  ))
})
tmp2 <- do.call(rbind, tmp2)
tmp2 <- data.frame(tmp2, stringsAsFactors = F)
tmp2$num <- as.numeric(tmp2$num); tmp2$perc <- as.numeric(tmp2$perc)
tmp2$gene <- gsub("__.*", "", tmp2$group)
tmp2$group2 <- gsub(".*__", "", tmp2$group)
tmp2$anno <- paste0(tmp2$num, " (", round(tmp2$perc), "%)")
tmp2 <- dcast(gene ~ group2, data = tmp2, value.var = c("anno"))
tmp2 <- tmp2[, c(1, 5, 4, 3, 2)]

# Oct 2024 update : table showing prop. of 
# clusters 1/2 that have 'upreg' of any cluster-genes
tmp <- pan3[, colnames(pan3) %in% c(netgcc, netgcc2)]
tmp$sample <- rownames(tmp); tmp <- melt(tmp, id.vars = c("sample"))
tmp$cluster <- anno$clust[match(tmp$sample, rownames(anno))]
tmp$cohort <- batches$p2[match(tmp$sample, batches$id)]
tmp$cohort[tmp$cohort %in% c("2.x", "3.x", "4.x")] <- "met"
tmp$cohort[tmp$cohort != "met"] <- "resected"
tmp$str <- paste0(tmp$variable, "__", tmp$cohort)
tmp <- split(tmp, f = tmp$str)
tmp <- lapply(tmp, FUN=function(x){
  x$up <- x$value >= quantile(x$value, 0.75)
  return(x)
})
tmp <- do.call(rbind, tmp)
tmp$cluster <- as.character(tmp$cluster)
tmp <- tmp[tmp$cluster %in% c("3", "1"), ]
tmp$str <- paste0(tmp$str, "_", tmp$cluster)
tmp <- split(tmp, f = tmp$str)
tmp2 <- lapply(tmp, FUN=function(x){
  return(cbind(group = unique(x$str),
               num = sum(x$up),
               perc = 100 * (sum(x$up) / nrow(x))
  ))
})
tmp2 <- do.call(rbind, tmp2)
tmp2 <- data.frame(tmp2, stringsAsFactors = F)
tmp2$num <- as.numeric(tmp2$num); tmp2$perc <- as.numeric(tmp2$perc)
tmp2$gene <- gsub("__.*", "", tmp2$group)
tmp2$group2 <- gsub(".*__", "", tmp2$group)
tmp2$anno <- paste0(tmp2$num, " (", round(tmp2$perc), "%)")
tmp2 <- dcast(gene ~ group2, data = tmp2, value.var = c("anno"))
tmp2 <- tmp2[, c(1, 5, 4, 3, 2)]

# ------------------------------------------------------------------------------
# protein data analysis
# ------------------------------------------------------------------------------
prot3 <- prot2[rownames(prot2) %in% c("ITGB1", "CCDC25", "ILK", "OLFM4", "CD177", "CD44"), ]
prot3$gene <- rownames(prot3); prot3 <- melt(prot3, id.vars = c("gene"))
prot3$clust <- anno$clust[match(prot3$variable, rownames(anno))]
prot3$clust2 <- prot3$clust; prot3$clust2[prot3$clust %in% c("2", "5", "6")] <- "other"
prot3$clust2 <- factor(prot3$clust2, levels = c("other", "3", "1", "4"))
prot3$gene <- factor(prot3$gene, levels = c("CD177", "OLFM4", "ILK", "CCDC25", "CD44", "ITGB1"))
prot3$anno <- prot3$clust %in% c("1", "3", "4")

# 8 x 6
p <- ggplot(prot3, aes(x = clust2, y = value, fill = clust2))
p + geom_boxplot(outlier.shape = NA, aes(group = clust2)) + geom_jitter(width = 0.1) +
  facet_wrap(~ gene, scales = "free")+
  scale_fill_manual(values = c("#dfdfde", "#86ad73", "#b15827", "#ccbea3")) +
  theme(panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(linetype = "dashed", color = "gray89"))

# anova
summary(aov(value ~ clust2, data = prot3[prot3$gene == "ITGB1", ]))
pairwise.wilcox.test(prot3$value[prot3$gene=="ITGB1"], prot3$clust2[prot3$gene=="ITGB1"],
                     p.adjust.method = "BH")

# ITGB1-specific
tmp <- pan3; tmp$sample <- rownames(tmp); tmp <- tmp[, colnames(tmp) %in% c("sample", "ITGB1")]
colnames(tmp)[1] <- "rna"
tmp$clust <- anno$clust[match(tmp$sample, rownames(anno))]
tmp$clust2 <- tmp$clust; tmp$clust2[tmp$clust %in% c("2", "5", "6")] <- "other"
tmp$clust2 <- factor(tmp$clust2, levels = c("other", "3", "1", "4"))
tmp$anno <- tmp$clust %in% c("1", "3", "4")

tt <- data.frame(t(prot2), stringsAsFactors = F)
tmp$protein <- tt$ITGB1[match(tmp$sample, rownames(tt))]
tmp <- melt(tmp, measure.vars = c("rna", "protein"))

# 9 x 4
p <- ggplot(tmp, aes(x = clust2, y = value, fill = anno))
p + geom_boxplot(outlier.shape = NA, aes(group = clust2)) + geom_jitter(width = 0.1) +
  facet_wrap(~ variable, scales = "free")+
  scale_fill_manual(values = c("#b3bad8", "#f28729")) +
  theme(panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(linetype = "dashed", color = "gray89"))

summary(aov(value ~ clust2, data = tmp[tmp$variable == "rna", ])) # 1.16e-11
pairwise.wilcox.test(tmp$value[tmp$variable == "rna"], tmp$clust2[tmp$variable == "rna"],
                     p.adjust.method = "BH")

# ------------------------------------------------------------------------------
# DEA b/w cluster pts 1 and 3
# ------------------------------------------------------------------------------
dea <- NULL
for(i in 1:nrow(pan2)){
  dea <- rbind(dea, 
               cbind(gene = rownames(pan2)[i],
                     pval = wilcox.test(unlist(pan2[i, colnames(pan2) %in% 
                                                      rownames(anno)[anno$clust == "3"]]),
                                        unlist(pan2[i, colnames(pan2) %in% 
                                                      rownames(anno)[anno$clust == "1"]]))$p.value,
                     dir = mean(unlist(pan2[i, colnames(pan2) %in% rownames(anno)[anno$clust == "3"]])) -
                       mean(unlist(pan2[i, colnames(pan2) %in% rownames(anno)[anno$clust == "1"]]))))
}
dea <- data.frame(dea, stringsAsFactors = F)
dea$dir <- as.numeric(dea$dir)
dea$padj <- p.adjust(as.numeric(dea$pval), method = "BH")
dea <- dea[order(dea$padj), ]

# GSEA
gsea.pos <- my.GSEA.hyper(head(dea$gene[dea$dir > 0], 500), 
                          anno.genes[grep("^HALLMARK|^KEGG|^REACTOME", anno.genes$msig), ], 
                          dea$gene)
gsea.neg <- my.GSEA.hyper(head(dea$gene[dea$dir < 0], 500), 
                          anno.genes[grep("^HALLMARK|^KEGG|^REACTOME", anno.genes$msig), ], 
                          dea$gene)

# demonstrate consistency across iterative cutoffs
gseaip <- NULL; gseain <- NULL
for(i in seq(50, 1000, 50)){
  gseaip <- rbind(gseaip, 
                  cbind(thres = i, my.GSEA.hyper(head(dea$gene[dea$dir > 0], i), 
                                      anno.genes[grep("^HALLMARK|^KEGG|^REACTOME",
                                                      anno.genes$msig), ], 
                                      dea$gene)))
  gseain <- rbind(gseain, 
                  cbind(thres = i, my.GSEA.hyper(head(dea$gene[dea$dir < 0], i), 
                                      anno.genes[grep("^HALLMARK|^KEGG|^REACTOME",
                                                                 anno.genes$msig), ], 
                                      dea$gene)))
}
gseaip$neg10 <- -log10(gseaip$padj)

tmp <- aggregate(neg10 ~ msig.annotation, data = gseaip, mean)
tmp <- tmp[order(tmp$neg10, decreasing = T), ]

p <- ggplot(gseaip[gseaip$msig.annotation %in% head(tmp$msig.annotation), ], 
            aes(x = thres, y = neg10, color = msig.annotation))
p + geom_point(size = 2) + geom_path(aes(group = msig.annotation)) +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  theme(panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(linetype = "dashed", color = "gray89"))

# bar plot for top 10 gene sets
gsea <- gsea.pos[1:10, ]
gsea$neg10 <- -log10(gsea$padj)

gsea$msig.annotation <- factor(gsea$msig.annotation, levels = rev(gsea$msig.annotation))
gsea$num.overlap <- as.numeric(gsea$num.overlap)
gsea$num.input <- as.numeric(gsea$num.input)
gsea$perc <- round((gsea$num.overlap / gsea$num.input) * 100, 0)
gsea$anno <- paste0(gsea$num.overlap, " / ", gsea$num.input, " (", gsea$perc, "%)")

# 10 x 4
p <- ggplot(gsea, aes(x = msig.annotation, y = neg10, fill = neg10))
p + geom_bar(stat = "identity") + coord_flip() +
  geom_text(data = gsea, aes(x = msig.annotation, y = neg10 - 1, label = anno),
            color = "white") +
  scale_fill_gradientn(colors = ygob_pal[120:256]) +
  theme(panel.background = element_rect(fill = NA, color = "black"),
        panel.grid.major = element_line(color = "gray90",
                                        linetype = "dashed"),
        axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank(), legend.position = "")

# protein DEA
# make design table
prot3 <- prot2[, colnames(prot2) %in% rownames(anno)[anno$clust %in% c("3", "1")]]
cond <- colnames(prot3) %in% rownames(anno)[anno$clust == "3"]
cond[cond != TRUE] <- "group_1_goodp"; cond[cond == "TRUE"] <- "group_3_poorp"
cond <- factor(cond, levels = c("group_1_goodp", "group_3_poorp"))
design <- model.matrix(~ cond)
colnames(design) <- gsub("cond", "", colnames(design))
colnames(design)[1] <- "Intercept"

# make contrast
contrast <-  makeContrasts(contrasts = c("group_3_poorp"), levels = design)
fit1 <- lmFit(prot3, design)
fit2 <- contrasts.fit(fit1, contrasts = contrast)
fit3 <- eBayes(fit2)

# correct bias of variance estimate based on minimum number of psms 
# per protein used for quantification
psm.count.table <- data.frame(count = prot$NumberPSM, row.names =  rownames(prot))
fit3$count <- psm.count.table[rownames(fit3$coefficients), "count"]
fit4 <- spectraCounteBayes(fit3)

dpa <- outputResult(fit4, coef_col = 1)
dpa$neg10 <- -log10(dpa$sca.adj.pval)

# ------------------------------------------------------------------------------
# NLR data
# ------------------------------------------------------------------------------
# May 2024
nlr <- read.delim("/projects/sftp/jtopham/incoming/NLR_pangen.csv",
                  sep=",", stringsAsFactors = F)
nlr <- nlr[nlr$Event.Name == "Enrollment & BL Visit", ]
nlr <- nlr[nlr$Study.timepoint == "Baseline", ]
clin2 <- clin[!is.na(clin$pangen), ]
clin2 <- clin2[clin2$pangen %in% nlr$Study.ID, ]
clin2$nlr <- nlr$NLR[match(clin2$pangen, nlr$Study.ID)]
clin2$nlr2 <- clin2$nlr > 6.75

ggsurv(survfit(Surv(os, censor) ~ nlr2, 
               data = clin2),
       surv.col = c("#fda500", "#cccccc"),
       size.est = 1.5, cens.size = 3.5, back.white = TRUE,
       order.legend = FALSE, xlab = "time (months)")
# log-rank p value
1 - pchisq(survdiff(Surv(os, censor) ~ nlr2,
                    data = clin2)$chisq,
           length(survdiff(Surv(os, censor) ~ nlr2,
                           data = clin2)$n) - 1)
summary(coxph(Surv(os, censor) ~ nlr2, data = clin2))

# updated NLR (Oct 2024)
nlr <- read.delim("/projects/sftp/jtopham/incoming/NLR_oct2024.tsv",
                  sep="\t", stringsAsFactors = F)
clin2 <- clin[!is.na(clin$pangen), ]
clin2 <- clin2[clin2$pangen %in% nlr$Study.ID, ]
clin2$nlr <- nlr$NLR[match(clin2$pangen, nlr$Study.ID)]
clin2$nlr2 <- clin2$nlr > 6.75
# 6 x 5
ggsurv(survfit(Surv(os, censor) ~ nlr2, 
               data = clin2),
       surv.col = c("#b07697", "#322639"),
       size.est = 1.5, cens.size = 3.5, back.white = TRUE,
       order.legend = FALSE, xlab = "time (months)")
# log-rank p value
1 - pchisq(survdiff(Surv(os, censor) ~ nlr2,
                    data = clin2)$chisq,
           length(survdiff(Surv(os, censor) ~ nlr2,
                           data = clin2)$n) - 1)
# HR and CI
summary(coxph(Surv(os, censor) ~ nlr2, data = clin2))

# ------------------------------------------------------------------------------
# circos for EMT genes
# ------------------------------------------------------------------------------
cres <- dea[dea$gene %in% 
        unlist(anno.genes[anno.genes$msig == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", 2]), ]
cres <- cres[cres$dir > 0 & cres$padj < 0.05, ]
cres$neg10 <- -log10(cres$padj)
cres2 <- cres
cres2$pl2 <- dpa$logFC[match(cres2$gene, dpa$gene)]
cres2$pl2[is.na(cres2$pl2)] <- 0

# add dummy to scale y axes
tmp <- cbind(gene = "dummy", pval = NA, dir = 0.55, padj = NA, neg10 = 15,
             pl2 = 1.5)
cres2 <- data.frame(rbind(cres2, tmp), stringsAsFactors = F)
cres2$gene <- as.character(cres2$gene)
for(i in 2:ncol(cres2)){cres2[, i] <- as.numeric(cres2[, i])}
cres2 <- cres2[order(cres2$neg10, decreasing = T), ]

# space between groups
spaces <- seq(1, (nrow(cres2) * 10), 10)
spaces <- spaces + 10

t <- max(spaces) + 10
names(t) <- "gg"

this.row <- GRanges(seqnames = factor(names(t)),
                    ranges = IRanges(start = spaces, width = 10),
                    strand = '+',
                    gene = cres2$gene,
                    dir = cres2$dir, 
                    neg10 = cres2$neg10,
                    prot = cres2$pl2,
                    gridline = 15,
                    score = 1,
                    seqlengths = t)

p <- ggplot()
p <- p +
  # vertical grid lines
  layout_circle(this.row, geom = 'rect', size = .1, color = "gray80",
                aes(x = start, score = 15), grid = T, fill = NA,
                radius = 16, trackWidth = 15, space.skip = .075,
                grid.background = NA, grid.line = "white") +
  # inside centre, for scaling. invisible
  layout_circle(this.row, geom = 'bar', size = .1, color = "white",
                aes(x = start, score = score), grid = T, fill = "white",
                radius = 1, trackWidth = 1, space.skip = .075,
                grid.background = NA, grid.line = "white") + 
  # bars and dots
  layout_circle(this.row, geom = 'bar', size = .1, color = "white",
                aes(x = start, y = neg10), grid = T, fill = "gray32",
                radius = 16, trackWidth = 15, space.skip = .075,
                grid.background = NA, grid.line = "white") +
  layout_circle(this.row, geom = 'bar', size = .1, color = "white",
                aes(x = start, score = score, fill = dir), grid = T,
                radius = 31, trackWidth = 2, space.skip = .075,
                grid.background = NA, grid.line = "white") +
  scale_fill_gradientn(colors = c("#f6d6d6", "#45102a")) + 
  layout_circle(this.row, geom = 'point', grid = T, 
                  aes(x = start, y = prot, size = prot, color = prot),
                  radius = 16, trackWidth = 15, space.skip = .075,
                  grid.background = NA, grid.line = "gray40") +
  # dot outer circles
  layout_circle(this.row, geom = 'point', grid = T, color = "black",
                aes(x = start, y = prot, size = prot), shape = 1,
                radius = 16, trackWidth = 15, space.skip = .075,
                grid.background = NA, grid.line = "gray40") +
  scale_color_gradientn(colors = c("#d4eef2", "#43787f")) +
  # gene labels
  layout_circle(this.row, geom = 'text', size = 3.2, color = "black",
                aes(label = gene), angle = 90, grid = F,
                radius = 34, trackWidth = 1, space.skip = .075,
                grid.background = NA, grid.line = "white") 

# add ribbons
load("/projects/pangen/analysis/focus/genefriends/gf_pergene.RData")
gfdata <- gfdata[gfdata$gene %in% cres2$gene, ]
gfdata <- split(gfdata, f = gfdata$gene)
gfdata <- lapply(gfdata, FUN=function(x){
  x <- data.frame(cbind(gene = x[, 1], friend = unlist(strsplit(x[, 2], "; ")),
                        rank = unlist(strsplit(x[, 3], "; "))), stringsAsFactors = F)
  x[, 3] <- as.numeric(x[, 3])
  x <- x[x$friend %in% cres2$gene, ]
  x <- x[order(x$rank), ]; #x <- head(x, 3)
  return(x)}
)
gfdata <- do.call(rbind, gfdata)
rib <- gfdata
rib <- rbind(rib, cbind(gene = "dummy", friend = "SLIT2",
                        rank = 3))
rib$gene <- as.character(rib$gene); rib$friend <- as.character(rib$friend)
rib$rank <- as.numeric(rib$rank)

rib$to.gr <- paste0("gg:", spaces[match(rib$friend, cres2$gene)], 
                    "-", (spaces[match(rib$friend, cres2$gene)] + 9),
                    ":+")

this.rib <- GRanges(seqnames = factor(names(t)),
                    ranges = IRanges(start = spaces[match(rib$gene,
                                                          cres2$gene)], width = 10),
                    strand = "+",
                    to.gr = GRanges(rib$to.gr),
                    strength = rib$rank,
                    seqlengths = t)

# coloring: save p, then save this p+, and paste ribbons from the latter. 8x8
p + 
  layout_circle(this.rib, geom = "link", linked.to = "to.gr", size = 1.2, 
                aes(x = start, color = strength, alpha = 0.8),
                radius = 14, trackWidth = 5, space.skip = .075) +
  scale_color_gradientn(colors = rev(c("#f5eedf", "#6a4917", "black")))
  