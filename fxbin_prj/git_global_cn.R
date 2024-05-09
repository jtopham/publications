# ------------------------------------------------------------------------------
#                           _      _      _
#                        __(.)< __(.)> __(.)=
#                        \___)  \___)  \___)
# ------------------------------------------------------------------------------
# FILENAME : fxn_bins_prj.R
#
# AUTHOR : James Topham
#
# DATE : May 9, 2024
#
# DESCRIPTION :
#
#   -  global copy number analysis in PDAC
#   -  patient identifiers have been censored (xxx)
# ------------------------------------------------------------------------------
#                           _      _      _
#                        __(.)< __(.)> __(.)=
#                        \___)  \___)  \___)
# ------------------------------------------------------------------------------
library(ggplot2)
library(reshape2)
library(survival)
library(sva)
library(ConsensusClusterPlus)
library(pheatmap)
library(GGally)
library(DESeq2)
library(DEqMS)
library(cowplot)

# PanGen cases (n=69)
samples <- read.delim("/projects/pangen/analysis/jtopham/data/seq_pangen.txt",
                      header = F, sep = '\t', stringsAsFactors = F)$V1
# remove PNET cases
samples <- samples[!samples %in% c("xxx", "xxx")]

# ensure sample mapping
load("/projects/pangen/analysis/jtopham/data/global_cn_segs.RData")
samples <- samples[samples %in% segs$sample]

# color palettes
source("/projects/pangen/analysis/jtopham/scripts/custom_pal.R")

# function for performing PCA to cluster samples
source("/projects/pangen/analysis/jtopham/scripts/my_pca.R")

# function for sorting columns in oncoprints
scoreCol <- function(x) {
  score <- 0;
  for(i in 1:length(x)) {
    if(x[i]==1) {score <- score + 2^(length(x)-i)}
  }
  return(score);
}

# gene sets and GSEA script
load("/projects/pangen/analysis/jtopham/data/msigdb_apr2021.RData")
source("/projects/pangen/analysis/jtopham/scripts/my_gsea.R")

# mapper for external IDs
dmap <- read.delim(paste0("/projects/pangen/analysis/jtopham/data/external_paad/",
                          "specimen.tsv"), header = T, sep = '\t',
                   stringsAsFactors = F)
dmap <- dmap[!duplicated(dmap$icgc_donor_id), ]

# ------------------------------------------------------------------------------
# data
# ------------------------------------------------------------------------------
# RNA-seq
load("/projects/pangen/analysis/jtopham/data/all_rna_oct2022.RData")
rownames(all.rna) <- all.rna$geneid
all.rna <- all.rna[, -c(1:6)]
all.rna <- all.rna[, colnames(all.rna) %in% samples]

# raw rnaseq for use with DESeq2
load("/projects/pangen/analysis/jtopham/data/raw_rna_aug2023.RData")
rownames(raw.rna) <- raw.rna$geneid
raw.rna <- raw.rna[, -c(1:6)]
raw.rna <- raw.rna[, colnames(raw.rna) %in% samples]
raw.rna <- raw.rna[!apply(raw.rna, 1, FUN=function(x){return(any(is.na(x)))}), ]

# tumor content values
load("/projects/pangen/analysis/jtopham/data/all_pdac_tcs_nov2021.RData")
tc <- tc[tc$id %in% samples, ]

# clinical metadata
load("/projects/pangen/analysis/jtopham/data/all_clin_may2023.RData")

# PDAC molecular subtypes
load("/projects/pangen/analysis/jtopham/data/all_pdac_subtypes_aug2022.RData")
pur <- pur[pur$sample %in% samples, ]

# mutation data
load("/projects/pangen/analysis/jtopham/data/all_pdac_snv_oct2022.RData")
snv <- snv[snv$sample %in% samples, ]

load("/projects/pangen/analysis/jtopham/data/all_pdac_cnv_oct2022.RData")
rownames(cnv) <- cnv$gene
cnv <- cnv[, colnames(cnv) %in% samples]

load("/projects/pangen/analysis/jtopham/data/kras_mutType.RData")
kt <- kt[kt$sample %in% samples, ]

# fusion data
load("/projects/pangen/analysis/jtopham/data/rna_fusions_jan2023.RData")
fus <- fus[fus$sample %in% samples, ]

# TMB
load("/projects/pangen/analysis/jtopham/data/tmb.RData")
tmb <- tmb[tmb$sample %in% samples, ]

# HRD
load("/projects/pangen/analysis/jtopham/data//hrd_nov2021.RData")
hrd <- hrd[hrd$sample %in% samples, ]

# RNAseq batches
batches <- data.frame(read.delim(paste0("/projects/pangen/analysis/jtopham/",
                                        "data/pog_batches_v2.tsv"),
                                 header = F, sep = '\t'), stringsAsFactors = F)
colnames(batches) <- c("id", "p2")
batches <- batches[batches$id %in% colnames(all.rna), ]
batches$p2 <- as.character(batches$p2); batches$id <- as.character(batches$id)

# resectable datasets
load("/projects/pangen/analysis/jtopham/data/ext_rnaseq.RData")

tmp <- cnv
load("/projects/pangen/analysis/jtopham/data/all_pdac_cnv_27APR20.RData")
rownames(cnv) <- cnv$gene
cnv <- cnv[, -c(1:5)]
cnv <- cnv[, colnames(cnv) %in% colnames(ext)]
ext <- ext[, colnames(ext) %in% colnames(cnv)]
rcnv <- cnv; cnv <- tmp

tmp <- snv
load("/projects/pangen/analysis/jtopham/data/all_pdac_snv_oct2022.RData")
rsnv <- snv[snv$sample %in% colnames(ext), ]
snv <- tmp

# ------------------------------------------------------------------------------
# functional (fxn) bin data
# ------------------------------------------------------------------------------
files <- paste0("/projects/pangen/analysis/jtopham/samples/", samples, 
                "/cnv/facets_stdout.txt")
pog.p <- NULL
for(i in 1:length(files)){
  if(is.na(file.info(files[i])$size) |
     file.info(files[i])$size == 0){next}
  
  sample <- strsplit(files[i], "/")[[1]][7]
  tmp <- read.delim(files[i], sep = '\t', header = F, stringsAsFactors = F)
  pog.p <- rbind(pog.p, cbind(sample = sample,
                              ploid = gsub(".*: ", "", tmp[4, ])))
}
pog.p <- data.frame(pog.p, stringsAsFactors = F)
pog.p$ploid2 <- round(as.numeric(pog.p$ploid))

if(file.exists("/projects/pangen/analysis/jtopham/data/global_cn_segs.RData")){
  load("/projects/pangen/analysis/jtopham/data/global_cn_segs.RData")
  load("/projects/pangen/analysis/jtopham/data/global_cn_segs2.RData")
}else{
  files <- paste0("/projects/pangen/analysis/jtopham/samples/", samples, 
                  "/cnv/hg19_all_cnvbins.bed")
  segs <- NULL
  for(i in files){
    tmp <- read.delim(i, sep = '\t', stringsAsFactors = F, header = F)
    segs <- rbind(segs, cbind(sample = strsplit(i, "/")[[1]][7], tmp))
  }
  
  segs$ploid <- pog.p$ploid2[match(segs$sample, pog.p$sample)]
  segs$cnv <- "neutral"
  segs$cnv[grep("NA|[.]", paste0(segs$V17, segs$V18))] <- "neutral" # NA
  segs$cnv[segs$V18 == "0"] <- "het loss"
  segs$cnv[segs$V17 == "0" & segs$V18 == "0"] <- "hom loss"
  # gains defined as being at least twice the ploidy (ie POG == 4)
  segs$cnv[as.numeric(segs$V17) >= segs$ploid* 2] <- "gain"
  segs$cnv[as.numeric(segs$V17) >= segs$ploid* 5] <- "gain_5"
  segs$cnv[as.numeric(segs$V17) >= segs$ploid* 10] <- "gain_10"
  
  segs <- segs[order(segs$V1, segs$V2), ]
  
  segs2 <- dcast(V4 ~ sample, data = segs, value.var = c("cnv"))
  segs2 <- segs2[match(unique(segs$V4), segs2$V4), ]
  rownames(segs2) <- segs2$V4; segs2$V4 <- NULL
  save(segs, file = "/projects/pangen/analysis/jtopham/data/global_cn_segs.RData")
  save(segs2, file = "/projects/pangen/analysis/jtopham/data/global_cn_segs2.RData")
}

# centromeric regions
cent <- read.delim("/projects/pangen/analysis/jtopham/data/hg19_cytobands.tsv",
                   sep = '\t', stringsAsFactors = F, header = F)
cent <- cent[cent$V5 == "acen", ]
cent$str <- paste0(cent$V1, ":", cent$V2, "_", cent$V3)
cent <- cent[!cent$V1 %in% c("X", "Y"), ]

# ------------------------------------------------------------------------------
# whole genome overview
# ------------------------------------------------------------------------------
# function to identify if coordinates (chr:start_stop) overlap with 
# any centromere regions
overlap <- function(x, y){
  # on same chr and one coord of x is within y
  if(any(((gsub(":.*", "", x) == gsub(":.*", "", y)) &
     (((as.numeric(gsub(".*:", "", gsub("_.*", "", x))) >
      as.numeric(gsub(".*:", "", gsub("_.*", "", y)))) &
      (as.numeric(gsub(".*:", "", gsub("_.*", "", x))) <
       as.numeric(gsub(".*_", "", y)))) |
      ((as.numeric(gsub(".*_", "", x))) >
        as.numeric(gsub(".*:", "", gsub("_.*", "", y)))) &
       (as.numeric(gsub(".*_", "", x))) <
        as.numeric(gsub(".*_", "", y)))) |
     # on same chr and one coord of y is within x
     ((gsub(":.*", "", y) == gsub(":.*", "", x)) &
      (((as.numeric(gsub(".*:", "", gsub("_.*", "", y))) >
         as.numeric(gsub(".*:", "", gsub("_.*", "", x)))) &
        (as.numeric(gsub(".*:", "", gsub("_.*", "", y))) <
         as.numeric(gsub(".*_", "", x)))) |
       ((as.numeric(gsub(".*_", "", y))) >
        as.numeric(gsub(".*:", "", gsub("_.*", "", x)))) &
       (as.numeric(gsub(".*_", "", y))) <
       as.numeric(gsub(".*_", "", x)))))){return(TRUE)
  }else{return(FALSE)}
}

# loop to walk along genome and aggregate invariable regions, while skipping
# centromeres
if(file.exists("/projects/pangen/analysis/jtopham/data/global_cn_segs3.RData")){
  load("/projects/pangen/analysis/jtopham/data/global_cn_segs3.RData")
}else{
  segs3 <- NA
  for(i in 1:nrow(segs2)){
    # if this is first entry, start first walk
    if(i == 1){
      segs3 <- rbind(segs3, segs2[i, ])
      next
    }
    # if this entry and previous entry are in centromeres, update start
    # coord and continue walking
    if(overlap(rownames(segs2[i-1, ]), cent$str)){
      if(overlap(rownames(segs2[i, ]), cent$str)){
        stop <- gsub(".*_", "", rownames(segs2)[i])
        rownames(segs3)[nrow(segs3)] <- paste0(gsub("_.*", "", rownames(segs3)[nrow(segs3)]), "_", stop)
      }else{
        # if previous entry within centromere and new entry is not, start new walk
        segs3 <- rbind(segs3, segs2[i, ])
        next
      }
      next
    }
    # if this entry overlaps a centromere, start new walk
    if(overlap(rownames(segs2[i, ]), cent$str)){
      segs3 <- rbind(segs3, segs2[i, ])
      next
    }
    # if this entry is different chr than latest entry, or non identical entry,
    # start new walk
    if(gsub(":.*", "", rownames(segs2)[i]) != gsub(":.*", "", rownames(segs2)[i-1])){
      segs3 <- rbind(segs3, segs2[i, ])
      next
    }
    if(!(all(segs2[i, ] == segs3[nrow(segs3), ]))){
      segs3 <- rbind(segs3, segs2[i, ])
      next
    }
    # if new entry is identical to latest entry, update stop coord and continue walking
    if(all(segs2[i, ] == segs3[nrow(segs3), ])){
      stop <- gsub(".*_", "", rownames(segs2)[i])
      rownames(segs3)[nrow(segs3)] <- paste0(gsub("_.*", "", rownames(segs3)[nrow(segs3)]), "_", stop)
      next
    }
  }
  segs3 <- segs3[-1, ]
  save(segs3, file = "/projects/pangen/analysis/jtopham/data/global_cn_segs3.RData")
}

# digitize
segs4 <- segs3
segs4[segs4 == "neutral"] <- "0"
segs4[segs4 == "het loss"] <- "-1"
segs4[segs4 == "hom loss"] <- "-2"
segs4[segs4 == "gain"] <- "3"
segs4[segs4 == "gain_5"] <- "3"
segs4[segs4 == "gain_10"] <- "3"
for(i in 1:ncol(segs4)){segs4[, i] <- as.numeric(segs4[, i])}

# disregard centromere regions
tmp <- NULL
for(i in rownames(segs4)){
  tmp <- c(tmp, overlap(i, cent$str))
}
for(i in 1:ncol(segs4)){segs4[tmp, i] <- 0}

anno.c <- segs4[, 1:2]
anno.c[, 1] <- gsub(":.*", "", rownames(anno.c))
anno.c[, 2] <- NULL; colnames(anno.c)[1] <- "chr"
anno.c$cent <- as.character(tmp)
anno.c$arm <- NA
anno.c <- split(anno.c, f = anno.c$chr)
anno.c <- lapply(anno.c, FUN = function(x){
  x$arm <- "q"
  x$arm[1:which(x$cent == TRUE)[1]] <- "p"
  x$arm[x$cent == TRUE] <- "cent"
  return(x)
})
anno.c <- do.call(rbind, anno.c)
rownames(anno.c) <- gsub(".*[.]", "", rownames(anno.c))

anno.r <- data.frame(t(segs4), stringsAsFactors = F)[, 1:2]
anno.r[, 1] <- as.character(rownames(anno.r) %in% c("xxx", "xxx", "xxx", "xxx",
                                                    "xxx", "xxx", "xxx", "xxx", 
                                                    "xxx"))
anno.r[, 2] <- NULL; colnames(anno.r)[1] <- "krasWT"
anno.r$ploid <- pog.p$ploid2[match(rownames(anno.r), pog.p$sample)]

clust <- hclust(dist(t(as.matrix(segs4[which(apply(segs4, 1,
         FUN=function(x){return(sum(x %in% c(3, -2))>2)})), ]))), method = "ward.D2")
segs4 <- segs4[, clust$order]

pheatmap(t(segs4), show_rownames = FALSE, show_colnames = FALSE, method = "ward.D2",
         cluster_cols = FALSE, annotation_col = anno.c, cluster_rows = FALSE,
         annotation_row = anno.r, color = c(echun_pal[c(18, 16)], "gray92", 
                                            "gray92", "gray92", "#bc4151"))

# number of events per bin
tmp <- apply(segs4, 1, FUN=function(x){return(sum(x == 3))})
tmp2 <- apply(segs4, 1, FUN=function(x){return(sum(x == -2))})
tmp <- data.frame(cbind(seg = rownames(segs4), gain = unname(tmp),
                        loss = unname(tmp2)), stringsAsFactors = F)
tmp$seg <- factor(tmp$seg, levels = rownames(segs4))
tmp$gain <- as.numeric(tmp$gain)
tmp$loss <- as.numeric(tmp$loss)
tmp <- melt(tmp, id.vars = c("seg"))

p <- ggplot(tmp, aes(x = seg, y = value, fill = variable))
p + geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("#bc4151", echun_pal[18])) +
  theme(axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(linetype = "dashed", color = "gray89"))

anno.r <- anno.r[match(colnames(segs4), rownames(anno.r)), ]
anno.r$group <- "3"
anno.r$group[1:25] <- "1"
anno.r$group[25:50] <- "2"
clin <- clin[clin$sample %in% rownames(anno.r), ]
clin$group <- anno.r$group[match(clin$sample, rownames(anno.r))]

# ------------------------------------------------------------------------------
# top hits and SNV density analyses
# ------------------------------------------------------------------------------
# genes with either hom loss or gain in at least 3 pts
tmp <- segs2
tmp[tmp == "gain_5"] <- "gain"
tmp[tmp == "gain_10"] <- "gain"
th <- apply(tmp, 1, table)
th <- lapply(th, FUN=function(x){
  return(data.frame(x))
})
th <- do.call(rbind, th)
th$Var1 <- as.character(th$Var1)

th <- th[th$Var1 %in% c("hom loss", "gain"), ]

th <- th[th$Freq >= 3, ]
th$seg <- gsub("[.].*", "", rownames(th))
th$chr <- gsub(":.*", "", th$seg)
th$mid <- as.numeric(gsub(".*:", "", gsub("_.*", "", th$seg))) + 10000
th <- th[order(th$Freq, decreasing = T), ]
# remove x instances where a seg has both gains+hom observed (keep most freq)
th <- th[!duplicated(th$seg ), ]

th <- th[order(th$chr, th$mid), ]

# group adjacent segments that have same Freq 
# start with unique groups
th$group <- 1:nrow(th)
th$done <- 0
# group nearby segs
for(i in 1:nrow(th)){
  # get segs on same chr and measure distance
  tmp <- th[th$chr == th$chr[i] & th$seg != th$seg[i], ]
  tmp$dist <- tmp$mid - th$mid[i]
  tmp <- tmp[order(abs(tmp$dist)), ]
  # only group segs within 1 mb
  tmp <- tmp[abs(tmp$dist) < 1e6, ]
  # if no nearby segs, keep it's unique group
  if(nrow(tmp) == 0){th$done[i] <- 1; next}
  # otherwise, set group to nearest neighbours group
  # give pref to those already iterated across
  if(any(tmp$done == 1)){tmp <- tmp[tmp$done == 1, ]}
  th$group[i] <- tmp$group[1]; th$done[i] <- 1
}
  
th <- th[order(th$Freq, decreasing = T), ]

gene <- read.delim("/projects/pangen/analysis/jtopham/data/allbins_vs_genes.tsv",
                   header = F, sep = '\t', stringsAsFactors = F)
gene$seg <- paste0(gene$V1, ":", gene$V2, "_", gene$V3)

th$gene <- gene$V14[match(th$seg, gene$seg)]
th2 <- th#[th$gene != ".", ]

# keep multi-gene info
# (segment "seg" IDs no longer meaningful)
th2$tmp <- paste0(th2$group, "_", th2$Freq)
th2 <- split(th2, f = th2$tmp)
th2 <- lapply(th2, FUN = function(x){
  x <- x[x$Freq == max(x$Freq), ]
  x$gene[1] <- paste(unique(x$gene), collapse = ",")
  x$aggstart <- min(x$mid)-10000
  x$aggend <- max(x$mid)+10000
  return(x[1, ])
})
th2 <- do.call(rbind, th2)
th2$tmp <- NULL

th2 <- th2[order(th2$Freq, decreasing = T), ]

# be 5 mb away from centromere
cent2 <- split(cent, f = cent$V1)
cent2 <- lapply(cent2, FUN=function(x){
  x$V2 <- min(x$V2); x$V3 <- max(x$V3)
  return(x[1, 1:3])
})
cent2 <- do.call(rbind, cent2)
th2$centstart <- cent2$V2[match(th2$chr, cent2$V1)]
th2$centend <- cent2$V3[match(th2$chr, cent2$V1)]
th2 <- th2[abs(th2$centend - th2$aggstart) > 5e6 &
           abs(th2$centstart - th2$aggend) > 5e6, ]

th2 <- th2[!duplicated(th2$group), ]

# Suppl table 1: Genomic bins with at least 3 samples with hom. loss or gains
# in the cohort. Bins with same mut. frequency that ar ewithin 1mb of each 
# other are merged
st1 <- th2[, c(4, 9:10, 1, 8, 2)]
st1$gene <- gsub(",[.]", "", st1$gene)
# write.table(st1, file = "/projects/sftp/jtopham/incoming/fxbin_paper_supptable1.tsv",
#             sep = '\t', quote = F, row.names = F)

th3 <- th2[1:20, ]
th3$perc <- (th3$Freq / ncol(segs3)) * 100
th3$Var1 <- factor(th3$Var1, levels = c("hom loss", "gain"))
th3 <- th3[order(th3$Var1, th3$Freq), ]
th3$seg <- factor(th3$seg, levels = th3$seg)

# 11 x 6
p <- ggplot(th3, aes(x = seg, y = perc, fill = Var1))
p + geom_bar(stat = "identity", size = 0.25, color = "white") + 
  coord_flip() + scale_fill_manual(values = c("#313695", "#bc4151")) +
  theme(axis.text = element_text(size = 8, color = "black"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(linetype = "dashed", color = "gray89")) +
  geom_text(data = th3, aes(x = seg, y = perc+15, label = gene))

# mini karyograms
chrsize <- read.delim("/projects/pangen/analysis/jtopham/data/hg19_chromsizes.tsv",
                      sep = '\t', header = F, stringsAsFactors = F)
th3$chrsize <- chrsize$V3[match(th3$chr, chrsize$V1)]
th3$centstartperc <- 100 * (th3$centstart / th3$chrsize)
th3$centendperc <- 100 * (th3$centend / th3$chrsize)

th3$aggstartperc <- 100 * (th3$aggstart / th3$chrsize)
th3$aggendperc <- 100 * (th3$aggend / th3$chrsize)

# 5 x 9
p <- ggplot(th3, aes(x = 100, y = seg))
p + geom_point() + geom_point(aes(x = 0, y = seg)) +
  geom_point(aes(x = centstartperc, y = seg), color = "red") +
  geom_point(aes(x = centendperc, y = seg), color = "red") +
  geom_point(aes(x = aggstartperc, y = seg), color = "#776ab9") +
  geom_point(aes(x = aggendperc, y = seg), color = "#776ab9") +
  theme(axis.text = element_text(size = 8, color = "black"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(linetype = "dashed", color = "gray89"))

# for each top-hit seg, aggregate # mutations for grouped-segs in each sample
# bins to # SNV
load("/projects/pangen/analysis/jtopham/data/allbins_to_snv.RData")

thp <- segs2; thp$seg <- rownames(thp); thp <- melt(thp, id.vars = c("seg"))
thp$str <- paste0(thp$variable, "_", thp$seg)
snvbin$str <- paste0(snvbin$sample, "_", snvbin$seg)
thp$ntotal <- snvbin$ntotal[match(thp$str, snvbin$str)]
thp$nlof <- snvbin$nLOF[match(thp$str, snvbin$str)]
thp$value[grepl("gain", thp$value)] <- "gain"
thp$anno <- thp$value
thp$anno[thp$nlof > 0] <- paste0(thp$value[thp$nlof > 0], " SNV")
thp$anno[thp$value == "hom loss"] <- "hom loss"

# hyper-mut
ord <- aggregate(ntotal ~ seg, data = thp, mean)
thp$mean_snv <- ord$ntotal[match(thp$seg, ord$seg)]
thp <- thp[order(thp$mean_snv, decreasing = T), ]
# top 5 most hyper-mut segs, save for later
boi <- head(unique(thp$seg), 6)
# remove one with stark outlier sample
boi <- boi[boi != "12:8360001_8380001"]

st2boi <- head(unique(thp$seg), 100)

# 2hits
thp$anno2 <- "other"
thp$anno2[grepl("het loss SNV|hom loss", thp$anno)] <- "2hit"
th4 <- data.frame(table(thp$seg, thp$anno2))
th4 <- th4[th4$Var2 == "2hit", ]
th4 <- th4[order(th4$Freq, decreasing = T), ]

th4$meansnv <- thp$mean_snv[match(th4$Var1, thp$seg)]
colnames(th4)[1] <- "seg"
th4$chr <- gsub(":.*", "", th4$seg)
th4$mid <- as.numeric(gsub(".*:", "", gsub("_.*", "", th4$seg))) + 10000
th4 <- th4[order(th4$chr, th4$mid), ]

# Suppl table 2: Top 100 highest SNV density bins
st2 <- th4[th4$seg %in% st2boi, -6]
st2 <- st2[order(st2$meansnv, decreasing = T), ]
st2$start <- gsub("_.*", "", gsub(".*:", "", st2$seg))
st2$end <- gsub(".*_", "", st2$seg)
st2$gene <- gene$V14[match(st2$seg, gene$seg)]
# write.table(st2[, c(5:8, 3, 4)], sep = '\t', quote = F, row.names = F,
#             file = "/projects/sftp/jtopham/incoming/fxbin_paper_supptable2.tsv")

th4 <- th4[th4$Freq >= 3 | th4$seg %in% boi, ]
th4$group <- 1:nrow(th4)
th4$done <- 0
th4$seg <- as.character(th4$seg)
# group nearby segs
for(i in 1:nrow(th4)){
  # get segs on same chr and measure distance
  tmp <- th4[th4$chr == th4$chr[i] & th4$seg != th4$seg[i], ]
  tmp$dist <- tmp$mid - th4$mid[i]
  tmp <- tmp[order(abs(tmp$dist)), ]
  # only group segs with4in 1 mb
  tmp <- tmp[abs(tmp$dist) < 1e6, ]
  # if no nearby segs, keep it's unique group
  if(nrow(tmp) == 0){th4$done[i] <- 1; next}
  # otherwise, set group to nearest neighbours group
  # give pref to those already iterated across
  if(any(tmp$done == 1)){tmp <- tmp[tmp$done == 1, ]}
  th4$group[i] <- tmp$group[1]; th4$done[i] <- 1
}

th4 <- th4[order(th4$Freq, decreasing = T), ]
th4$gene <- gene$V14[match(th4$seg, gene$seg)]
tmp <- th4[!duplicated(th4$group) | th4$seg %in% boi, ]
tmp <- tmp[(tmp$meansnv >= 0.13043478 & tmp$Freq >=5) | tmp$seg %in% boi, ]
tmp <- tmp[order(tmp$meansnv, decreasing = T), ]

thp2 <- thp[thp$seg %in% tmp$seg, ]
thp2$seg <- factor(thp2$seg, levels = tmp$seg)

# 10 x 5
p <- ggplot(thp2, aes(x = seg, y = ntotal))
p + geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  theme(axis.text = element_text(size = 8, color = "black"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(linetype = "dashed", color = "gray89"))

thp2$anno <- factor(thp2$anno, levels = rev(c("hom loss", "het loss SNV", "het loss", 
                                              "neutral SNV", "neutral", "gain", 
                                              "gain SNV")))

p <- ggplot(thp2, aes(x = seg, fill = anno))
p + geom_bar(position = "stack", color = "black", size = 0.25) +
  scale_fill_manual(values = rev(c("#003e86", "#9288c7", "#b2c5da", "#999999", 
                                   "gray92", "#c94c4c", "#6f1e1b"))) +
  theme(axis.text = element_text(size = 8, color = "black"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(linetype = "dashed", color = "gray89"))

# ------------------------------------------------------------------------------
# functionality analysis
# ------------------------------------------------------------------------------
# get # read centers (RNA) overlapping with each bin
if(file.exists("/projects/pangen/analysis/jtopham/data/fbin_rbins.RData")){
  load("/projects/pangen/analysis/jtopham/data/fbin_rbins.RData")
}else{
  files <- paste0("/projects/pangen/analysis/jtopham/samples/", samples, 
                  "/cnv/cnvbin_RNAcov.bed")
  rbins <- NULL
  for(i in files){
    if(!file.exists(i) | file.size(i)==0){next}
    tmp <- read.delim(i, sep = '\t', stringsAsFactors = F, header = F)
    if(nrow(tmp) != 144061){next}
    rbins <- cbind(rbins, tmp[, 5])
    colnames(rbins)[ncol(rbins)] <- strsplit(i, "/")[[1]][7]
  }
  rownames(rbins) <- tmp$V4
  rbins <- data.frame(rbins, stringsAsFactors = F)
  
  save(rbins, file = "/projects/pangen/analysis/jtopham/data/fbin_rbins.RData")
}

# normalize for depth and log-transform 
# R(center)PKM (log10); bins are 20kb
load("/projects/pangen/analysis/jtopham/data/bamstats.RData")
for(i in 1:ncol(rbins)){
  rbins[, i] <- log10(((rbins[, i] * 1e6 * 0.05) / 
                        stats$mapped[stats$lib == colnames(rbins)[i]]) + 1)
}

segs3 <- segs
segs3$cnv <- "neutral"
segs3$cnv[grep("NA|[.]", paste0(segs3$V17, segs3$V18))] <- "neutral"
segs3$cnv[segs3$V18 == "0"] <- "het loss"
segs3$cnv[segs3$V17 == "0" & segs3$V18 == "0"] <- "hom loss"
# gains defined as being at least twice the ploidy (ie POG == 4)
segs3$cnv[as.numeric(segs3$V17) >= segs3$ploid* 2] <- "gain"

segs3 <- dcast(V4 ~ sample, data = segs3, value.var = c("cnv"))
rownames(segs3) <- segs3[, 1]; segs3 <- segs3[, -1]

segs3 <- segs3[, colnames(segs3) %in% colnames(rbins)]
segs3 <- segs3[match(rownames(rbins), rownames(segs3)),
               match(colnames(rbins), colnames(segs3))]

# ------------------------------------------------------------------------------
# determine functional bins
# ------------------------------------------------------------------------------
if(file.exists("/projects/pangen/analysis/jtopham/data/fbin_res_object.RData")){
  load("/projects/pangen/analysis/jtopham/data/fbin_res_object.RData")
}else{
  res <- NULL
  for(i in 1:nrow(segs3)){
    if(sum(unlist(segs3[i, ])=="gain") < 3){this.g <- NA; this.gp <- NA
    }else{
      this.g <- median(unlist(rbins[i, ])[unlist(segs3[i, ]) == "gain"]) -
        median(unlist(rbins[i, ])[unlist(segs3[i, ]) == "neutral"])
      this.gp <- wilcox.test(unlist(rbins[i, ])[unlist(segs3[i, ]) == "gain"],
                             unlist(rbins[i, ])[unlist(segs3[i, ]) == "neutral"])$p.value
    }
    if(sum(unlist(segs3[i, ])=="hom loss") < 3){this.l <- NA; this.lp <- NA
    }else{
      this.l <- median(unlist(rbins[i, ])[unlist(segs3[i, ]) == "hom loss"]) -
        median(unlist(rbins[i, ])[unlist(segs3[i, ]) == "neutral"])
      this.lp <- wilcox.test(unlist(rbins[i, ])[unlist(segs3[i, ]) == "hom loss"],
                             unlist(rbins[i, ])[unlist(segs3[i, ]) == "neutral"])$p.value
    }
    res <- rbind(res, cbind(seg = rownames(segs3)[i],
                            gain.diff = this.g,
                            gain.pval = this.gp,
                            loss.diff = this.l,
                            loss.pval = this.lp))
  }
  res <- data.frame(res, stringsAsFactors = F)
  for(i in 2:ncol(res)){res[, i] <- as.numeric(res[, i])}
  res$gain.padj <- p.adjust(res$gain.pval, method = "BH")
  res$loss.padj <- p.adjust(res$loss.pval, method = "BH")
  
  res$gbin <- res$gain.diff > 0 & res$gain.padj < 0.05
  res$lbin <- res$loss.diff < 0 & res$loss.padj < 0.05
  
  res$gbin2 <- res$gbin; res$gbin2[is.na(res$gbin)] <- FALSE
  res$lbin2 <- res$lbin; res$lbin2[is.na(res$lbin)] <- FALSE
  
  save(res, file = "/projects/pangen/analysis/jtopham/data/fbin_res_object.RData")
}

tmp <- res[!is.na(res$gain.diff) | !is.na(res$loss.diff), c(1, 2, 4, 6:7, 10:11)]
tmp$gain.neg10 <- -log10(tmp$gain.padj)
tmp$loss.neg10 <- -log10(tmp$loss.padj)
tmp$bin <- "ns"
tmp$bin[!is.na(tmp$gain.diff)] <- "putative gain"
tmp$bin[!is.na(tmp$loss.diff)] <- "putative loss"
tmp$bin[!is.na(tmp$gain.diff) & !is.na(tmp$loss.diff)] <- "putative both"
tmp$bin[tmp$lbin2 == TRUE] <- "loss"
tmp$bin[tmp$gbin2 == TRUE] <- "gain"
tmp$neg10 <- tmp$gain.neg10
tmp$neg10[is.na(tmp$neg10)] <- tmp$loss.neg10[is.na(tmp$neg10)]
tmp$diff <- tmp$gain.diff
tmp$diff[is.na(tmp$diff)] <- tmp$loss.diff[is.na(tmp$diff)]

tmp$anno <- "other"
tmp$anno[tmp$seg == "7:98640001_98660001" & tmp$bin == "gain"] <- "SMURF1"
tmp$anno[tmp$seg == "1:200020001_200040001" & tmp$bin == "gain"] <- "NR5A2"
tmp$anno[tmp$seg == "1:214180001_214200001" & tmp$bin == "gain"] <- "PROX1"
tmp$anno[tmp$seg == "12:15780001_15800001" & tmp$bin == "gain"] <- "EPS8"
tmp$anno[tmp$seg == "19:45280001_45300001" & tmp$bin == "gain"] <- "CBLC"
tmp$anno[tmp$seg == "9:140120001_140140001" & tmp$bin == "loss"] <- "SLC34A3"
tmp$anno[tmp$seg == "9:21960001_21980001" & tmp$bin == "loss"] <- "CDKN2A"
tmp$anno[tmp$seg == "9:21860001_21880001" & tmp$bin == "loss"] <- "MTAP"
tmp$anno[tmp$seg == "9:21320001_21340001" & tmp$bin == "loss"]<- "KLHL9"

# get which genes are affected by gain/loss bins
gene <- read.delim("/projects/pangen/analysis/jtopham/data/allbins_vs_genes.tsv",
                   header = F, sep = '\t', stringsAsFactors = F)
gene$seg <- paste0(gene$V1, ":", gene$V2, "_", gene$V3)
gene$gbin <- res$gbin2[match(gene$seg, res$seg)]
gene$lbin <- res$lbin2[match(gene$seg, res$seg)]

tmp$gene <- gene$V14[match(tmp$seg, gene$seg)]
res2 <- tmp

# 7 x 5
p <- ggplot(tmp, aes(x = neg10, y = diff, color = bin))
p + geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("#bc4151", "#313695", "#532a66", "#e4b3b9", "#abcde3")) +
  geom_vline(xintercept = 1.3, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(data = tmp[tmp$anno != "other", ],
            aes(x = neg10, y = diff, label = anno), color = "black") +
  theme(axis.text = element_text(size = 8, color = "black"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(linetype = "dashed", color = "gray89"))

# Suppl table 3: Genomic bins (20kb) with significant (p<0.05) difference in
# expression (hom. loss or gain versus neutral) and absolute difference in
# median expression > 0.25 (log10 RPKM); n=274
st3 <- tmp[tmp$neg10>1.3 & abs(tmp$diff)>0.25,]
st3$gene <- gene$V14[match(st3$seg, gene$seg)]
st3 <- st3[order(st3$neg10, decreasing = T), ]
st3$chr <- gsub(":.*", "", st3$seg)
st3$start <- gsub("_.*", "", gsub(".*:", "", st3$seg))
st3$end <- gsub(".*_", "", st3$seg)
st3 <- st3[, c(15:17, 14, 10, 12, 11)]
st3 <- st3[!grepl("putative", st3$bin), ]
colnames(st3)[ncol(st3)] <- "padj"
st3$padj <- 10^(st3$padj*-1)
# write.table(st3, sep = '\t', quote = F, row.names = F,
#             file = "/projects/sftp/jtopham/incoming/fxbin_paper_supptable3.tsv")

goi <- unique(tmp$anno)[unique(tmp$anno)!="other"]

# donuts+
th$fxn <- res2$bin[match(th$seg, res2$seg)]
th$fxn[grepl("put", th$fxn)] <- "not fxn"
th$fxn[is.na(th$fxn)] <- "not fxn"
th$fxn2 <- paste0(th$Var1, "_", th$fxn)

donut <- data.frame(table(th$Var1, th$fxn), stringsAsFactors = F)
donut <- lapply(split(donut, f = donut$Var1), FUN = function(x){
  x$Freq2 <- (x$Freq / sum(x$Freq)) * 100
  return(x)
}); donut <- do.call(rbind, donut)

# overall 2444/28847 (8.47%) fxgain bins, 269/3252 (8.27%) fxloss bins
donut$Var2 <- as.character(donut$Var2)
donut$Var2[donut$Var2 == "not fxn" & donut$Var1 == "gain"] <- "gain; not fxn"
donut$Var2[donut$Var2 == "not fxn" & donut$Var1 == "hom loss"] <- "loss; not fxn"
p <- ggplot(donut, aes(x = 1, y = Freq2, fill = Var2))
p + geom_bar(stat = "identity", position = "stack") +
  coord_polar("y", start = 0) + facet_wrap(~ Var1, nrow = 1) +
  scale_fill_manual(values = c("#963440", "#d07a85", "#003e86", "#668bb6"))

# ------------------------------------------------------------------------------
# hits
# ------------------------------------------------------------------------------
tmp2 <- segs3; tmp2$seg <- rownames(tmp2); tmp2 <- melt(tmp2, id.vars = c("seg"))
tmp2$variable <- as.character(tmp2$variable)

goi2 <- data.frame(rbind(cbind(gene = "smurf1", 
        pt = tmp2$variable[tmp2$seg == "7:98640001_98660001" & tmp2$value == "gain"]),
        cbind(gene = "nr5a2", 
        pt = tmp2$variable[tmp2$seg == "1:200020001_200040001" & tmp2$value == "gain"]),
        cbind(gene = "prox1", 
        pt = tmp2$variable[tmp2$seg == "1:214180001_214200001" & tmp2$value == "gain"]),
        cbind(gene = "eps8", 
        pt = tmp2$variable[tmp2$seg == "12:15780001_15800001" & tmp2$value == "gain"]),
        cbind(gene = "cblc", 
        pt = tmp2$variable[tmp2$seg == "19:45280001_45300001" & tmp2$value == "gain"]),
        cbind(gene = "slc34a3", 
        pt = tmp2$variable[tmp2$seg == "9:140120001_140140001" & grepl("loss", tmp2$value)]),
        cbind(gene = "mtap", 
        pt = tmp2$variable[tmp2$seg == "9:21860001_21880001" & grepl("loss", tmp2$value)]),
        cbind(gene = "cdkn2a", 
        pt = tmp2$variable[tmp2$seg == "9:21960001_21980001" & grepl("loss", tmp2$value)]),
        cbind(gene = "klhl9", 
        pt = tmp2$variable[tmp2$seg == "9:21320001_21340001" & grepl("loss", tmp2$value)])),
        stringsAsFactors = F)

# ------------------------------------------------------------------------------
# (PanGen) oncoprint
# ------------------------------------------------------------------------------
pgenes <- c("CBLC", "EPS8", "SMURF1", "PROX1", "NR5A2", "KLHL9", "MTAP", 
            "CDKN2A", "SLC34A3", "KRAS", "TP53", "SMAD4", "ARID1A", 
            "MYC", "KMT2D", "KMT2C", "RNF43")
# cnv
cnv2 <- cnv[rownames(cnv) %in% pgenes, ]

cnv2$gene <- rownames(cnv2)
cnv2 <- melt(cnv2, id.vars = c("gene"))
cnv2$value <- gsub("[.]", "NA", cnv2$value)
cnv2$cnv <- "neutral"
cnv2$cnv[grep("NA", cnv2$value)] <- NA
cnv2$cnv[gsub(".*_", "", cnv2$value) == "0"] <- "het loss"
cnv2$cnv[gsub(".*_", "", cnv2$value) == "0" &
           gsub(".*-", "", gsub("_.*", "", cnv2$value)) == "0"] <- "hom loss"
# gains defined as being at least twice the ploidy (ie POG == 4)
cnv2$cnv[as.numeric(gsub(".*-", "", gsub("_.*", "", cnv2$value))) >= 
           as.numeric(gsub("-.*", "", cnv2$value)) * 2] <- "gain"
cnv2$variable <- as.character(cnv2$variable)

# row order
ro <- rev(c("CBLC", "EPS8", "SMURF1", "PROX1", "NR5A2", "KLHL9", "MTAP", 
            "CDKN2A", "SLC34A3", "KRAS", "TP53", "SMAD4", "ARID1A", 
            "MYC", "KMT2D", "KMT2C", "RNF43", "PABPC3"))
cnv2$gene <- factor(cnv2$gene, levels = ro)

# col order
tmp <- dcast(gene ~ variable, data = cnv2, value.var = c("cnv"))
rownames(tmp) <- tmp$gene; tmp <- tmp[, -1]
tmp[is.na(tmp)] <- "0"
for(i in 1:ncol(tmp)){
  for(j in 1:nrow(tmp)){
    if(rownames(tmp)[j] %in% c("SMURF1", "EPS8", "CBLC") &
       tmp[j, i] != "gain"){tmp[j, i] <- 0}
  }
}
tmp[tmp == "neutral"] <- "0"; tmp[tmp == "het loss"] <- "0"
tmp[tmp != "0"] <- "1"
for(i in 1:ncol(tmp)){tmp[, i] <- as.numeric(tmp[, i])}
scores <- apply(na.omit(tmp[rev(match(ro, rownames(tmp))), ]), 2, scoreCol)

cnv2$score <- unname(scores)[match(cnv2$variable, names(scores))]
cnv2 <- cnv2[order(cnv2$score, decreasing = T), ]
cnv2$variable <- factor(cnv2$variable, levels = unique(cnv2$variable))

p <- ggplot(cnv2, aes(x = variable, y = gene, fill = cnv))
p + geom_tile(color = "white", size = .25) +
  scale_fill_manual(values = c("#b20000", "#b2c5da", "#003e86", "gray92")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 14, color = "black"),
        panel.background = element_rect(fill = NA))

snv2 <- snv[!is.na(snv$gene), ]
snv2 <- snv2[snv2$gene != "", ]
snv2 <- snv2[snv2$gene %in% pgenes, ]
snv2$tmp <- factor(snv2$varType2, levels = c("SI/PP", "frameshift_INDEL", "nonsense",
                                             "missense", "inframe_INDEL", "splice", 
                                             "other", "other_silent"))
snv2$str <- paste0(snv2$sample, "__", snv2$gene)
snv2 <- snv2[order(snv2$str, snv2$tmp), ]
snv2 <- snv2[!duplicated(snv2$str), ]

snv2 <- snv2[!snv2$varType2 %in% c("splice", "other", "other_silent"), ]

# add wiltype rows
tmp <- NULL
for(i in pgenes){
  for(j in unique(snv$sample)){
    tmp <- c(tmp, paste0(j, "__", i))
  }
}

tmp <- data.frame(tmp, stringsAsFactors = F)
tmp <- cbind(sample = gsub("__.*", "", tmp$tmp),
             gene = gsub(".*__", "", tmp$tmp),
             vaf = NA, varType2 = "wildtype", tmp)
colnames(tmp)[5] <- "str"
tmp <- tmp[!tmp$str %in% snv2$str, ]

snv2 <- rbind(snv2[, c(1, 7, 10, 11, 13)], tmp)

# row order
snv2$gene <- factor(snv2$gene, levels = levels(cnv2$gene))

# col order
snv2$sample <- factor(snv2$sample, levels = levels(cnv2$variable))

# 12 x 4
p <- ggplot(snv2, aes(x = sample, y = gene, fill = varType2))
p + geom_tile(color = "white", size = .25) +
  scale_fill_manual(values = c("#d57d00", "#7a58b2", "#626262", 
                               "#8bc24c", "#626262", "gray92")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 14, color = "black"),
        panel.background = element_rect(fill = NA))

# fusions
fus2 <- fus[fus$gene2 %in% c("NRG1", "FGFR2", "BRAF", "NTRK2") | 
              fus$X.gene1 %in% c("NRG1", "FGFR2", "BRAF", "NTRK2"), c(1:3, 13:18, 23, 26:27)]
fus2 <- fus2[fus2$confidence == "high" & fus2$reading_frame == "in-frame", ]
fus2$gene <- "NRG1"
fus2$gene[fus2$gene2 == "BRAF" | fus2$X.gene1 == "BRAF"] <- "BRAF"
fus2$gene[fus2$gene2 == "FGFR2" | fus2$X.gene1 == "FGFR2"] <- "FGFR2"
fus2$gene[fus2$gene2 == "NTRK2" | fus2$X.gene1 == "NTRK2"] <- "NTRK2"
fus3 <- fus2[, c(1, ncol(fus2))]

# add WT status
fus3$str <- paste0(fus3$sample, "__", fus3$gene)
fus3$anno <- "fusion"
tmp <- NULL
for(i in c("NRG1", "FGFR2", "BRAF", "NTRK2")){
  for(j in intersect(fus$sample, snv2$sample)){
    tmp <- c(tmp, paste0(j, "__", i))
  }
}

tmp <- data.frame(tmp, stringsAsFactors = F)
tmp <- cbind(sample = gsub("__.*", "", tmp$tmp),
             gene = gsub(".*__", "", tmp$tmp), tmp,
             anno = "wildtype")
colnames(tmp)[3] <- "str"
tmp <- tmp[!tmp$str %in% fus3$str, ]
fus3 <- rbind(fus3, tmp)

fus3$sample <- factor(fus3$sample, levels = levels(snv2$sample))
fus3$gene <- factor(fus3$gene, levels = c("BRAF", "NTRK2", "FGFR2", "NRG1"))

p <- ggplot(fus3, aes(x = sample, y = gene, fill = anno))
p + geom_tile(color = "white", size = .25) +
  scale_fill_manual(values = c("#7b81aa", "#dddeff")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 14, color = "black"),
        panel.background = element_rect(fill = NA))

# upper bars (TMB)
pur2 <- pur[pur$sample %in% snv2$sample, ]
pur2$sample <- factor(pur2$sample, levels = levels(snv2$sample))

pur2$tmb <- tmb$tmb_total[match(pur2$sample, tmb$sample)]
pur2$hrd <- hrd$hr_status[match(pur2$sample, hrd$sample)]

p <- ggplot(pur2, aes(x = sample, y = tmb, fill = hrd))
p + geom_bar(stat = "identity", alpha = 0.9) +
  scale_fill_manual(values = c("#9c522f", "#d3b574")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(linetype = "dashed", color = "gray89"))

# lower subtype heatmap
pur3 <- melt(pur2[, c(1, 6, 3, 13, 14)], id.vars = c("sample"))
pur3$variable <- factor(pur3$variable, levels = rev(c("purist_call", "collisson",
                                                      "bailey", "metab")))

p <- ggplot(pur3, aes(x = sample, y = variable, fill = value))
p + geom_tile(color = "white", size = .25) +
  scale_fill_manual(values = c("#944c5e", "#463b66", "#99af5d",
                               "#a7b1e6", "#944c5e", "#62a1a9", "gray55",
                               "#8c492a", "#a7b1e6", "#463b66", "#d3b574", "#463b66")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 14, color = "black"),
        panel.background = element_rect(fill = NA))

clin2 <- clin[clin$sample %in% colnames(all.rna), ]
tmp <- clin2
tmp$gbinpt <- tmp$sample %in% goi2$pt[goi2$gene %in% c("cblc", "eps8", "smurf1")]
tmp$myc <- tmp$sample %in% cnv2$variable[cnv2$gene == "MYC" & cnv2$cnv == "gain"]
tmp$kras <- tmp$sample %in% cnv2$variable[cnv2$gene == "KRAS" & cnv2$cnv == "gain"]
tmp$cdkn2a <- tmp$sample %in% cnv2$variable[cnv2$gene == "CDKN2A" & cnv2$cnv == "hom loss"]

# MYC: 0.06155
fisher.test(table(tmp$gbinpt, tmp$myc))
# KRAS: 0.03154
# fisher.test(table(tmp$gbinpt, tmp$kras))
# # CDKN2A: 0.0577
# fisher.test(table(tmp$gbinpt, tmp$cdkn2a))

tmp$moff <- pur$purist_call[match(tmp$sample, pur$sample)]
fisher.test(table(tmp$gbinpt, tmp$moff)) # 0.2486
tmp$metab <- pur$metab[match(tmp$sample, pur$sample)]
fisher.test(table(tmp$gbinpt, tmp$metab)) # 0.2168
tmp$bail <- pur$bailey[match(tmp$sample, pur$sample)]
fisher.test(table(tmp$gbinpt, tmp$bail)) # 0.6487
tmp$coll <- pur$collisson[match(tmp$sample, pur$sample)]
fisher.test(table(tmp$gbinpt, tmp$coll)) # 0.8337

# vs TMB / HRD
tmp$tmb <- tmb$tmb_total[match(tmp$sample, tmb$sample)]
wilcox.test(tmp$tmb[tmp$gbinpt == TRUE], tmp$tmb[tmp$gbinpt == FALSE])$p.value # 1.0
tmp$hrd <- hrd$hr_status[match(tmp$sample, hrd$sample)]
fisher.test(table(tmp$gbinpt, tmp$hrd)) # 1.0

# mutual exclusivity (0.7314)
tmp <- cnv2[!is.na(cnv2$cnv), ]
tmp$variable <- as.character(tmp$variable)

fisher.test(matrix(c(length(intersect(tmp$variable[tmp$cnv == "gain" & tmp$gene == "EPS8"],
                        tmp$variable[tmp$cnv == "gain" & tmp$gene == "CBLC"])),
       length(tmp$variable[tmp$cnv == "gain" & tmp$gene == "EPS8"]),
       length(tmp$variable[tmp$cnv == "gain" & tmp$gene == "CBLC"]),
       length(intersect(tmp$variable[tmp$cnv != "gain" & tmp$gene == "EPS8"],
                        tmp$variable[tmp$cnv != "gain" & tmp$gene == "CBLC"]))),
       ncol = 2, nrow = 2), alternative = "less")

# ------------------------------------------------------------------------------
# vs clin variates
# ------------------------------------------------------------------------------
clin2 <- clin[clin$sample %in% colnames(all.rna), ]

cres <- NULL
for(i in unique(goi2$gene)){
  clin2$group <- clin2$sample %in% goi2$pt[goi2$gene == i]
  this.p <- 1 - pchisq(survdiff(Surv(os, censor) ~ group, 
                                data = clin2)$chisq, 
                       length(survdiff(Surv(os, censor) ~ group, 
                                       data = clin2)$n) - 1)
  this.sex <- fisher.test(table(clin2$group, clin2$sex))$p.value
  this.diab <- fisher.test(table(clin2$group, clin2$diabetes))$p.value
  this.eopc <- fisher.test(table(clin2$group, clin2$age<=55))$p.value
  this.bx <- fisher.test(table(clin2$group, clin2$bx=="liver"))$p.value
  cres <- rbind(cres, cbind(goi = i, pval = this.p,
                            sex = this.sex, diab = this.diab,
                            eopc = this.eopc, bx = this.bx))
}
cres <- data.frame(cres, stringsAsFactors = F)
for(i in 2:ncol(cres)){cres[, i] <- as.numeric(cres[, i])}

tmp <- clin2
plots <- NULL
for(i in 1:length(c("cblc", "eps8", "smurf1"))){
  tmp$group <- tmp$sample %in% goi2$pt[goi2$gene == c("cblc", "eps8", "smurf1")[i]]
  plots[[i]] <- ggsurv(survfit(Surv(os, censor) ~ group, 
                               data = tmp),
                       surv.col = c("#7d9fcb", "gray75"),
                       main = c("cblc", "eps8", "smurf1")[i],
                       size.est = 1.5, cens.size = 3.5, back.white = TRUE,
                       order.legend = FALSE, xlab = "time (months)")
}
plot_grid(plots[[1]], plots[[2]], plots[[3]], ncol = 3)

# CBLC yes/no = 5,64 - 2,13 - 0,3 - 0,1
# EPS8 yes/no = 4,65 - 1,14 - 0,3 - 0,1
# SMURF1 yes/no = 4,65 - 0,15 - 0,3 - 0,1
# tmp$group <- tmp$sample %in% goi2$pt[goi2$gene == "smurf1"]
# table(tmp$group[tmp$os >= 0])

# CBLC: med 17.1 vs 12.9; HR=0.77 95%CI=[0.28-2.1], p=0.61
# EPS8: med 11.7 vs 12.9; HR=0.94 95%CI=[0.29-3.0], p=0.91
# SMURF1: med 12.2 vs 13.0; HR=1.7 95%CI=[0.62-4.9], p=0.29
tmp$group <- tmp$sample %in% goi2$pt[goi2$gene == "cblc"]
res.cox <- coxph(Surv(os, censor) ~ group, 
                 data = tmp)
summary(res.cox)
1 - pchisq(survdiff(Surv(os, censor) ~ group, 
                    data = tmp)$chisq,
           length(survdiff(Surv(os, censor) ~ group,
                           data = tmp)$n) - 1)

# ------------------------------------------------------------------------------
# visualize focality and amplitude
# ------------------------------------------------------------------------------
foc <- res[, c(1, 10, 11)]
foc$gene <- gene$V14[match(foc$seg, gene$seg)]

foc$chr <- gsub(":.*", "", foc$seg)
foc$start <- as.numeric(gsub(".*:", "", gsub("_.*", "", foc$seg)))
foc$end <- as.numeric(gsub(".*_", "", foc$seg))

foc2 <- foc[foc$gene %in% goi & (foc$gbin2 == TRUE | foc$lbin2 == TRUE), ]
foc2 <- foc2[!duplicated(foc2$gene), ]

foc3 <- NULL
for(i in foc2$gene){
  these.segs <- foc$seg[foc$chr == foc2$chr[foc2$gene == i] &
                        foc$start > (foc2$start[foc2$gene == i] - 2e6) &
                        foc$end < (foc2$start[foc2$gene == i] + 2e6)]
  foc3 <- rbind(foc3, cbind(gene = i,
                            segs[segs$V4 %in% these.segs, c(1:5, 18:19, 21:22)]))
}
foc3 <- data.frame(foc3, stringsAsFactors = F)
foc3 <- foc3[order(foc3$V1, foc3$V2), ]
foc3$V4 <- factor(foc3$V4, levels = unique(foc3$V4))
colnames(foc3) <- c("gene", "sample", "chr", "start", "end", "seg",
                    "tcn", "mcn", "ploid", "cnv")
foc3$tcn <- as.numeric(foc3$tcn)
foc3$mcn <- as.numeric(foc3$mcn)

foc3$cnv2 <- foc3$cnv
foc3$cnv2[foc3$tcn >= (foc3$ploid * 2)] <- "gain 2x"
foc3$cnv2[foc3$tcn >= (foc3$ploid * 3)] <- "gain 3x"
foc3$cnv2[foc3$tcn >= (foc3$ploid * 4)] <- "gain 4x"
foc3$cnv2[foc3$tcn >= (foc3$ploid * 6)] <- "gain 6x"
foc3$cnv2[foc3$tcn >= (foc3$ploid * 8)] <- "gain 8x"
foc3$cnv2[foc3$tcn >= (foc3$ploid * 10)] <- "gain 10x"

tmp <- foc3; tmp$score <- 0
tmp$score[tmp$cnv2 == "hom loss" & tmp$gene %in% c("KLHL9", "MTAP")] <- 5
tmp$score[tmp$cnv2 == "het loss" & tmp$gene %in% c("KLHL9", "MTAP")] <- 1
tmp$score[grepl("gain", tmp$cnv2) & tmp$gene %in% c("KLHL9", "MTAP")] <- -1

tmp$score[tmp$cnv2 == "gain 2x" & !tmp$gene %in% c("KLHL9", "MTAP")] <- 2
tmp$score[tmp$cnv2 == "gain 3x" & !tmp$gene %in% c("KLHL9", "MTAP")] <- 5
tmp$score[tmp$cnv2 == "gain 4x" & !tmp$gene %in% c("KLHL9", "MTAP")] <- 8
tmp$score[tmp$cnv2 == "gain 6x" & !tmp$gene %in% c("KLHL9", "MTAP")] <- 12
tmp$score[tmp$cnv2 == "gain 8x" & !tmp$gene %in% c("KLHL9", "MTAP")] <- 18
tmp$score[tmp$cnv2 == "gain 10x" & !tmp$gene %in% c("KLHL9", "MTAP")] <- 25
tmp$score[tmp$cnv2 == "het loss" & !tmp$gene %in% c("KLHL9", "MTAP")] <- -1

tmp <- aggregate(score ~ sample + gene, data = tmp, sum)
tmp$str <- paste0(tmp$sample, "_", tmp$gene)
tmp <- tmp[order(tmp$score, decreasing = T), ]
tmp$score[tmp$str %in% paste0(goi2$pt, "_", toupper(goi2$gene))] <- 1e6

foc3$str <- paste0(foc3$sample, "_", foc3$gene)
foc3$str <- factor(foc3$str, levels = rev(tmp$str))
foc3$cnv2 <- factor(foc3$cnv2, levels = c("neutral", "het loss", "gain 2x",
                                          "gain 8x", "gain 4x", "gain 10x",
                                          "gain 6x", "hom loss", "gain 3x"))

# 20 x 6
p <- ggplot(foc3[foc3$gene != "MTAP", ], aes(x = seg, y = str, fill = cnv2))
p + geom_tile() + facet_wrap(~ gene, scales = "free", nrow = 1) +
  scale_fill_manual(values = c("gray92", "#b2c5da", "#e0b7b7", "#5b0909",
                               "#990f0f", "#1e0303", "#7a0c0c", "#003e86",
                               "#ad3e3e")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 7, color = "black"),
        panel.background = element_rect(fill = NA))

# upper tracks
foc3$gp <- -log10(res$gain.padj[match(foc3$seg, res$seg)])
foc3$lp <- -log10(res$loss.padj[match(foc3$seg, res$seg)])

foc3$binp <- foc3$gp
foc3$binp[foc3$gene %in% c("KLHL9", "MTAP")] <- foc3$lp[foc3$gene %in% c("KLHL9", "MTAP")]
foc3$bin[foc3$gene %in% c("KLHL9", "MTAP")] <- "put loss"
foc3$bin[!foc3$gene %in% c("KLHL9", "MTAP")] <- "put gain"
foc3$bin[foc3$seg %in% res2$seg[res2$bin == "gain"] &
         !foc3$gene %in% c("KLHL9", "MTAP")  ] <- "gain bin"
foc3$bin[foc3$seg %in% res2$seg[res2$bin == "loss"] &
         foc3$gene %in% c("KLHL9", "MTAP")] <- "loss bin"

tmp <- foc3[!duplicated(foc3$sample), ]
tmp$sample <- NULL

# 20 x 4
p <- ggplot(foc3[foc3$gene != "MTAP", ], aes(x = seg, y = binp, fill = bin))
p + geom_bar(stat = "identity") + facet_wrap(~ gene, scales = "free", nrow = 1) +
  scale_fill_manual(values = c("#5b0909", "#2d5eae", "#e0b7b7", "#b2c5da")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 7, color = "black"),
        axis.text.y = element_text(size = 7, color = "black"),
        panel.background = element_rect(fill = NA))

# show gene location
tmp <- foc3[foc3$gene %in% c("CBLC", "EPS8", "SMURF1"), ]
tmp$gene2 <- gene$V14[match(tmp$seg, gene$seg)]
tmp$gene2[!tmp$gene2 %in% c("CBLC", "EPS8", "SMURF1")] <- "other"

p <- ggplot(tmp, aes(x = seg, y = 1, fill = gene2))
p + geom_tile() + facet_wrap(~ gene, scales = "free", nrow = 1) +
  scale_fill_manual(values = c("#b76e79", "#999999", "#ff8200", "#58228b")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 7, color = "black"),
        panel.background = element_rect(fill = NA))

# side tracks
tmp <- all.rna
# remove genes with zero variance in any batch
for(b in unique(batches$p2)){
  vars <- apply(tmp[, colnames(tmp) %in%
                      batches$id[batches$p2 == b]], 1, var)
  vars <- vars[!is.na(vars)]
  tmp <- tmp[rownames(tmp) %in% names(vars)[vars > 0], ]
}

# normalize RNAseq data, correcting for RNAseq batch
pan2 <- ComBat(dat = as.matrix(tmp[, match(batches$id, colnames(tmp))]), 
               batch = batches$p2, par.prior = TRUE)
pan2 <- data.frame(pan2, stringsAsFactors = F)

tmp <- pan2; tmp$gene <- rownames(tmp); tmp <- melt(tmp, id.vars = c("gene"))
tmp$str <- paste0(tmp$variable, "_", tmp$gene)

foc3$exp <- tmp$value[match(foc3$str, tmp$str)]
tmp <- foc3[!duplicated(foc3$str), ]

# 11 x 6
p <- ggplot(tmp[tmp$gene != "MTAP", ], aes(x = str, y = exp, fill = exp))
p + geom_bar(stat = "identity", size = 0.25, color = "white") + 
  facet_wrap(~ gene, scales = "free", nrow = 1) +
  coord_flip() + scale_fill_gradientn(colors = ygob_pal[100:256]) +
  theme(axis.text.x = element_text(size = 7, color = "black"),
        axis.text.y = element_blank(),
        panel.background = element_rect(fill = NA))

# ------------------------------------------------------------------------------
# CDKN2A deletion breadth
# ------------------------------------------------------------------------------
cdb <- foc3[foc3$gene=="KLHL9", ]
cdb$gene <- NULL
cdb$gene <- gene$V14[match(cdb$seg, gene$seg)]

tmp <- segs[segs$V4 %in% gene$seg[gene$V14 %in% 
                                  c("CDKN2A", "CDKN2B", "MTAP", "KLHL9")], ]
tmp$gene <- gene$V14[match(tmp$V4, gene$seg)]
tmp$gene[grepl("CDKN2A", tmp$gene)] <- "CDKN2A"
tmp$gene[grepl("CDKN2B", tmp$gene)] <- "CDKN2B"
tmp$sample <- as.character(tmp$sample)

clin2 <- clin
clin2$klh <- clin2$sample %in% tmp$sample[tmp$gene == "KLHL9" & tmp$cnv == "hom loss"]
clin2$mtap <- clin2$sample %in% tmp$sample[tmp$gene == "MTAP" & tmp$cnv == "hom loss"]
clin2$cdkn2a <- clin2$sample %in% tmp$sample[tmp$gene == "CDKN2A" & tmp$cnv == "hom loss"]
clin2$cdkn2b <- clin2$sample %in% tmp$sample[tmp$gene == "CDKN2B" & tmp$cnv == "hom loss"]
for(i in 1:nrow(clin2)){
  clin2$cdb[i] <- paste(clin2$klh[i], clin2$mtap[i], clin2$cdkn2a[i], clin2$cdkn2b[i], sep = "_")
}
clin2$cdb[clin2$cdb == "FALSE_FALSE_FALSE_FALSE"] <- "none"
clin2$cdb[clin2$cdb == "FALSE_TRUE_TRUE_TRUE"] <- "N2A+N2B+MTAP"
clin2$cdb[clin2$cdb == "TRUE_TRUE_TRUE_TRUE"] <- "N2A+N2B+MTAP+KLH"
clin2$cdb[!clin2$cdb %in% c("none", "N2A+N2B+MTAP", "N2A+N2B+MTAP+KLH")] <- "other"

# 8 x 5
ggsurv(survfit(Surv(os, censor) ~ cdb, 
               data = clin2[clin2$cdb %in% c("N2A+N2B+MTAP",
                                             "N2A+N2B+MTAP+KLH"), ]),
       surv.col = c("#6ea6b2", "#ffa366"),
       size.est = 1.5, cens.size = 3.5, back.white = TRUE,
       order.legend = FALSE, xlab = "time (months)")

# HR=0.10 95%CI=[0.027-0.380]
res.cox <- coxph(Surv(os, censor) ~ cdb, 
                 data = clin2[clin2$cdb %in% c("N2A+N2B+MTAP",
                                               "N2A+N2B+MTAP+KLH"), ])
summary(res.cox)

# p=6.5e-5
1 - pchisq(survdiff(Surv(os, censor) ~ cdb, 
                    data = clin2[clin2$cdb %in% c("N2A+N2B+MTAP",
                                                  "N2A+N2B+MTAP+KLH"), ])$chisq,
           length(survdiff(Surv(os, censor) ~ cdb,
                           data = clin2[clin2$cdb %in% c("N2A+N2B+MTAP",
                                                         "N2A+N2B+MTAP+KLH"), ])$n) - 1)
cdb$group <- clin2$cdb[match(cdb$sample, clin2$sample)]
cdb$group <- factor(cdb$group, levels = c("none", "other", "N2A+N2B+MTAP",
                                          "N2A+N2B+MTAP+KLH"))
# fix removal of prev factors
cdb <- cdb[order(cdb$str), ]
cdb$str <- as.character(cdb$str)
cdb$str <- factor(cdb$str, levels = unique(cdb$str))
# reset factors again
cdb <- cdb[order(cdb$group, cdb$str), ]
cdb$str <- as.character(cdb$str)
cdb$str <- factor(cdb$str, levels = unique(cdb$str))

p <- ggplot(cdb, aes(x = seg, y = str, fill = cnv2))
p + geom_tile() +
  scale_fill_manual(values = c("gray92", "#b2c5da", "#003e86")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 7, color = "black"),
        panel.background = element_rect(fill = NA))

p <- ggplot(cdb, aes(x = seg, y = binp, fill = bin))
p + geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#2d5eae", "#b2c5da")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 7, color = "black"),
        axis.text.y = element_text(size = 7, color = "black"),
        panel.background = element_rect(fill = NA))

tmp <- cdb[!duplicated(cdb$str), ]

# 11 x 6
p <- ggplot(tmp, aes(x = str, y = exp, fill = exp))
p + geom_bar(stat = "identity", size = 0.25, color = "white") + 
  coord_flip() + scale_fill_gradientn(colors = ygob_pal[100:256]) +
  theme(axis.text.x = element_text(size = 7, color = "black"),
        axis.text.y = element_blank(),
        panel.background = element_rect(fill = NA))

# show where each gene lies
tmp <- cdb[!duplicated(cdb$seg), ]
tmp$gene[grepl("CDKN2A", tmp$gene)] <- "CDKN2A"
tmp$gene[grepl("CDKN2B", tmp$gene)] <- "CDKN2B"
tmp$gene[!tmp$gene %in% c("CDKN2A", "CDKN2B", "MTAP", "KLHL9")] <- "other"

p <- ggplot(tmp, aes(x = seg, y = 1, fill = gene))
p + geom_tile() +
  scale_fill_manual(values = c("#b76e79", "#999999", "#ff8200", "#58228b", "#b5ffa4")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 7, color = "black"),
        panel.background = element_rect(fill = NA))

pur$cbd <- clin2$cdb[match(pur$sample, clin2$sample)]

# ------------------------------------------------------------------------------
# external cohort validate
# ------------------------------------------------------------------------------
# TC
tmp <- tc
load("/projects/pangen/analysis/jtopham/data/all_pdac_tcs_nov2021.RData")
v.tc <- tc; tc <- tmp

# Hartwig
hpdac <- read.delim("/projects/hartwig/jtopham_prj/PDAC/samples.tsv",
                    sep = '\t', header = F, stringsAsFactors = F)$V1
load("/projects/hartwig/jtopham_prj/r_objects/cnv.RData")
rownames(hcnv) <- hcnv$gene
hcnv <- hcnv[, -c(1:4)]
hcnv <- hcnv[, colnames(hcnv) %in% hpdac]
load("/projects/hartwig/jtopham_prj/r_objects/norm_rna.RData")
rownames(hrna) <- hrna$geneid
hrna <- hrna[, -c(1:6)]
hrna <- hrna[, colnames(hrna) %in% hpdac]
load("/projects/hartwig/jtopham_prj/r_objects/snv.RData")
hsnv <- hsnv[hsnv$sample %in% hpdac, ]

# COMPASS
tmp <- all.rna
load("/projects/pangen/analysis/jtopham/data/all_rna_oct2022.RData")
rownames(all.rna) <- all.rna$geneid
all.rna <- all.rna[, -c(1:6)]
all.rna <- all.rna[, grepl("COMP", colnames(all.rna))]
comp.rna <- all.rna; all.rna <- tmp

tmp <- cnv
load("/projects/pangen/analysis/jtopham/data/all_pdac_cnv_oct2022.RData")
rownames(cnv) <- cnv$gene
cnv <- cnv[, grepl("COMP", colnames(cnv))]
comp.cnv <- cnv; cnv <- tmp

tmp <- snv
load("/projects/pangen/analysis/jtopham/data/all_pdac_snv_oct2022.RData")
snv <- snv[grepl("COMP", snv$sample), ]
comp.snv <- snv; snv <- tmp

# combine
# rna
hrna <- hrna[rownames(hrna) %in% intersect(rownames(hrna), rownames(comp.rna)), ]
comp.rna <- comp.rna[rownames(comp.rna) %in% intersect(rownames(hrna), rownames(comp.rna)), ]
hrna <- hrna[match(rownames(comp.rna), rownames(hrna)), ]
v.rna <- cbind(hrna, comp.rna)

ext <- ext[rownames(ext) %in% intersect(rownames(ext), rownames(v.rna)), ]
v.rna <- v.rna[rownames(v.rna) %in% intersect(rownames(ext), rownames(v.rna)), ]
ext <- ext[match(rownames(v.rna), rownames(ext)), ]
v.rna <- cbind(v.rna, ext)

# cnv
hcnv <- hcnv[rownames(hcnv) %in% intersect(rownames(hcnv), rownames(comp.cnv)), ]
comp.cnv <- comp.cnv[rownames(comp.cnv) %in% intersect(rownames(hcnv), rownames(comp.cnv)), ]
hcnv <- hcnv[match(rownames(comp.cnv), rownames(hcnv)), ]
v.cnv <- cbind(hcnv, comp.cnv)

rcnv <- rcnv[rownames(rcnv) %in% intersect(rownames(rcnv), rownames(v.cnv)), ]
v.cnv <- v.cnv[rownames(v.cnv) %in% intersect(rownames(rcnv), rownames(v.cnv)), ]
rcnv <- rcnv[match(rownames(v.cnv), rownames(rcnv)), ]
v.cnv <- cbind(v.cnv, rcnv)

v.cnv <- v.cnv[, colnames(v.cnv) %in% colnames(v.rna)]

# snv
v.snv <- rbind(hsnv, rsnv[, c(1:8, 11)], comp.snv[, c(1:8, 11)])

v.cnv2 <- v.cnv[rownames(v.cnv) %in% c(goi, "CDKN2B"), ]
v.cnv2$gene <- rownames(v.cnv2)
for(i in 1:ncol(v.cnv2)){v.cnv2[, i] <- as.character(v.cnv2[, i])}
v.cnv2 <- melt(v.cnv2, id.vars = c("gene"))
v.cnv2$value <- gsub("[.]", "NA", v.cnv2$value)
v.cnv2$cnv <- "neutral"
v.cnv2$cnv[grep("NA", v.cnv2$value)] <- NA
v.cnv2$cnv[gsub(".*_", "", v.cnv2$value) == "0"] <- "het loss"
v.cnv2$cnv[gsub(".*_", "", v.cnv2$value) == "0" &
           gsub(".*-", "", gsub("_.*", "", v.cnv2$value)) == "0"] <- "hom loss"
# gains defined as being at least twice the ploidy (ie POG == 4)
v.cnv2$variable <- as.character(v.cnv2$variable)
tmp <- v.cnv2[!v.cnv2$variable %in% dmap$icgc_donor_id[dmap$project_code=="PAAD-US"], ]
tmp$cnv[as.numeric(gsub(".*-", "", gsub("_.*", "", tmp$value))) >= 
           as.numeric(gsub("-.*", "", tmp$value)) * 2] <- "gain"
# for TCGA only, require gains to be greater than twice ploidy
tmp2 <- v.cnv2[v.cnv2$variable %in% dmap$icgc_donor_id[dmap$project_code=="PAAD-US"], ]
tmp2$cnv[as.numeric(gsub(".*-", "", gsub("_.*", "", tmp2$value))) > 
          as.numeric(gsub("-.*", "", tmp2$value)) * 2] <- "gain"
v.cnv2 <- rbind(tmp, tmp2)

tmp <- v.rna; tmp$gene <- rownames(tmp); tmp <- melt(tmp, id.vars = c("gene"))
tmp$str <- paste0(tmp$variable, "_", tmp$gene)
v.cnv2$str <- paste0(v.cnv2$variable, "_", v.cnv2$gene)
v.cnv2$exp <- tmp$value[match(v.cnv2$str, tmp$str)]
v.cnv2$cohort <- "ICGC"
v.cnv2$cohort[grepl("COMP", v.cnv2$variable)] <- "COMPASS"
v.cnv2$cohort[v.cnv2$variable %in% hpdac] <- "Hartwig"
v.cnv2$cohort[v.cnv2$variable %in% dmap$icgc_donor_id[dmap$project_code == "PAAD-US"]] <- "TCGA"

# ------------------------------------------------------------------------------
# validation with KLHL9
# ------------------------------------------------------------------------------
p <- ggplot(v.cnv2[v.cnv2$gene %in% c("KLHL9") &
                     !is.na(v.cnv2$cnv), ], aes(x = cnv, y = exp, fill = cnv))
p + geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1) + facet_wrap(gene ~ cohort, scales = "free") +
  scale_fill_manual(values = c("#ad3e3e", "#b2c5da", "#2d5aea", "gray92")) +
  theme(axis.text = element_text(size = 8, color = "black"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(linetype = "dashed", color = "gray89"))

# clinical
tmp <- clin
load("/projects/pangen/analysis/jtopham/data/all_clin_may2023.RData")
v.clin <- clin; clin <- tmp
v.clin <- v.clin[v.clin$sample %in% gsub("-[A-Z]", "", v.cnv2$variable), ]
v.clin <- v.clin[v.clin$sample %in% gsub("-[A-Z]", "", v.snv$sample), ]

v.clin$klh <- v.clin$sample %in% v.cnv2$variable[v.cnv2$gene == "KLHL9" & v.cnv2$cnv == "hom loss"]
v.clin$mtap <- v.clin$sample %in% v.cnv2$variable[v.cnv2$gene == "MTAP" & v.cnv2$cnv == "hom loss"]
v.clin$cdkn2a <- v.clin$sample %in% v.cnv2$variable[v.cnv2$gene == "CDKN2A" & v.cnv2$cnv == "hom loss"]
v.clin$cdkn2b <- v.clin$sample %in% v.cnv2$variable[v.cnv2$gene == "CDKN2B" & v.cnv2$cnv == "hom loss"]
for(i in 1:nrow(v.clin)){
  v.clin$cdb[i] <- paste(v.clin$klh[i], v.clin$mtap[i], v.clin$cdkn2a[i], v.clin$cdkn2b[i], sep = "_")
}
v.clin$cdb[v.clin$cdb == "FALSE_FALSE_FALSE_FALSE"] <- "none"
v.clin$cdb[v.clin$cdb == "FALSE_TRUE_TRUE_TRUE"] <- "N2A+N2B+MTAP"
v.clin$cdb[v.clin$cdb == "TRUE_TRUE_TRUE_TRUE"] <- "N2A+N2B+MTAP+KLH"
v.clin$cdb[!v.clin$cdb %in% c("none", "N2A+N2B+MTAP", "N2A+N2B+MTAP+KLH")] <- "other"

# 14 x 9
plots <- NULL
for(i in 1:length(unique(v.clin$cohort))){
  plots[[i]] <- ggsurv(survfit(Surv(os, censor) ~ cdb, 
                             data = v.clin[v.clin$cdb %in% c("N2A+N2B+MTAP",
                                                             "N2A+N2B+MTAP+KLH") &
                                           v.clin$cohort == unique(v.clin$cohort)[i], ]),
                     surv.col = c("#6ea6b2", "#ffa366"),
                     main = unique(v.clin$cohort)[i],
                     size.est = 1.5, cens.size = 3.5, back.white = TRUE,
                     order.legend = FALSE, xlab = "time (months)")
}
plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], ncol = 2)

# TCGA KLH/no: 17,6 - 8,4 - 3,3 - 2,2 - 1,0
# ICGC KLH/no: 11,20 - 1,4 - 0,1 - 0,0
# COMPASS KLH/no: 14,18 - 7,7 - 4,2 - 0,1
# Hartwig KLH/no: 7,4 - 6,3 - 3,1 - 3,1 - 0,1
#table(v.clin$cdb[v.clin$cohort == "ICGC" & v.clin$os >= 0])

# TCGA: 17 (13.6% - mOS=9.1) vs 6 (4.8% - mOS=15.6); HR=1.35, CI=0.33-5.6, p=0.68
# ICGC: 11 (12.0% - mOS=20.0) vs 20 (21.7% - mOS=18.9); HR=1.1, CI=0.49-2.6, p=0.79
# COMPASS: 14 (7.2% - mOS=9.7) vs 18 (9.2% - mOS=8.4); HR=0.77, CI=0.37-1.6, p=0.49
# Hartwig: 7 (15.2% - mOS=7.8) vs 4 (8.7% - mOS=4.0); HR=1.0, CI=0.25-4.2, p=0.98
res.cox <- coxph(Surv(os, censor) ~ cdb,
                 data = v.clin[v.clin$cdb %in% c("N2A+N2B+MTAP",
                                                 "N2A+N2B+MTAP+KLH")&
                                 v.clin$cohort == "Hartwig", ])
summary(res.cox)

1 - pchisq(survdiff(Surv(os, censor) ~ cdb,
                    data = v.clin[v.clin$cdb %in% c("N2A+N2B+MTAP",
                                                    "N2A+N2B+MTAP+KLH")&
                                    v.clin$cohort == "Hartwig", ])$chisq,
           length(survdiff(Surv(os, censor) ~ cdb,
                           data = v.clin[v.clin$cdb %in% c("N2A+N2B+MTAP",
                                                           "N2A+N2B+MTAP+KLH")&
                                           v.clin$cohort == "Hartwig", ])$n) - 1)

# ------------------------------------------------------------------------------
# mPDAC external oncoprints
# ------------------------------------------------------------------------------
# CNV + SNV data:
# COMPASS Hartwig    ICGC    TCGA 
#   195      46      92     125 
v.cnv3 <- v.cnv[rownames(v.cnv) %in% 
                  c("CBLC", "EPS8", "SMURF1", "PROX1", "NR5A2", "KLHL9", "MTAP", 
                    "CDKN2A", "KRAS", "TP53", "SMAD4", "MYC"), ]
v.cnv3 <- v.cnv3[, colnames(v.cnv3) %in% intersect(v.snv$sample, colnames(v.cnv3))]

v.cnv3$gene <- rownames(v.cnv3)
v.cnv3 <- melt(v.cnv3, id.vars = c("gene"))
v.cnv3$value <- gsub("[.]", "NA", v.cnv3$value)
v.cnv3$cnv <- "neutral"
v.cnv3$cnv[grep("NA", v.cnv3$value)] <- NA
v.cnv3$cnv[gsub(".*_", "", v.cnv3$value) == "0"] <- "het loss"
v.cnv3$cnv[gsub(".*_", "", v.cnv3$value) == "0" &
           gsub(".*-", "", gsub("_.*", "", v.cnv3$value)) == "0"] <- "hom loss"
# gains defined as being at least twice the ploidy (ie POG == 4)
v.cnv3$variable <- as.character(v.cnv3$variable)
tmp <- v.cnv3[!v.cnv3$variable %in% dmap$icgc_donor_id[dmap$project_code=="PAAD-US"], ]
tmp$cnv[as.numeric(gsub(".*-", "", gsub("_.*", "", tmp$value))) >= 
          as.numeric(gsub("-.*", "", tmp$value)) * 2] <- "gain"
# for TCGA only, require gains to be greater than twice ploidy
tmp2 <- v.cnv3[v.cnv3$variable %in% dmap$icgc_donor_id[dmap$project_code=="PAAD-US"], ]
tmp2$cnv[as.numeric(gsub(".*-", "", gsub("_.*", "", tmp2$value))) > 
           as.numeric(gsub("-.*", "", tmp2$value)) * 2] <- "gain"
v.cnv3 <- rbind(tmp, tmp2)

# row order
ro <- rev(c("CBLC", "EPS8", "SMURF1", "PROX1", "NR5A2", "KLHL9", "MTAP", 
            "CDKN2A", "KRAS", "TP53", "SMAD4", "MYC"))
v.cnv3$gene <- factor(v.cnv3$gene, levels = ro)

# col order
tmp <- dcast(gene ~ variable, data = v.cnv3, value.var = c("cnv"))
rownames(tmp) <- tmp$gene; tmp <- tmp[, -1]
tmp[is.na(tmp)] <- "0"
for(i in 1:ncol(tmp)){
  for(j in 1:nrow(tmp)){
    if(rownames(tmp)[j] %in% c("SMURF1", "EPS8", "CBLC") &
       tmp[j, i] != "gain"){tmp[j, i] <- 0}
  }
}
tmp[tmp == "neutral"] <- "0"; tmp[tmp == "het loss"] <- "0"
tmp[tmp != "0"] <- "1"
for(i in 1:ncol(tmp)){tmp[, i] <- as.numeric(tmp[, i])}
scores <- apply(na.omit(tmp[rev(match(ro, rownames(tmp))), ]), 2, scoreCol)

v.cnv3$score <- unname(scores)[match(v.cnv3$variable, names(scores))]
v.cnv3 <- v.cnv3[order(v.cnv3$score, decreasing = T), ]
v.cnv3$variable <- factor(v.cnv3$variable, levels = unique(v.cnv3$variable))
v.cnv3$cohort <- v.cnv2$cohort[match(v.cnv3$variable, v.cnv2$variable)]

# 18 x 10
p <- ggplot(v.cnv3, aes(x = variable, y = gene, fill = cnv))
p + geom_tile(color = "white", size = .25) +
  facet_wrap(~ cohort, scales = "free") +
  scale_fill_manual(values = c("#b20000", "#b2c5da", "#003e86", "gray92")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 14, color = "black"),
        panel.background = element_rect(fill = NA))

v.snv2 <- v.snv[!is.na(v.snv$gene), ]
v.snv2 <- v.snv2[v.snv2$sample %in% intersect(v.snv2$sample, v.cnv3$variable), ]
pts <- intersect(v.snv2$sample, v.cnv3$variable)
v.snv2 <- v.snv2[v.snv2$gene != "", ]
v.snv2 <- v.snv2[v.snv2$gene %in% 
                   c("CBLC", "EPS8", "SMURF1", "PROX1", "NR5A2", "KLHL9", "MTAP", 
                     "CDKN2A", "KRAS", "TP53", "SMAD4", "MYC"), ]
v.snv2$tmp <- factor(v.snv2$varType2, levels = c("SI/PP", "frameshift_INDEL", "nonsense",
                                             "missense", "inframe_INDEL", "splice", 
                                             "other", "other_silent"))
v.snv2$str <- paste0(v.snv2$sample, "__", v.snv2$gene)
v.snv2 <- v.snv2[order(v.snv2$str, v.snv2$tmp), ]
v.snv2 <- v.snv2[!duplicated(v.snv2$str), ]

v.snv2 <- v.snv2[!v.snv2$varType2 %in% c("splice", "other", "other_silent"), ]

# add wiltype rows
tmp <- NULL
for(i in c("CBLC", "EPS8", "SMURF1", "PROX1", "NR5A2", "KLHL9", "MTAP", 
           "CDKN2A", "KRAS", "TP53", "SMAD4", "MYC")){
  for(j in pts){
    tmp <- c(tmp, paste0(j, "__", i))
  }
}

tmp <- data.frame(tmp, stringsAsFactors = F)
tmp <- cbind(sample = gsub("__.*", "", tmp$tmp),
             gene = gsub(".*__", "", tmp$tmp),
             varType2 = "wildtype", tmp)
colnames(tmp)[4] <- "str"
tmp <- tmp[!tmp$str %in% v.snv2$str, ]

v.snv2 <- rbind(v.snv2[, c(1, 7, 9, 11)], tmp)

# row order
v.snv2$gene <- factor(v.snv2$gene, levels = levels(v.cnv3$gene))

# col order
v.snv2$sample <- factor(v.snv2$sample, levels = levels(v.cnv3$variable))
v.snv2$cohort <- "ICGC"
v.snv2$cohort[grepl("COMP", v.snv2$sample)] <- "COMPASS"
v.snv2$cohort[v.snv2$sample %in% hpdac] <- "Hartwig"
v.snv2$cohort[v.snv2$sample %in% dmap$icgc_donor_id[dmap$project_code == "PAAD-US"]] <- "TCGA"

# COMPASS Hartwig    ICGC    TCGA 
#   195      46      92     125 
# 18 x 10
p <- ggplot(v.snv2, aes(x = sample, y = gene, fill = varType2))
p + geom_tile(color = "white", size = .25) +
  facet_wrap(~ cohort, scales = "free") +
  scale_fill_manual(values = c("#d57d00", "#7a58b2", "#626262", 
                               "#8bc24c", "#626262", "gray92")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 14, color = "black"),
        panel.background = element_rect(fill = NA))

# freqs
tmp <- v.cnv3; tmp$met <- tmp$cohort %in% c("Hartwig", "COMPASS")
tmp$gain2 <- tmp$cnv == "gain"
fisher.test(table(tmp$met[tmp$gene=="CBLC"], tmp$gain2[tmp$gene=="CBLC"])) # 0.2075
fisher.test(table(tmp$met[tmp$gene=="EPS8"], tmp$gain2[tmp$gene=="EPS8"])) # 0.03394
fisher.test(table(tmp$met[tmp$gene=="SMURF1"], tmp$gain2[tmp$gene=="SMURF1"])) # 0.1566
table(v.cnv3$cohort[v.cnv3$gene=="CBLC"], v.cnv3$cnv[v.cnv3$gene=="CBLC"])

# versus MYC amp
tmp <- v.cnv3[!duplicated(v.cnv3$variable), c(2,6)]
tmp$gbinpt <- tmp$variable %in% v.cnv3$variable[v.cnv3$gene %in% c("CBLC", "EPS8", "SMURF1") &
                                                  v.cnv3$cnv == "gain"]
tmp$myc <- tmp$variable %in% v.cnv3$variable[v.cnv3$gene %in% c("MYC") &
                                               v.cnv3$cnv == "gain"]

# TCGA (0.01768), ICGC (4.4e-4), COMPASS (0.6642), Hartwig (1.0)
fisher.test(table(tmp$gbinpt[tmp$cohort=="Hartwig"], tmp$myc[tmp$cohort=="Hartwig"]))

# mutual exclusivity
tmp <- v.cnv3[!duplicated(v.cnv3$variable), c(2,6)]
tmp$cblc <- "none"; tmp$eps8 <- "none"; tmp$smurf1 <- "none"
tmp$cblc[tmp$variable %in% v.cnv3$variable[v.cnv3$gene %in% c("CBLC") & v.cnv3$cnv == "gain"]] <- "cblc"
tmp$eps8[tmp$variable %in% v.cnv3$variable[v.cnv3$gene %in% c("EPS8") & v.cnv3$cnv == "gain"]] <- "eps8"
tmp$smurf1[tmp$variable %in% v.cnv3$variable[v.cnv3$gene %in% c("SMURF1") & v.cnv3$cnv == "gain"]] <- "smurf1"
tmp$goi <- NA
for(i in 1:nrow(tmp)){tmp$goi[i] <- paste0(tmp$cblc[i], "_", tmp$eps8[i], "_", tmp$smurf1[i])}
tmp$tmp <- grepl("none_none_none", tmp$goi)
table(tmp$cohort, grepl("none_eps8_none|none_none_s|c_none_none", tmp$goi))

# simplified CNV oncoprints for vali cohorts
# 18 x 10
tmp <- v.cnv3[v.cnv3$gene %in% c("CBLC", "EPS8", "SMURF1"), ]
tmp <- tmp[tmp$variable %in% tmp$variable[tmp$cnv == "gain" &
                                          tmp$gene %in% c("CBLC", "EPS8", "SMURF1")], ]
tmp$gene <- factor(as.character(tmp$gene), levels = rev(c("CBLC", "EPS8", "SMURF1")))
tmp$var2 <- factor(as.character(tmp$variable),
                   levels = unique(tmp$variable[order(tmp$cohort, tmp$variable)]))
p <- ggplot(tmp, aes(x = var2, y = gene, fill = cnv))
p + geom_tile(color = "white", size = .25) +
  scale_fill_manual(values = c("#b20000", "#b2c5da", "gray92")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 14, color = "black"),
        panel.background = element_rect(fill = NA))

# donuts
donut <- v.cnv3[!duplicated(v.cnv3$variable), c(2, 6)]
donut$goi <- donut$variable %in% v.cnv3$variable[v.cnv3$cnv == "gain" &
                                                 v.cnv3$gene %in% c("CBLC", "EPS8", "SMURF1")]
donut$myc <- donut$variable %in% v.cnv3$variable[v.cnv3$cnv == "gain" &
                                                 v.cnv3$gene == "MYC"]

donut <- rbind(donut, cbind(variable = unique(cnv2$variable),
                            cohort = "PanGen",
                            goi = unique(cnv2$variable) %in% 
                              goi2$pt[goi2$gene %in% c("cblc", "eps8", "smurf1")],
                            myc = unique(cnv2$variable) %in% cnv2$variable[cnv2$cnv == "gain" &
                                                                             cnv2$gene == "MYC"]))
donut$variable <- as.character(donut$variable)

donut <- data.frame(table(donut$cohort, donut$goi, donut$myc), stringsAsFactors = F)
colnames(donut)[2:3] <- c("goi", "myc")

donut <- lapply(split(donut, f = paste0(donut$Var1, donut$goi)), FUN = function(x){
  x$Freq2 <- (x$Freq / sum(x$Freq)) * 100
  return(x)
}); donut <- do.call(rbind, donut)

donut$goi <- as.character(donut$goi); donut$myc <- as.character(donut$myc)
donut$goi[donut$goi == TRUE] <- "goi"; donut$goi[donut$goi == "FALSE"] <- "wt"
# 9 x 4
p <- ggplot(donut[donut$goi == "goi", ], aes(x = 1, y = Freq2, fill = myc))
p + geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Var1) +
  coord_polar("y", start = 0) + facet_wrap(~ Var1, nrow = 1) +
  scale_fill_manual(values = c("#545454", "#c58d9a"))

p <- ggplot(donut[donut$goi == "wt", ], aes(x = 1, y = Freq2, fill = myc))
p + geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Var1) +
  coord_polar("y", start = 0) + facet_wrap(~ Var1, nrow = 1) +
  scale_fill_manual(values = c("#545454", "#c58d9a"))

# vali boxplots (fig 2 E)
# -- these are the ones that are used --
tmp <- v.cnv2
tmp$value <- NULL
tmp2 <- pan2; tmp2$gene <- rownames(tmp2); tmp2 <- melt(tmp2, id.vars = c("gene"))
tmp2$str <- paste0(tmp2$variable, "_", tmp2$gene)

tmp <- rbind(tmp, cbind(gene = as.character(cnv2$gene), variable = as.character(cnv2$variable),
                        cnv = cnv2$cnv, str = paste0(as.character(cnv2$variable),
                                                     "_", as.character(cnv2$gene)),
                        exp = tmp2$value[match(paste0(as.character(cnv2$variable),
                                                      "_", as.character(cnv2$gene)), tmp2$str)],
                        cohort = "PanGen"))
tmp$exp <- as.numeric(tmp$exp)
tmp <- tmp[tmp$gene %in% c("CBLC", "EPS8", "SMURF1") &
           !is.na(tmp$cnv), ]
tmp$cohort <- factor(tmp$cohort, levels = c("PanGen", "TCGA", "ICGC", "COMPASS", "Hartwig"))
tmp$cnv <- factor(tmp$cnv, levels = c("gain", "neutral", "het loss", "hom loss"))

# 9 x 8
p <- ggplot(tmp, aes(x = cnv, y = exp, fill = cnv))
p + geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1) + facet_wrap(cohort ~ gene, scales = "free", ncol = 3) +
  scale_fill_manual(values = c("#b20000", "gray92", "#b2c5da", "#003e86")) +
  theme(axis.text = element_text(size = 8, color = "black"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(linetype = "dashed", color = "gray89"))

# CBLC: 5.5e-4, 0.042, 0.008506, 0.0013, 0.0359
# EPS8: 0.0023, 0.11656, 0.7100 , 0.982, 0.408
# SMURF1: 8.8e-4, 0.001252, 0.1049, 3.24e-6,  0.00884
wilcox.test(tmp$exp[tmp$cohort == "TCGA" & tmp$cnv == "gain" &
                    tmp$gene == "EPS8"],
            tmp$exp[tmp$cohort == "TCGA" & tmp$cnv != "gain" &
                    tmp$gene == "EPS8"])$p.value

# RNA differences at a glance
v.rna2 <- v.rna; v.rna2$gene <- rownames(v.rna2); v.rna2 <- melt(v.rna2, id.vars = c("gene"))
v.rna2$cohort <- "ICGC"
v.rna2$cohort[grepl("COMP", v.rna2$variable)] <- "COMPASS"
v.rna2$cohort[v.rna2$variable %in% hpdac] <- "Hartwig"
v.rna2$cohort[v.rna2$variable %in% dmap$icgc_donor_id[dmap$project_code == "PAAD-US"]] <- "TCGA"

v.cnv4 <- v.cnv[rownames(v.cnv) %in% c("EGFR", "EPS8", "CBLC"), ]
v.cnv4 <- v.cnv4[, colnames(v.cnv4) %in% intersect(v.snv$sample, colnames(v.cnv4))]

v.cnv4$gene <- rownames(v.cnv4)
v.cnv4 <- melt(v.cnv4, id.vars = c("gene"))
v.cnv4$value <- gsub("[.]", "NA", v.cnv4$value)
v.cnv4$cnv <- "neutral"
v.cnv4$cnv[grep("NA", v.cnv4$value)] <- NA
v.cnv4$cnv[gsub(".*_", "", v.cnv4$value) == "0"] <- "het loss"
v.cnv4$cnv[gsub(".*_", "", v.cnv4$value) == "0" &
             gsub(".*-", "", gsub("_.*", "", v.cnv4$value)) == "0"] <- "hom loss"
# gains defined as being at least twice the ploidy (ie POG == 4)
v.cnv4$variable <- as.character(v.cnv4$variable)
tmp <- v.cnv4[!v.cnv4$variable %in% dmap$icgc_donor_id[dmap$project_code=="PAAD-US"], ]
tmp$cnv[as.numeric(gsub(".*-", "", gsub("_.*", "", tmp$value))) >= 
          as.numeric(gsub("-.*", "", tmp$value)) * 2] <- "gain"
# for TCGA only, require gains to be greater than twice ploidy
tmp2 <- v.cnv4[v.cnv4$variable %in% dmap$icgc_donor_id[dmap$project_code=="PAAD-US"], ]
tmp2$cnv[as.numeric(gsub(".*-", "", gsub("_.*", "", tmp2$value))) > 
           as.numeric(gsub("-.*", "", tmp2$value)) * 2] <- "gain"
v.cnv4 <- rbind(tmp, tmp2)

v.rna2 <- v.rna2[v.rna2$gene %in% c("EGFR", "EPS8", "CBLC"), ]
tmp <- dcast(variable + cohort ~ gene, data = v.rna2, value.var = c("value"))

v.cnv4$str <- paste0(v.cnv4$gene, "_", v.cnv4$variable)
tmp$str <- paste0("EPS8", "_", tmp$variable)
tmp$eps8_cnv <- v.cnv4$cnv[match(tmp$str, v.cnv4$str)]
tmp$str <- paste0("CBLC", "_", tmp$variable)
tmp$cblc_cnv <- v.cnv4$cnv[match(tmp$str, v.cnv4$str)]
tmp$str <- NULL

# COMPASS Hartwig    ICGC    TCGA 
#    195      46     149     129 

p <- ggplot(tmp, aes(x = CBLC, y = EGFR, color = cblc_cnv))
p + geom_point(size = 2) + stat_smooth(method = "lm", color = "black", se = F) +
  facet_wrap(~ cohort, scales = "free") +
  scale_color_manual(values = c("#b20000", "#b2c5da", "gray92", "black")) +
  theme(axis.text = element_text(size = 8, color = "black"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(linetype = "dashed", color = "gray89"))

p <- ggplot(tmp, aes(x = EPS8, y = EGFR, color = eps8_cnv))
p + geom_point(size = 2) + stat_smooth(method = "lm", color = "black", se = F) +
  facet_wrap(~ cohort, scales = "free") +
  scale_color_manual(values = c("#b20000", "#b2c5da", "gray92", "black")) +
  theme(axis.text = element_text(size = 8, color = "black"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(linetype = "dashed", color = "gray89"))

# ------------------------------------------------------------------------------
# clinical (all cohorts)
# ------------------------------------------------------------------------------
tmp <- clin
load("/projects/pangen/analysis/jtopham/data/all_clin_may2023.RData")
a.clin <- clin; clin <- tmp
a.clin <- a.clin[a.clin$sample %in% c(gsub("-[A-Z]", "", v.cnv2$variable),
                                      colnames(all.rna)), ]

# 13 x 15
plots <- NULL; k = 1
for(i in 1:length(unique(a.clin$cohort))){
  for(j in c("CBLC", "EPS8", "SMURF1")){
    a.clin$group <- NA; a.clin$group2 <- NA
    a.clin$group[a.clin$sample %in% cnv2$variable[cnv2$cnv == "gain" & cnv2$gene == j]] <- "gain"
    a.clin$group[a.clin$sample %in% cnv2$variable[cnv2$cnv != "gain" & cnv2$gene == j]] <- "no gain"
    a.clin$group2[a.clin$sample %in% v.cnv2$variable[v.cnv2$cnv == "gain" & v.cnv2$gene == j]] <- "gain"
    a.clin$group2[a.clin$sample %in% v.cnv2$variable[v.cnv2$cnv != "gain" & v.cnv2$gene == j]] <- "no gain"
    a.clin$group[is.na(a.clin$group)] <- a.clin$group2[is.na(a.clin$group)]
    this.p <- 1 - pchisq(survdiff(Surv(os, censor) ~ group, 
                                  data = a.clin[a.clin$cohort == unique(a.clin$cohort)[i], ])$chisq, 
                         length(survdiff(Surv(os, censor) ~ group, 
                                         data = a.clin[a.clin$cohort == unique(a.clin$cohort)[i], ])$n) - 1)
    plots[[k]] <- ggsurv(survfit(Surv(os, censor) ~ group, 
                                 data = a.clin[a.clin$cohort == unique(a.clin$cohort)[i], ]),
                         surv.col = c("#6ea6b2", "#ffa366"),
                         size.est = 1.5, cens.size = 3.5, back.white = TRUE,
                         order.legend = FALSE, xlab = "time (months)",
                         main = paste0(unique(a.clin$cohort)[i], "_", j, "_", 
                                       sum(a.clin$group[a.clin$cohort == unique(a.clin$cohort)[i]] == "gain",
                                           na.rm = T),
                                       "_vs_",
                                       sum(a.clin$group[a.clin$cohort == unique(a.clin$cohort)[i]] != "gain",
                                           na.rm = T), "_", signif(this.p, 3)))
    k <- k + 1
  }
}
plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]],
          plots[[6]], plots[[7]], plots[[8]], plots[[9]], plots[[10]],
          plots[[11]], plots[[12]], plots[[13]], plots[[14]], plots[[15]], ncol = 3)

# ------------------------------------------------------------------------------
# pan-cancer analysis
# ------------------------------------------------------------------------------
# TCGA
# meta data
meta.t <- read.delim(paste0("/projects/pangen/analysis/jtopham/data/tcga/",
                            "complete_TCGA_metadata.tsv"),
                     stringsAsFactors = F, sep = '\t', header = T)

# get all TCGA RNAseq (RSEM) data
load("/projects/pangen/analysis/jtopham/data/tcga/RNAseq_data_RSEM.RData")
tt <- tcga.all
colnames(tcga.all) <- gsub("[.].*", "", colnames(tcga.all))

# filter for cancer samples
tcga.all$sample_id <- as.character(tcga.all$sample_id)
tcga.all <- tcga.all[!do.call(rbind, strsplit(tcga.all$sample_id, "-"))[,4] %in% 
                       c("11A", "11B", "11C"), ]

# filter duplicated samples
tcga.all <- tcga.all[order(tcga.all$sample_id), ]
tmp <- paste(do.call(rbind, strsplit(tcga.all$sample_id, "-"))[, 1],
             do.call(rbind, strsplit(tcga.all$sample_id, "-"))[, 2],
             do.call(rbind, strsplit(tcga.all$sample_id, "-"))[, 3], 
             sep = "-")

tcga.all <- tcga.all[!duplicated(tmp), ]
tcga.all$cancer_type <- as.character(tcga.all$cancer_type)

# remove non-PAAD PAADs
tcga.all$donor <- paste(do.call(rbind, strsplit(tcga.all$sample_id, "-"))[, 1],
                        do.call(rbind, strsplit(tcga.all$sample_id, "-"))[, 2],
                        do.call(rbind, strsplit(tcga.all$sample_id, "-"))[, 3],
                        sep = "-")
tt <- tcga.all
tcga.all <- tcga.all[tcga.all$cancer_type != "PAAD" |
                     tcga.all$donor %in%
                     dmap$submitted_donor_id[dmap$icgc_donor_id %in%
                                             colnames(v.rna)], ]
tcga.all$donor <- NULL

# filter for CA types with at least 100 samples
filt <- table(tcga.all$cancer_type)

# removes "ACC"  "CHOL" "DLBC" "KICH" "MESO" "READ" "UCS"  "UVM"
tcga.all <- tcga.all[!tcga.all$cancer_type %in% names(filt)[filt < 100], ]

# normalization
tcga.all2 <- tcga.all[, colnames(tcga.all) %in% c("cancer_type", "sample_id",
                                                  "CBLC", "EPS8", "SMURF1")]

# 24 different cancer types, each with >= 120 samples
# log10 transform
for(i in 3:ncol(tcga.all2)){tcga.all2[, i] <- log10(tcga.all2[, i] + 1)}

if(file.exists("/projects/pangen/analysis/jtopham/data/tcga/pan_cnv_mapped_goi.RData")){
  load("/projects/pangen/analysis/jtopham/data/tcga/pan_cnv_mapped_goi.RData")
}else{
  files <- list.files("/projects/pangen/analysis/jtopham/data/tcga/pancan_cnv",
                     full.names = T)
  files <- grep("mapped", files, value = T)
  pancnv <- NULL
  for(i in 1:length(files)){
    if(file.info(files[i])$size == 0){next}
    tmp <- read.delim(files[i], stringsAsFactors = F, sep = '\t', header = F)
    # hefty filter as dframe is too large to handle
    tmp <- tmp[tmp$V5 %in% goi, ]
    pancnv <- rbind(pancnv, tmp[, c(9, 10, 5, 11, 12)])
  }
  colnames(pancnv) <- c("sample", "catype", "gene", "numprobes", "meanseg")
  save(pancnv, file = "/projects/pangen/analysis/jtopham/data/tcga/pan_cnv_mapped_goi.RData")
}

pancnv <- pancnv[pancnv$sample %in% tcga.all2$sample_id, ]
tcga.all2 <- tcga.all2[tcga.all2$sample_id %in% pancnv$sample, ] # 8163

# threshold of segmean = +/- 0.2 and # probes >= 10 (note diff gain thresh used later)
pancnv$numprobes <- as.numeric(pancnv$numprobes)
pancnv$meanseg <- as.numeric(pancnv$meanseg)
pancnv$cnv <- "neutral"
pancnv$cnv[pancnv$numprobes >= 10 & pancnv$meanseg < -0.2] <- "loss"
pancnv$cnv[pancnv$numprobes >= 10 & pancnv$meanseg > 0.2] <- "gain"
pancnv <- pancnv[pancnv$gene %in% c("EPS8", "CBLC", "SMURF1"), ]
pancnv$str <- paste0(pancnv$sample, "_", pancnv$gene)

tmp <- melt(tcga.all2, id.vars = c("cancer_type", "sample_id"))
tmp$str <- paste0(tmp$sample_id, "_", tmp$variable)
pancnv$exp <- tmp$value[match(pancnv$str, tmp$str)]
pancnv$paad <- pancnv$catype == "PAAD"

pancnv$catype2 <- "other"
pancnv$catype2[pancnv$catype %in% c("COAD", "READ", "ESCA", "STAD", "LIHC")] <- "GI"
pancnv$catype2[pancnv$catype %in% c("BLCA", "KIRC", "KIRP", "PRAD", "TGCT")] <- "urologic"
pancnv$catype2[pancnv$catype %in% c("CESC", "OV", "UCEC")] <- "gyna"
pancnv$catype2[pancnv$catype %in% c("LUAD", "LUSC")] <- "resp"
pancnv$catype2[pancnv$catype %in% c("GBM", "LGG")] <- "brain"
pancnv$catype2[pancnv$catype %in% c("HNSC", "PCPG", "THCA")] <- "head neck"
pancnv$catype2[pancnv$catype %in% c("BRCA")] <- "breast"
pancnv$catype2[pancnv$catype %in% c("PAAD")] <- "pancreas"
pair_pal <- pair_pal[-c(5,7,11)]
pair_pal[5] <- "gray12"
pair_pal[7] <- "#0a4a4a"
pancnv$paad <- factor(pancnv$paad, levels = c("TRUE", "FALSE"))

# supp figure to show selection of 0.49
p <- ggplot(pancnv, aes(x = meanseg))
p + geom_histogram() +
  geom_vline(xintercept = 0.49, linetype = "dashed") +
  theme(axis.text = element_text(size = 8, color = "black"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(linetype = "dashed", color = "gray89"))

# fig 4a. 7 x 8.5
p <- ggplot(pancnv, aes(x = meanseg, y = exp, color = catype2))
p + geom_point(size = 2, alpha = 0.8) + facet_wrap(paad ~ gene, ncol = 2) +
  geom_vline(xintercept = 0.49, linetype = "dashed") +
  scale_color_manual(values = pair_pal) +
  theme(axis.text = element_text(size = 8, color = "black"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(linetype = "dashed", color = "gray89"))

# fig 4b
pancnv$str <- paste0(pancnv$catype, "_", pancnv$gene)
cve <- NULL
for(i in unique(pancnv$str)){
  tmp <- pancnv[pancnv$str == i & pancnv$numprobes >= 10, ]
  cve <- rbind(cve, cbind(catype = unique(tmp$catype),
                          gene = unique(tmp$gene),
                          p = cor.test(tmp$meanseg, tmp$exp,
                                       method = "spearman")$p.value,
                          cor = cor(tmp$meanseg, tmp$exp, method = "spearman")))
}
cve <- data.frame(cve, stringsAsFactors = F)
for(i in 3:4){cve[, i] <- as.numeric(cve[, i])}
cve$padj <- NA
for(i in unique(cve$gene)){
  cve$padj[cve$gene == i] <- p.adjust(cve$p[cve$gene == i], method = "BH")
}
cve$neg10 <- -log10(cve$padj)
cve$gene <- factor(cve$gene, levels = rev(c("CBLC", "EPS8", "SMURF1")))
cve <- cve[order(cve$gene, cve$cor, decreasing = T), ]
cve$catype <- factor(cve$catype, levels = rev(unique(cve$catype)))
cve$gene <- factor(cve$gene, levels = c("CBLC", "EPS8", "SMURF1"))
cve$col <- "ns"
cve$col[cve$padj < 0.05] <- "<0.05"
cve$col[cve$padj < 0.01] <- "<0.01"
cve$col[cve$padj < 1e-5] <- "<1e-5"
cve$col[cve$padj < 1e-10] <- "<1e-10"
cve$col[cve$padj < 1e-20] <- "<1e-20"
cve$col <- factor(cve$col, levels = c("ns", "<0.05", "<0.01", "<1e-5", "<1e-10", "<1e-20"))

# fig 4b. 7 x 8.5
p <- ggplot(cve, aes(x = catype, y = cor, fill = col))
p + geom_bar(stat = "identity", size = 0.25, color = "white") +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = blue_pal[c(7,8,9,10,11,12)]) +
  facet_wrap(~ gene, nrow = 1) +
  theme(axis.text = element_text(size = 8, color = "black"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(linetype = "dashed", color = "gray89"))

# Suppl table 4: CBLC, EPS8 and SMURF1 segment mean vs. expression
# correlations (Spearman rho) in 24 cancer types
tmp <- cve[, c(1,2,4:5)]
st4 <- dcast(catype ~ gene, data = tmp, value.var = c("cor"))
tmp2 <- dcast(catype ~ gene, data = tmp, value.var = c("padj"))
colnames(st4)[2:4] <- paste0(colnames(st4)[2:4], "_cor")
colnames(tmp2)[2:4] <- paste0(colnames(tmp2)[2:4], "_padj")
st4 <- cbind(st4, tmp2[, -1])
st4 <- st4[order(st4$CBLC_cor, decreasing = T), ]
tmp <- data.frame(table(pancnv$catype[pancnv$gene == "CBLC"]), stringsAsFactors = F)
st4$n <- tmp$Freq[match(st4$catype, tmp$Var1)]
# write.table(st4[, c(1, 8, 2:7)], sep = '\t', quote = F, row.names = F,
#             file = "/projects/sftp/jtopham/incoming/fxbin_paper_supptable4.tsv")

# fig 4c
pancnv$tmp <- "non gain"
pancnv$tmp[pancnv$meanseg > 0.49 & pancnv$numprobes >= 10] <- "gain"

pandon <- data.frame(table(pancnv$tmp, pancnv$catype, pancnv$gene), 
                     stringsAsFactors = F)
colnames(pandon) <- c("cnv", "catype", "gene", "freq")
pandon$str <- paste0(pandon$catype, "_", pandon$gene)
pandon <- split(pandon, f = pandon$str)
pandon <- lapply(pandon, FUN = function(x){
  x$perc <- (x$freq / sum(x$freq)) * 100
  return(x)
})
pandon <- do.call(rbind, pandon)
pandon$str <- NULL

# upper bars
ubars <- NULL
for(i in unique(pancnv$catype)){
  for(j in unique(pancnv$gene)){
    tmp <- pancnv[pancnv$catype == i & pancnv$numprobes >= 10 & pancnv$gene == j, ]
    this.p <- tryCatch(wilcox.test(tmp$exp[tmp$tmp == "gain"],
                                   tmp$exp[tmp$tmp == "non gain"])$p.value,
                       error=function(err) NA)
    this.diff <- mean(tmp$exp[tmp$tmp == "gain"]) -
      mean(tmp$exp[tmp$tmp == "non gain"])
    ubars <- rbind(ubars, cbind(catype = i, gene = j, p = this.p, diff = this.diff))
  }
}
ubars <- data.frame(ubars, stringsAsFactors = F)
ubars$p <- as.numeric(ubars$p)
ubars$diff <- as.numeric(ubars$diff)

#
ord <- ubars[ubars$gene == "CBLC", ]
ord <- ord[order(ord$diff, decreasing = T), ]
ord$catype <- as.character(ord$catype)

ubars$catype <- factor(ubars$catype, levels = ord$catype)
ubars$gene <- factor(ubars$gene, levels = c("CBLC", "EPS8", "SMURF1"))
ubars$padj <- NA
for(i in unique(ubars$gene)){
  ubars$padj[ubars$gene == i] <- p.adjust(ubars$p[ubars$gene == i], method = "BH")
}
ubars$p2 <- "ns"
ubars$p2[ubars$padj < 0.05] <- "<0.05"
ubars$p2[ubars$padj < 0.01] <- "<0.01"
ubars$p2[ubars$padj < 1e-5] <- "<1e-5"
ubars$p2[ubars$padj < 1e-10] <- "<1e-10"
ubars$p2 <- factor(ubars$p2, levels = c("ns", "<0.05", "<0.01", "<1e-5", "<1e-10"))

# 15 x 9
p <- ggplot(ubars, aes(x = catype, y = diff, fill = p2))
p + geom_bar(stat = "identity", size = 0.4, color = "white") +
  scale_fill_manual(values = rev(purples_pal[c(12,10,8,6,4)])) +
  facet_wrap(~ gene, ncol = 1) +
  theme(axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8, 
                                   color = "black"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(linetype = "dashed", color = "gray89"))

# lower
pandon$catype <- factor(pandon$catype, levels = levels(ubars$catype))
pandon$gene <- factor(pandon$gene, levels = c("SMURF1", "EPS8", "CBLC"))
pandon$perc2 <- "0"
pandon$perc2[pandon$perc > 0 & pandon$perc < 2] <- "<2"
pandon$perc2[pandon$perc > 2 & pandon$perc < 4] <- "2-4"
pandon$perc2[pandon$perc > 4 & pandon$perc < 6] <- "4-6"
pandon$perc2[pandon$perc > 6 & pandon$perc < 8] <- "6-8"
pandon$perc2[pandon$perc > 8 & pandon$perc < 10] <- "8-10"
pandon$perc2[pandon$perc > 10 & pandon$perc < 20] <- "10-20"
pandon$perc2[pandon$perc > 20 & pandon$perc < 50] <- "20-50"
pandon$perc2[pandon$perc > 50] <- "50+"
pandon$perc2 <- factor(pandon$perc2, levels = c("0", "<2", "2-4", "4-6", "6-8", "8-10",
                                                "10-20", "20-50", "50+"))

# 15 x 4
p <- ggplot(pandon[pandon$cnv == "gain", ], aes(x = catype, y = gene, fill = perc2))
p + geom_tile(size = 0.4, color = "white") +
  scale_fill_manual(values = c("white", blue_pal[6], "#e4d7d9" ,
                               blue_pal[8:12], "gray15")) +
  theme(axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8, 
                                   color = "black"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(linetype = "dashed", color = "gray89"))

# Suppl table 5
st5 <- pandon[pandon$cnv == "gain", ]; st5$str <- paste0(st5$catype, "_", st5$gene)
st5$diff <- ubars$diff[match(st5$str, paste0(ubars$catype, "_", ubars$gene))]
st5$diff_padj <- ubars$padj[match(st5$str, paste0(ubars$catype, "_", ubars$gene))]
st5 <- st5[order(st5$diff, decreasing = T), ]
# write.table(st5[, c(2, 3, 5, 8:9)], sep = '\t', quote = F, row.names = F,
#             file = "/projects/sftp/jtopham/incoming/fxbin_paper_supptable5.tsv")

# fig 4d
# meta data
meta.t <- read.delim(paste0("/projects/pangen/analysis/jtopham/data/tcga/",
                            "complete_TCGA_metadata.tsv"),
                     stringsAsFactors = F, sep = '\t', header = T)
meta.t$death_days_to <- as.numeric(meta.t$death_days_to)
meta.t$last_contact_days_to <- as.numeric(meta.t$last_contact_days_to)
meta.t$days_to_last_followup <- as.numeric(meta.t$days_to_last_followup)

to <- data.frame(cbind(id = tcga.all2$sample_id, catype = tcga.all2$cancer_type),
                 stringsAsFactors = F)
to$id2 <- paste(do.call(rbind, strsplit(to$id, "-"))[, 1],
                do.call(rbind, strsplit(to$id, "-"))[, 2],
                do.call(rbind, strsplit(to$id, "-"))[, 3], sep = "-")
to$to_death <- meta.t$death_days_to[match(to$id2, meta.t$bcr_patient_barcode)]
to$last_follow <- meta.t$last_contact_days_to[match(to$id2, meta.t$bcr_patient_barcode)]

to <- to[!(is.na(to$to_death) & is.na(to$last_follow)), ]
to$censor <- 1
to$censor[is.na(to$to_death)] <- 0

to$days <- to$to_death
to$days[is.na(to$to_death)] <- to$last_follow[is.na(to$to_death)]
to$days <- to$days / 30.42
to <- to[to$days > 0, ]

# filter for CA types with at least 10% of samples not censored
filt <- table(to$censor, to$catype)[2, ] / table(to$censor, to$catype)[1, ]

# no PRAD, TGCT, THCA, THYM or UCEC
to <- to[!to$catype %in% names(filt)[filt < .1], ]

# catypes with at least 10 samples of interest
to$group <- to$id %in% pancnv$sample[pancnv$gene == "CBLC" & pancnv$tmp == "gain"]
to <- to[to$catype %in% 
         names(table(to$catype[to$group == T])[table(to$catype[to$group == T]) >= 10]), ]

# testing 10 ca types, 3927 samples
tcres <- NULL
for(t in unique(to$catype)){
  tmp <- to[to$catype == t, ]
  res.cox <- coxph(Surv(days, censor) ~ group, data = tmp)
  tcres <- rbind(tcres,
                 cbind(catype = t,
                       p = 1 - pchisq(survdiff(Surv(days, censor) ~ group, 
                           data = tmp)$chisq, 
                           length(survdiff(Surv(days, censor) ~ group, 
                           data = tmp)$n) - 1),
                       n.yes = sum(tmp$group == T),
                       n.no = sum(tmp$group == F),
                       med.yes = unname(summary(survfit(Surv(days, censor) ~ group, 
                                                        data = tmp))$table[,7])[2],
                       med.no = unname(summary(survfit(Surv(days, censor) ~ group, 
                                                       data = tmp))$table[,7])[1],
                       hr = data.frame(summary(res.cox)$conf.int, stringsAsFactors = F)[, 1],
                       ci95 = paste0(data.frame(summary(res.cox)$conf.int, stringsAsFactors = F)[,3],
                                     "-", data.frame(summary(res.cox)$conf.int, stringsAsFactors = F)[,4])))
}
tcres <- data.frame(tcres, stringsAsFactors = F)
for(i in 2:7){tcres[, i] <- as.numeric(tcres[, i])}

ggsurv(survfit(Surv(days, censor) ~ group, 
               data = to[to$catype=="CESC", ]),
       surv.col = c("#7289da", "#ffa31a"),
       size.est = 1.5, cens.size = 3.5, back.white = TRUE,
       order.legend = FALSE, xlab = "time (months)")
