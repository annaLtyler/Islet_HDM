library(here)
library(clusterProfiler)
library(fgsea)
library(pheatmap)
library(gprofiler2)
library(pathview)

tissue.cols <- as.matrix(read.delim(here("Data", "general", "tissue_colors.txt"), header = FALSE, row.names = 1))[,1]

transcript_loadings <- readRDS(here("Results", "DO", "High_Dim_Med_no_MGE", 
    "tissue_together-_-complete_mediation-germline_kinship", "Loadings_Transcripts_0.RDS"))
tissue.names <- names(transcript_loadings)

kegg.file <- here("Data", "general", "KEGG.Mouse.RDS")
if(!file.exists(kegg.file)){
    all.kegg <- download_KEGG("mmu", "KEGG", "kegg")
    saveRDS(all.kegg, kegg.file)
}else{
    all.kegg <- readRDS(kegg.file)
}

u_path <- gsub(" - Mus musculus (house mouse)", "", all.kegg[[2]][,"to"], fixed = TRUE)

# build GSEA list
path.locale <- lapply(u_path, function(x) grep(x, all.kegg[[2]][,2], fixed = TRUE))
path.id <- lapply(path.locale, function(x) all.kegg[[2]][x,][1])
path.idx <- lapply(path.id, function(x) which(all.kegg[[1]][,1] == unlist(x)[1]))
path.gene.id <- lapply(path.idx, function(x) all.kegg[[1]][x,2])
names(path.gene.id) <- u_path

#load gene conversion table
gene.table <- read.table(here("Data", "general", "mouse_gene_info.txt"), sep = "\t", header = TRUE)
path.gene.ensembl <- lapply(path.gene.id, function(x) gene.table[match(x, gene.table[,"entrezgene_id"]), "ensembl_gene_id"])


#search.term <- "mTOR"
#search.term <- "Valine, leucine and isoleucine degradation"
#search.term <- "Valine, leucine and isoleucine biosynthesis"
#search.term <- "Apelin signaling pathway"
#search.term <- "Fatty acid biosynthesis"

path_enrich <- matrix(NA, nrow = length(u_path), ncol = length(tissue.names))
rownames(path_enrich) <- u_path
colnames(path_enrich) <- tissue.names

all_min <- min(sapply(transcript_loadings, min))
all_max <- max(sapply(transcript_loadings, max))


#get enrichment scores for all pathways and all tissues
for(tx in 1:length(tissue.names)){

    sorted.tx.loadings <- transcript_loadings[[tx]][order(transcript_loadings[[tx]], decreasing = TRUE),1]
    gsea.enrich <- fgsea::fgsea(path.gene.ensembl, sorted.tx.loadings, scoreType = "std")
    norm.es <- as.numeric(as.matrix(gsea.enrich[,"NES"]))
    names(norm.es) <- gsea.enrich$pathway

    name.idx <- match(names(norm.es), rownames(path_enrich))
    path_enrich[name.idx,tx] <- norm.es

    #pdf("~/Desktop/NES.pdf", width = 10, height = 100)
    #par(mar = c(4,30, 2, 2))
    #barplot(sort(norm.es), las = 2, horiz = TRUE)
    #dev.off()
}

pdf("~/Desktop/NES_by_tissue.pdf", width = 9, height = 80)
par(mar = c(4,25, 2, 2))
for(tx in 1:length(tissue.names)){
    barplot(path_enrich[order(path_enrich[,tx], decreasing = FALSE),tx], 
        horiz = TRUE, las = 2, main = tissue.names[tx], 
        xlim = c(min(path_enrich), max(path_enrich)))
}
dev.off()

pdf("~/Desktop/NES_heat.pdf", width = 7, height = 50)
pheatmap(path_enrich, display_numbers = TRUE)
dev.off()



pdf("~/Desktop/Kegg_Pathway_Loadings.pdf", width = 12, height = 12)
par(mfrow = c(3,3))
for(i in 1:length(u_path)){
    report.progress(i, length(u_path))
    search.term <- u_path[i]

    path.locale <- grep(search.term, all.kegg[[2]][,2], fixed = TRUE)
    path.id <- all.kegg[[2]][path.locale,][1]
    path.name <- all.kegg[[2]][path.locale,][2]

    path.idx <- which(all.kegg[[1]][,1] == unlist(path.id))
    path.gene.id <- all.kegg[[1]][path.idx,2]
    path.gene.table <- gconvert(as.numeric(path.gene.id), "mmusculus", numeric_ns = "ENTREZGENE_ACC")

    tissue.path.loadings <- vector(mode = "list", length = length(transcript_loadings))
    names(tissue.path.loadings) <- tissue.names
    for(tx in 1:length(transcript_loadings)){
        gene.idx <- match(path.gene.table[,"target"], rownames(transcript_loadings[[tx]]))
        tissue.path.loadings[[tx]] <- transcript_loadings[[tx]][gene.idx]
    }

    #test <- aov_list(tissue.path.loadings)
    test <- lapply(tissue.path.loadings, function(x) t.test(x))
    pval_p <- sapply(test, function(x) x$p.value)
    path_p[i,] <- pval_p

    vioplot(tissue.path.loadings, main = paste0(search.term, "\n", length(path.gene.id), " genes"), 
        ylim = c(all_min, all_max), col = tissue.cols)
    text(x = 1:4, y = all_max, labels = paste0("p = ", signif(pval_p, 2)))
    abline(h = 0)
}
dev.off()

pdf("~/Desktop/Kegg_path_p.pdf", width = 7, height = 44)
pheatmap(signif(-log10(path_p), 2), display_numbers = TRUE)
dev.off()

#q is the size of the overlap between pathway and eQTL genes
#m is number of qtls
#n is number of non qtls
#k number of genes in the pathway

local.eqtl <- readRDS(here("Results", "DO", "Transcriptomes", "eQTL_coef_local.RDS"))

search.term <- "Thermogenesis"
#search.term <- "Phagosome"
#search.term <- "Metabolic pathways"
search.term <- "Citrate cycle (TCA cycle)"

lod.thresh = 8
num.eQTL <- sapply(local.eqtl, function(x) length(which(x$lod > lod.thresh)))

path.locale <- grep(search.term, all.kegg[[2]][,2], fixed = TRUE)
path.id <- all.kegg[[2]][path.locale,][1]
path.name <- all.kegg[[2]][path.locale,][2]

path.idx <- which(all.kegg[[1]][,1] == unlist(path.id))
path.gene.id <- all.kegg[[1]][path.idx,2]
path.gene.table <- gconvert(as.numeric(path.gene.id), "mmusculus", numeric_ns = "ENTREZGENE_ACC")

eqtl.overlap <- lapply(local.eqtl, function(x) path.gene.table[,"target"][which(path.gene.table[,"target"] %in% x$gene.id)])
num.overlap <- sapply(eqtl.overlap, length)

non.eqtl.overlap <- lapply(local.eqtl, function(x) setdiff(x$gene.id, path.gene.table[,"target"]))
num.non.overlap <- sapply(non.eqtl.overlap, length)


num.in.path <- nrow(path.gene.table)

test.p <- sapply(1:length(num.overlap), function(x) phyper(num.overlap[x], num.eQTL[x], 
    num.non.overlap[x], num.in.path))

num.overlap
test.p

#================================================
# KEGG pathway viewing
#================================================

#search.term = "mTOR"
#search.term = "Thermogenesis"
#search.term = "Phagosome"
#search.term <- "Citrate cycle (TCA cycle)"
#search.term <- "Oxidative phosphorylation"
#search.term <- "2-Oxocarboxylic acid metabolism"
#search.term <- "Lysosome"
#search.term <- "Valine, leucine and isoleucine degradation"
#search.term <- "Valine, leucine and isoleucine biosynthesis"
search.term <- "Spliceosome"

path.locale <- grep(search.term, all.kegg[[2]][,2], fixed = TRUE)
path.id <- all.kegg[[2]][path.locale,][1]
path.name <- all.kegg[[2]][path.locale,][2]

path.idx <- which(all.kegg[[1]][,1] == unlist(path.id))
path.gene.id <- all.kegg[[1]][path.idx,2]
path.gene.table <- gconvert(as.numeric(path.gene.id), "mmusculus", numeric_ns = "ENTREZGENE_ACC")

tissue.path.loadings <- vector(mode = "list", length = length(transcript_loadings))
names(tissue.path.loadings) <- tissue.names
for(tx in 1:length(transcript_loadings)){
    gene.idx <- match(path.gene.table[,"target"], rownames(transcript_loadings[[tx]]))
    gene.set <- transcript_loadings[[tx]][gene.idx]
    names(gene.set) <- path.gene.table[,"input"]
    gene.idx <- gene.idx[which(!is.na(gene.idx))]
    tissue.path.loadings[[tx]] <- gene.set
}

scale.factor <- max(sapply(transcript_loadings, function(x) max(abs(x), na.rm = TRUE)))

tissue <- "Islet"
tissue.idx <- which(tissue.names == tissue)
pv.out <- pathview(gene.data = tissue.path.loadings[[tissue.idx]]/scale.factor, 
    pathway.id = path.id, species = "mmu", out.suffix = paste(tissue, search.term, sep = "_"),
    kegg.dir = here("Results", "Kegg"))



