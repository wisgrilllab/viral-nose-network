
###################
#### Libraries ####
###################
library(tidyverse) 
library(DESeq2)
library(biomaRt)
library(genefilter)
library(ComplexHeatmap)
library(circlize)
library(msigdbr)
library(clusterProfiler)
library(ggcorrplot)
library(RColorBrewer)
library(EnhancedVolcano)
library(edgeR)
library(limma)

###################
#### Functions ####
###################
col_fun = colorRamp2(c(-2, 0, 2), c("#5E4FA2", "#FFFFBF", "#9E0142"))
deg.analysis <- function(group1, group2) {
  results.name <- paste("DEG", group1, group2, "deseq2.txt", sep="_")
  dds.current <- results(dds, contrast = c("overgroup", group1, group2), alpha = 0.05) %>% 
    as.data.frame() %>% subset(padj < 0.05) %>% 
    mutate(Gene.Name = gene.list[match(rownames(.), gene.list$Gene.ID),]$Gene.Name) %>% 
    rownames_to_column(var = "Gene.ID") %>% mutate(Gene.ID = sub("\\..*", "", Gene.ID)) %>% 
    left_join(genes.ensembl.entrez, by = c("Gene.ID" = "ensembl_gene_id")) %>% arrange(desc(log2FoldChange)) %>% 
    filter(log2FoldChange > 0.58 | log2FoldChange < -0.58)
  write.table(dds.current, results.name, sep="\t", quote = F, row.names = FALSE)
  return(dds.current)
}

########################################################################
#### Initial Data Processing and differentially expression analysis ####
########################################################################

#### Import raw data & create objects for downstream analysis ####
rna.raw <- as.matrix(read.delim("rna_raw.tsv", row.names = 1))
rna.meta <- read.delim("rna_metadata.txt", stringsAsFactors=TRUE)
gene.list <- read.delim("gene_list.txt")

# Create a DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = rna.raw, colData = rna.meta, design = ~ overgroup)
# Filter genes that are only expressed in >4 samples with at least 10 counts
dds <- dds[rowSums(counts(dds) >= 4) >= 10,]
dds <- DESeq(dds)
summary(results(dds))

# Create expression dataset
data.rna <- as.data.frame(assay(vst(dds, blind = F)))
data.rna$genes <- gene.list[match(rownames(data.rna), gene.list$Gene.ID),]$Gene.Name
data.rna <- data.rna[!duplicated(data.rna$genes),]
rownames(data.rna) <- data.rna$genes
data.rna$genes <- NULL
# Export gene expression matrix
#write.table(data.rna, "gene_expression_matrix_NEC.txt", sep="\t", quote = F)


#### Differentially Expressed Genes (DEG) Analysis ####

# Generate Biomart object
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "https://useast.ensembl.org")
genes.ensembl.entrez <- getBM(attributes = c('ensembl_gene_id', 'entrezgene_id'), mart = mart)

#### Unstimulated ####
deg.pvsa.unstim <- deg.analysis("Preterm_MOCK", "Adult_MOCK")
deg.tvsa.unstim <- deg.analysis("Term_MOCK", "Adult_MOCK")
deg.pvst.unstim <- deg.analysis("Preterm_MOCK", "Term_MOCK")

#### IAV vs. Unstimulated ####
deg.p.iav <- deg.analysis("Preterm_INF","Preterm_MOCK")
deg.t.iav <- deg.analysis("Term_INF", "Term_MOCK")
deg.a.iav <- deg.analysis("Adult_INF", "Adult_MOCK")

#### RSV vs. Unstimulated ####
deg.p.rsv <- deg.analysis("Preterm_RSV","Preterm_MOCK")
deg.t.rsv <- deg.analysis("Term_RSV", "Term_MOCK")
deg.a.rsv <- deg.analysis("Adult_RSV", "Adult_MOCK")

#### IAV between the groups ####
deg.pvsa.iav <- deg.analysis("Preterm_INF", "Adult_INF")
deg.tvsa.iav <- deg.analysis("Term_INF", "Adult_INF")
deg.pvst.iav <- deg.analysis("Preterm_INF", "Term_INF")

#### RSV between the groups ####
deg.pvsa.rsv <- deg.analysis("Preterm_RSV", "Adult_RSV")
deg.tvsa.rsv <- deg.analysis("Term_RSV", "Adult_RSV")
deg.pvst.rsv <- deg.analysis("Preterm_RSV", "Term_RSV")

#### RSV vs. IAV ####
deg.p.rsv.iav <- deg.analysis("Preterm_RSV","Preterm_INF")
deg.t.rsv.iav <- deg.analysis("Term_RSV", "Term_INF")
deg.a.rsv.iav <- deg.analysis("Adult_RSV", "Adult_INF")


# Using edgeR and limma-voom for DEG comparison
x <- DGEList(counts = rna.raw, group = rna.meta$overgroup, genes = gene.list[match(rownames(rna.raw), gene.list$Gene.ID), ])
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)
L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6

keep.exprs <- filterByExpr(x, group=x$samples$group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

design <- model.matrix(~0+group, data = x$samples)
colnames(design) <- gsub("group", "", colnames(design))
design

contr.matrix <- makeContrasts(
  deg.p.iav = Preterm_INF - Preterm_MOCK, 
  deg.t.iav = Term_INF - Term_MOCK, 
  deg.a.iav = Adult_INF - Adult_MOCK,
  deg.p.rsv = Preterm_RSV - Preterm_MOCK, 
  deg.t.rsv = Term_RSV - Term_MOCK, 
  deg.a.rsv = Adult_RSV - Adult_MOCK, 
  levels = colnames(design))
contr.matrix

par(mfrow=c(1,2))
v <- voom(x, design, plot=TRUE)
v
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
summary(decideTests(efit, p.value = 0.05, lfc = 0.58))

tfit <- treat(vfit, lfc=0)
dt <- decideTests(tfit, p.value = 0.05)
summary(dt)


#############################################
#### Downstream Analysis & Visualization ####
#############################################

# Heatmap of highly differentially expressed genes as in Figure 3A
eset <- ExpressionSet(assayData = as.matrix(data.rna), phenoData = AnnotatedDataFrame(rna.meta))
#condition <- eset$overgroup
data.sd <- rowSds(exprs(eset))
top.genes <- names(sort(data.sd, decreasing = TRUE))[1:200]
data.var <- eset[top.genes, ]
sh <- shorth(data.sd)
data.sh <- eset[data.sd >= sh, ]
dim(data.sh)
tt <- rowFtests(data.sh, data.sh$overgroup, var.equal = T)
tt$p.adj <- p.adjust(tt$p.value, method = "BH")
heat.sig <- data.sh[tt$p.adj <= 1e-11, ]
dim(heat.sig)


# Annotation setup
col_hm <- list('group' = c('Adult' = '#3087BC', 'Term' = '#F7AB61', 'Preterm' = '#D9464D'),
               'condition' = c('MOCK' = '#E6F598', 'INF' = '#D53E4F', 'RSV' = '#ABDDA4'))

colAnn <- HeatmapAnnotation(df = rna.meta[c('group', 'condition')],
                            which = 'col',
                            col = col_hm,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'))

ha <- rowAnnotation(mark_genes = anno_mark(at = c(9,14,62,71,72,73,74,76,82,84,94),labels = rownames(heat.sig)[c(9,14,62,71,72,73,74,76,82,84,94)]))

Heatmap(t(scale(t(exprs(heat.sig)))), 
        name = "Z-score", 
        col = col_fun,
        show_row_names = F,
        row_dend_gp = gpar(col = "gray"),
        show_column_names = F, 
        right_annotation = ha, 
        row_names_gp = gpar(fontsize = 6),
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "complete",
        top_annotation = colAnn
)

# The Venn Diagram in Figure 3B was drawn using the jvenn application: https://jvenn.toulouse.inrae.fr/app/index.html

# Overview of the differentially expressed genes as seen in Figure 3B
fc.rna.a.iav <- data.frame(AIAV = deg.a.iav$Gene.Name, AIAV_lfc = deg.a.iav$log2FoldChange)
fc.rna.t.iav <- data.frame(TIAV = deg.t.iav$Gene.Name, TIAV_lfc = deg.t.iav$log2FoldChange)
fc.rna.p.iav <- data.frame(PIAV = deg.p.iav$Gene.Name, PIAV_lfc = deg.p.iav$log2FoldChange)
fc.rna.a.rsv <- data.frame(ARSV = deg.a.rsv$Gene.Name, ARSV_lfc = deg.a.rsv$log2FoldChange)
fc.rna.t.rsv <- data.frame(TRSV = deg.t.rsv$Gene.Name, TRSV_lfc = deg.t.rsv$log2FoldChange)
fc.rna.p.rsv <- data.frame(PRSV = deg.p.rsv$Gene.Name, PRSV_lfc = deg.p.rsv$log2FoldChange)

fc.rna.deg <- full_join(fc.rna.a.iav, fc.rna.t.iav, by = c("AIAV" = "TIAV"))
fc.rna.deg <- full_join(fc.rna.deg, fc.rna.p.iav, by = c("AIAV" = "PIAV"))
fc.rna.deg <- full_join(fc.rna.deg, fc.rna.a.rsv, by = c("AIAV" = "ARSV"))
fc.rna.deg <- full_join(fc.rna.deg, fc.rna.t.rsv, by = c("AIAV" = "TRSV"))
fc.rna.deg <- full_join(fc.rna.deg, fc.rna.p.rsv, by = c("AIAV" = "PRSV"))
fc.rna.deg[is.na(fc.rna.deg)] <- 0
fc.rna.deg <- fc.rna.deg[!duplicated(fc.rna.deg$AIAV),]
rownames(fc.rna.deg) <- fc.rna.deg$AIAV
fc.rna.deg$AIAV=NULL

deg.plot <- fc.rna.deg %>% gather() %>%
            group_by(key) %>%  
            summarize(pos = count(value>0), neg = count(value<0)) %>%
            gather("comp", "dir", 2:3)
deg.plot$key <- factor(deg.plot$key, levels = c("AIAV_lfc", "TIAV_lfc", "PIAV_lfc", "ARSV_lfc", "TRSV_lfc", "PRSV_lfc"))

ggplot(deg.plot, aes(x = comp, y = key)) +
  geom_point(aes(color = comp, size = dir)) +
  scale_color_manual(values = c("#FDAE61", "#66C2A5")) +
  scale_size(range = c(0.5, 10)) +
  geom_text(aes(label = dir), size = 2) +
  theme_minimal() +
  theme(axis.title.y = element_blank(), legend.position = "none")


# Hallmark Enrichment Analysis of Virus Stimulation as seen in Figure 3C
lst.deg <- list(
  "Adult_IAV" = na.omit(deg.a.iav$entrezgene_id),
  "Term_IAV" = na.omit(deg.t.iav$entrezgene_id),
  "Preterm_IAV" = na.omit(deg.p.iav$entrezgene_id),
  "Adult_RSV" = na.omit(deg.a.rsv$entrezgene_id),
  "Term_RSV" = na.omit(deg.t.rsv$entrezgene_id),
  "Preterm_RSV" = na.omit(deg.p.rsv$entrezgene_id)
)

msig.hall <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)

hall.comp.all.deg <- compareCluster(lst.deg, fun="enricher", TERM2GENE=msig.hall)
hall.all <- hall.comp.all.deg@compareClusterResult[,c(1,3,7)]
hall.all <- hall.all[order(hall.all$p.adjust),]
hall.all <- as.data.frame(hall.all %>% pivot_wider(names_from = Cluster, values_from = p.adjust))
rownames(hall.all) <- substring(hall.all$Description, 10)
hall.all$Description <- NULL
hall.all <- hall.all[order(rowSums(is.na(hall.all)),decreasing=FALSE),]

ggcorrplot(hall.all[-13,c(5,6,4:1)], hc.order = F,
           outline.col = "white",
           ggtheme = ggplot2::theme_minimal,
           lab = F, lab_size = 3) + 
           scale_fill_gradientn("mean.value", colours = brewer.pal(9, "Spectral")) +
           theme(text = element_text(size=6), axis.text.x = element_text(size=6),axis.text.y = element_text(size=6))

# Heatmap of Hallmark enrichment Analysis as seen in Figure 3D
vec.hall <- data.frame(genes = c("MX1", "OAS1", "OASL","BATF2","IL7","IL15","CXCL11","IRF2","IRF7","IRF9",
                               "FAS","ATF3","CASP1","CASP3","CASP4","CASP7","CASP8","IL6","TNF","TNFSF10","IFNB1","IRF1",
                               "BATF","IRF4","IRF6","IL1RL1","IL1R2","IL10","IL10RA","IL18R1","CXCL10","SOCS1","SOCS2",
                               "DTX1","WNT5A","LFNG","NOTCH1","NOTCH3","PRKCA"),
                               cat = c(rep("IFN response", 10),rep("Apoptosis", 12), rep("IL2/STAT5 signaling", 11), rep("NOTCH signaling", 6)))

hm.hall <- data.rna[match(vec.hall$genes,rownames(data.rna)),]
hm.hall <- cbind(data.frame(t(hm.hall[,1:59])), "overgroup" = rna.meta$overgroup)
hm.hall.mod <- hm.hall %>% group_by(overgroup) %>% summarise_all(mean) %>% t() %>% as.data.frame()
colnames(hm.hall.mod) <- hm.hall.mod[1,]
hm.hall.mod <- hm.hall.mod[-1,]
hm.hall.mod <- as.matrix(sapply(hm.hall.mod, as.numeric))
rownames(hm.hall.mod) <- rownames(t(hm.hall[,-40]))

col_fun = colorRamp2(c(-1, 0, 1), c("#5E4FA2", "#FFFFBF", "#9E0142"))

Heatmap(scale(t(hm.hall.mod)), name = "Z-score", col = col_fun, column_split = vec.hall$cat,
        column_title_gp = gpar(fontsize = 8),
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 6),
        column_names_rot = 45)

# Volcano Plots as seen in Supplemental Figure 3
# jvenn was used to create the Venn diagrams in the Figure
EnhancedVolcano(deg.a.iav,
                lab = deg.a.iav$Gene.Name,
                selectLab = c('IFNL1', 'IFNL3', 'IFNB1', 'AIM2', 'IRF7','IRF9', 'ICAM4', 'CXCL13','STAT1', 'MX2'),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-10, 10),
                ylim = c(0, -log10(10e-75)),
                pCutoff = 0.05,
                FCcutoff = 0.58,
                axisLabSize = 8,
                pointSize = 1,
                labSize = 3.0,
                colAlpha = 1,
                legendPosition = 'none',
                legendLabSize = 8,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                arrowheads = FALSE,
                max.overlaps = 10)

EnhancedVolcano(deg.t.iav,
                lab = deg.t.iav$Gene.Name,
                selectLab = c('IFNL1', 'IFNL3', 'IFNB1', 'AIM2', 'IRF7','IRF9','STAT1', 'MX2','IL6','TLR7','IL10','MTATP6P1', 'NOTCH3', 'MT-ND5', 'MT-CYB', 'HDAC11'),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-10, 10),
                ylim = c(0, -log10(10e-75)),
                pCutoff = 0.05,
                FCcutoff = 0.58,
                axisLabSize = 8,
                pointSize = 1,
                labSize = 3.0,
                colAlpha = 1,
                legendPosition = 'none',
                legendLabSize = 8,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                arrowheads = FALSE,
                max.overlaps = 100)

EnhancedVolcano(deg.p.iav,
                lab = deg.p.iav$Gene.Name,
                selectLab = c('IFNL1', 'IFNL3', 'IFNB1', 'AIM2', 'IRF7','IRF9','STAT1', 'MX2','IL6','TLR7','IL10','IRAK2','IRAK3', 'IRF4','CASP5','IL20RA','ATP6V1E2'),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-10, 10),
                ylim = c(0, -log10(10e-75)),
                pCutoff = 0.05,
                FCcutoff = 0.58,
                axisLabSize = 8,
                pointSize = 1,
                labSize = 4.0,
                colAlpha = 1,
                legendPosition = 'none',
                legendLabSize = 8,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                arrowheads = FALSE,
                max.overlaps = 100)

EnhancedVolcano(deg.a.rsv,
                lab = deg.a.rsv$Gene.Name,
                selectLab = c('IFNL1', 'IFNL3', 'IFNB1', 'AIM2', 'IRF7','IRF9', 'ICAM4', 'CXCL13','L1CAM', 'MX2', 'VCAM1'),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-10, 10),
                ylim = c(0, -log10(10e-75)),
                pCutoff = 0.05,
                FCcutoff = 0.58,
                axisLabSize = 8,
                pointSize = 1,
                labSize = 3.0,
                colAlpha = 1,
                legendPosition = 'none',
                legendLabSize = 8,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                arrowheads = FALSE,
                max.overlaps = 10)

EnhancedVolcano(deg.t.rsv,
                lab = deg.t.rsv$Gene.Name,
                selectLab = c('IFNL1', 'IL6', 'IFNB1', 'AIM2', 'IRF7','IRF9', 'IL17C', 'HIF3A','CLEC11A', 'MX2', 'VCAM1'),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-10, 10),
                ylim = c(0, -log10(10e-75)),
                pCutoff = 0.05,
                FCcutoff = 0.58,
                axisLabSize = 8,
                pointSize = 1,
                labSize = 3.0,
                colAlpha = 1,
                legendPosition = 'none',
                legendLabSize = 8,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                arrowheads = FALSE,
                max.overlaps = 10)

EnhancedVolcano(deg.p.iav,
                lab = deg.p.iav$Gene.Name,
                selectLab = c('IFNL1', 'IL6', 'IFNB1', 'AIM2', 'IRF7','IRF9', 'IL17C', 'IL17D','CLEC11A', 'MX2', 'ICAM3'),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-10, 10),
                ylim = c(0, -log10(10e-75)),
                pCutoff = 0.05,
                FCcutoff = 0.58,
                axisLabSize = 8,
                pointSize = 1,
                labSize = 3.0,
                colAlpha = 1,
                legendPosition = 'none',
                legendLabSize = 8,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                arrowheads = FALSE,
                max.overlaps = 10)

# Heatmap of expressed cytokines a seen in Supplemental Figure 3E
cyt.vec = c("IL1A","IL1B","IL6","IL7","IL10","IL11","IL12A","IL15","IL16","IL17C","IL17D","IL18","IL19","IL23A","IL27","IL32","IL33","IL34","IL36G","IL37",
            "IFNB1","IFNL1","IFNL2","IFNL3","IFNE","CCL2","CCL3","CCL4","CCL5","CCL14","CCL15","CCL17","CCL18","CCL19","CCL20","CCL22","CCL24","CCL25","CCL26","CCL28",
            "CXCL1","CXCL2","CXCL3","CXCL5","CXCL6","CXCL8","CXCL9","CXCL10","CXCL11","CXCL13","CXCL14","CXCL16","CXCL17", "TNF", "TGFB2", "TGFB3", "TGFB1")
merge.cyto = data.rna[match(cyt.vec,rownames(data.rna)),]
hm.cyto = cbind(as.data.frame(t(merge.cyto)), rna.meta$overgroup)
hm.cyto = hm.cyto %>% group_by(rna.meta$overgroup) %>% summarise_all(mean) %>% t() %>% as.data.frame()
colnames(hm.cyto) = hm.cyto[1,]
hm.cyto = hm.cyto[-1,]
hm.cyto = as.matrix(sapply(hm.cyto, as.numeric))
rownames(hm.cyto) = cyt.vec

Heatmap((scale(t(hm.cyto))), name = "Z-score", col = col_fun,
        show_row_names = T, cluster_rows = T, show_column_names = T, #row_split = hm_hall$group,
        cluster_columns = T, column_km = 3,
        row_names_gp = gpar(fontsize = c(8)),
        column_names_gp = gpar(fontsize = c(7)),
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:4),
                                                            labels = c("Cluster 1", "Cluster 2", "Cluster 3"), 
                                                            labels_gp = gpar(col = "white", fontsize = 8))))

