###################
#### Libraries ####
###################

library(tidyverse) 
library(mixOmics)
library(ggpubr)
library(biomaRt)
library(clusterProfiler)
library(ReactomePA)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)

# Initate biomaRt object
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "https://useast.ensembl.org")

#### Import Proteomics & RNAseq datasets and identify overlap ####
X <- read.delim("protein_expression_matrix_NEC.txt") # Proteom data
Y <- read.delim("gene_expression_matrix_NEC.txt") # RNAseq data
M <- tibble::tibble(group = factor(rep(c("Adult", "Term", "Preterm"), each = 4), levels = c("Adult", "Term", "Preterm"))) %>% 
     tibble::rownames_to_column("sample") %>% 
     mutate(sample = colnames(X)) %>% 
     as.data.frame() %>% 
     tibble::column_to_rownames("sample")

# Identify overlapping proteins/genes and filter both dataframes
common_rows <- intersect(rownames(X), rownames(Y))
X <- X[common_rows,]
Y <- Y[common_rows,]

# Select the RNAseq samples for the model & adapt the colnames
rna.samples <- c("A1_RSV_M","A2_RSV_M","A3_RSV_M","A4_RSV_M","T1_RSV_M","T2_RSV_M","T3_RSV_M","T4_RSV_M","P1_RSV_M","P2_RSV_M","P3_RSV_M","P4_RSV_M")
Y <- Y %>% dplyr::select(all_of(rna.samples))
colnames(Y) = colnames(X)

# Check the two dataframes
list(dim(X), dim(Y))

#### Check correlation between Proteome & RNAseq data
X.corplot = as.data.frame(rowMeans(X))
Y.corplot = as.data.frame(rowMeans(Y))

#Illustration as seen in Figure 2D
ggplot(data = data.frame(x = X.corplot$`rowMeans(X)`, y = Y.corplot$`rowMeans(Y)`), aes(x = x, y = y)) +
  geom_point(color="darkgrey", size=0.1)+
  geom_smooth(method=lm, se=FALSE, color="firebrick")+
  stat_cor(method = "spearman")+
  theme_classic()

#### Generate sPLS-DA models ####
set.seed(1234)

# Calculate MAD and filter the dataset
rmads.X <- apply(X, 1, mad)
X = t(X[rmads.X>0.2,])
rmads.Y <- apply(Y, 1, mad)
Y = t(Y[rmads.Y>0.2,])

# Look at the PCA
pca.prot <- pca(X, ncomp = 10, center = TRUE, scale = TRUE)
pca.rna <- pca(Y, ncomp = 10, center = TRUE, scale = TRUE)
plot(pca.prot)
plot(pca.rna)

plotIndiv(pca.prot, comp = c(1, 2), 
          group = M$group, 
          ind.names = M$group, 
          legend = TRUE, title = 'Proteomics, PCA comp 1 - 2')

plotIndiv(pca.rna, comp = c(1, 2), 
          group = M$group, 
          ind.names = M$group, 
          legend = TRUE, title = 'RNAseq, PCA comp 1 - 2')

# Initial sPLS model
spls <- spls(X = X, Y = Y, ncomp = 2, mode = 'canonical')
# Parameter tuning
perf.spls <- perf(spls, validation = 'Mfold', folds = 10, nrepeat = 20, progressBar = TRUE, cpus = 8)
plot(perf.spls, criterion = 'Q2.total')

# Selecting the number of variables
# set range of test values for number of variables to use from X dataframe
list.keepX <- c(seq(50, 200, 10))
# set range of test values for number of variables to use from Y dataframe
list.keepY <- c(seq(50, 200, 10)) 

# Parameter tuning
# This takes a while - down there are the number of keepX & keepY
tune.spls <- tune.spls(X, Y, ncomp = 2,
                      test.keepX = list.keepX,
                      test.keepY = list.keepY,
                      nrepeat = 50, folds = 10, # use 10 folds
                      mode = 'canonical', measure = 'cor', limQ2=0.095,
                      progressBar = TRUE) 
plot(tune.spls)
optimal.keepX <- tune.spls$choice.keepX #keepX = c(180, 50)
optimal.keepY <- tune.spls$choice.keepY #keepY = c(80, 50)

# Final spls model
final.spls <- spls(X, Y, ncomp = 2, 
                  keepX = optimal.keepX,
                  keepY = optimal.keepY,
                  mode = "canonical")

# Plot as seen in Figure 2D
plotIndiv(final.spls, ind.names = FALSE, 
          rep.space = "XY-variate",
          group = M$group,
          pch = as.factor(M$group),
          col.per.group = color.mixo(1:3),
          title = 'sPLS model',
          legend = TRUE)

cim.res <- cim(final.spls, save = 'jpeg', name.save = 'PLS_CIM_image')

# Get EntrezIDs from identified proteins
spls.entrez <- getBM(attributes = c('entrezgene_id', 'hgnc_symbol'), 
                     filters = 'hgnc_symbol', 
                     values = data.frame("genes" = cim.res$row.names), 
                     mart = mart)

# ReactomePA enrichment & plotting
spls.reactome = enrichPathway(spls.entrez$entrezgene_id, readable = TRUE)

# Plot interesting pathways
ggplot(spls.reactome[c(18,22,40)], showCategory = 80, 
       aes(Count, fct_reorder(Description, Count))) + 
       geom_segment(aes(xend=0, yend = Description)) +
       geom_point(aes(color=p.adjust, size = 0.3)) +
       scale_color_distiller(palette = "Spectral", direction = 1) +
       theme_minimal() + 
       xlab("Gene Ratio") +
       ylab(NULL) +
       theme(text = element_text(size=6), axis.text.x = element_text(size=6),axis.text.y = element_text(size=6))

# Selected proteins enriched in hedgehog & ROBO signaling
prot.targets <- c("EFCAB7","GAS8","IFT122","IFT140","TTC21B","WDR19",
                  "RPL14","RPL23A","RPL24","RPL26", "RPL27A","RPL28","RPL3","RPL36A","RPL37A","RPL4","RPL6","RPL7",
                  "RPL7A","RPL8","RPS24","RPS6","RPS8", "RPS13")
prot.spls <- t(X[,prot.targets])

# Prepare data from protein expression data
prot.spls <-  prot.spls %>% as.data.frame() %>%
              rownames_to_column(var = "Protein") %>%
              pivot_longer(cols = -Protein, names_to = "Individual", values_to = "Expression") %>%
              mutate(Group = case_when(str_detect(Individual, "Adult") ~ "Adult",
                                       str_detect(Individual, "Term") ~ "Term",
                                       str_detect(Individual, "Preterm") ~ "Preterm")) %>%
              group_by(Protein, Group) %>%
              summarize(MeanExpression = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
              pivot_wider(names_from = Group, values_from = MeanExpression) %>% as.data.frame()
rownames(prot.spls) <- prot.spls$Protein
prot.spls$Protein <- NULL

# Define color scheme
col.fun = colorRamp2(c(-1, 0, 1), c("#5E4FA2", "#FFFFBF", "#9E0142"))
col.group = c("Adult" = "#3087BC", "Term"="#F7AB61","Preterm"="#D9464D")
# Generate Heatmap
ha.prot.spls = HeatmapAnnotation(group = anno_simple(colnames(prot.spls), col = col.group))
Heatmap(t(scale(t(prot.spls))),name = "Z-score", col = col.fun,
        row_dend_gp = gpar(col = "gray"), column_dend_gp = gpar(col = "gray"),
        top_annotation = ha.prot.spls)







