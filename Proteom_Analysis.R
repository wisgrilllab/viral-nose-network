###################
#### Libraries ####
###################
library(tidyverse) 
library(edgeR)
library(limma)
library(ggfortify)
library(circlize)
library(ComplexHeatmap)
library(biomaRt)
library(Biobase)
library(clusterProfiler)
library(ReactomePA)

###################
#### Functions ####
###################

apply_tmm_factors <- function(y, color = NULL, plot = TRUE) {
  # computes the tmm normalized data from the DGEList object
  # y - DGEList object
  # returns a dataframe with normalized intensities
  
  # compute grand total (library size) scalings
  lib_facs <- mean(y$samples$lib.size) / y$samples$lib.size
  
  # the TMM factors are library adjustment factors (so divide by them)
  norm_facs <- lib_facs / y$samples$norm.factors
  cat("Overall Factors (lib.size+TMM):\n", sprintf("%-5s -> %f\n", colnames(y$counts), norm_facs))
  
  # compute the normalized data as a new data frame
  tmt_tmm <- as.data.frame(sweep(y$counts, 2, norm_facs, FUN = "*"))
  colnames(tmt_tmm) <- str_c(colnames(y$counts), "")
  
  # visualize results and return data frame
  if(plot == TRUE) {
    boxplot(log10(tmt_tmm), col = color, notch = TRUE, main = "TMM Normalized data")
  }
  tmt_tmm
}
CV <- function(df) {
  # Computes CVs of data frame rows
  # df - data frame, 
  # returns vector of CVs (%)
  ave <- rowMeans(df)    # compute averages
  sd <- apply(df, 1, sd) # compute standard deviations
  cv <- 100 * sd / ave   # compute CVs in percent (last thing gets returned)
}
labeled_boxplot <- function(df, ylim, title) {
  # Makes a box plot with the median value labeled
  # df - data frame with data to compute CVs of
  # ylim - upper limit for y-axis
  # title - plot title
  cv = CV(df)
  boxplot(cv, ylim = c(0, ylim), notch = TRUE, main = title)
  text(x = 0.65, y = boxplot.stats(cv)$stats[3], 
       labels = round(boxplot.stats(cv)$stats[3], 1))
}
plot_func <- function(plot_params) {
  labeled_boxplot(plot_params$data, ylim = plot_params$ylim, title = plot_params$title)
}

##################################################################
#### Initial Data Processing and different abundance analysis ####
##################################################################

#### Import raw data & create objects for downstream analysis ####
# This table has 12 samples, with the first column containing the Uniprot accession numbers
# Rows containing NAs have been removed
prot.raw = read.delim("proteome_raw.txt")
# Extract the accession numbers for later functions
accession <- data.frame(uniprot = prot.raw[,1])
# Keep only the reporter ion intensities
prot.tmt <- prot.raw[,-1]
# Import protein list for later annotation of the dataset
prot.list = read.csv("protein_list_final.csv", header=TRUE, sep=";")
# Create Metadata table
prot.meta <- tibble::tibble(group = factor(rep(c("Adult", "Term", "Preterm"), each = 4), levels = c("Adult", "Term", "Preterm"))) %>% 
             tibble::rownames_to_column("sample") %>% 
             mutate(sample = colnames(prot.tmt)) %>% 
             as.data.frame() %>% 
             tibble::column_to_rownames("sample")

#### TMM Normalization ####
 # Visulaize data before TMM Normalization
color.box = c(rep("#3087BC", 4), rep("#F7AB61", 4), rep("#D9464D", 4))
boxplot(log10(prot.tmt), col = color.box, notch = TRUE, main = "Starting Proteome Data before TMM")

# load data into DGEList object
group <- c(rep("Adult", 4), rep("Term", 4), rep("Preterm", 4))
y <- DGEList(counts = prot.tmt, group = group, genes = accession)

# Calculate Normalization factor
y <- calcNormFactors(y)

# Apply normalized values and plot afterwards
prot.norm <- apply_tmm_factors(y, color.box)

# Visualize CVs before and after TMM transformation
box_limit=150
plots <- list(
  list(data = prot.tmt[1:4], ylim = box_limit, title = "Starting Adult CVs"),
  list(data = prot.tmt[5:8], ylim = box_limit, title = "Starting Term CVs"),
  list(data = prot.tmt[9:12], ylim = box_limit, title = "Starting Preterm CVs"),
  list(data = prot.norm[1:4], ylim = box_limit, title = "Adult CVs after TMM"),
  list(data = prot.norm[5:8], ylim = box_limit, title = "Term CVs after TMM"),
  list(data = prot.norm[9:12], ylim = box_limit, title = "Preterm CVs after TMM")
)
# create plots
par(mfrow = c(2, 3))
lapply(plots, plot_func)
par(mfrow = c(1, 1))

#### Prepare final dataset ####

# log2 transformation and switch uniprot to HUGO symbols
prot.data <- cbind(uniprot = accession$uniprot, log2(prot.norm))
prot.data <- prot.data %>% 
             inner_join(prot.list, by = c("uniprot" = "accession")) %>%
             distinct(gene_symbol, .keep_all = TRUE) %>%
             dplyr::select(-uniprot) %>%  # This line removes the 'uniprot' column
             column_to_rownames("gene_symbol")

#### Differential Abundant Protein (DAP) Analysis ####

# Create an Expression Set     
pset <- ExpressionSet(assayData = as.matrix(prot.data), phenoData = AnnotatedDataFrame(prot.meta))

# Create design and contrast matrix
design <- model.matrix(~0+group, prot.meta)
contrast.matrix <- makeContrasts(T_vs_A = groupTerm-groupAdult,
                                 P_vs_A = groupPreterm-groupAdult,
                                 P_vs_T = groupPreterm-groupTerm,
                                 levels=design)
# Fit limma model
fit <- lmFit(pset, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend = TRUE)         
summary(decideTests(fit2, p.value = 0.05))

# Load Biomart database before continuing
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "https://useast.ensembl.org")

# Make group comparison, filter, and export DAP
p.vs.a <- topTable(fit2, coef = "P_vs_A", number = nrow(fit2), sort.by = "P") %>%
                   filter(adj.P.Val < 0.05) %>%
                   mutate(prot = rownames(.)) %>%
                   rowwise() %>%
                   mutate(entrez = {symbols = getBM(attributes = c('entrezgene_id', 'hgnc_symbol'), 
                                   filters = 'hgnc_symbol', 
                                   values = prot, 
                                   mart = mart) 
                                   symbols <- symbols %>% distinct(hgnc_symbol, .keep_all = TRUE)  # Keep only first mapping
                                   if (nrow(symbols) > 0) symbols$entrezgene[1] else NA  # If there are multiple mappings, keep the first one
                                   }) %>% as.data.frame()
#write.table(p.vs.a, "LIMMA_PRETERM_ADULT_PROT.txt", sep="\t", quote = F)

t.vs.a <- topTable(fit2, coef = "T_vs_A", number = nrow(fit2), sort.by = "P") %>%
                  filter(adj.P.Val < 0.05) %>%
                  mutate(prot = rownames(.)) %>%
                  rowwise() %>%
                  mutate(entrez = {symbols = getBM(attributes = c('entrezgene_id', 'hgnc_symbol'), 
                                  filters = 'hgnc_symbol', 
                                  values = prot, 
                                  mart = mart) 
                                  symbols <- symbols %>% distinct(hgnc_symbol, .keep_all = TRUE)  # Keep only first mapping
                                  if (nrow(symbols) > 0) symbols$entrezgene[1] else NA  # If there are multiple mappings, keep the first one
                                  }) %>% as.data.frame()

#write.table(t.vs.a, "LIMMA_TERM_ADULT_PROT.txt", sep="\t", quote = F)

#p.vs.t has no DAPs


#######################################
#### Data Analysis & Visualization ####
#######################################

#### PCA ####
#As seen in Figure 2A
prot.PCA <- prcomp(t(prot.data))
summary(prot.PCA)
autoplot(prot.PCA, label = F, loadings = F, shape = 'group', loadings.label = F,
         data = prot.meta, colour = 'group', size=3, frame=TRUE, frame.type = 'norm') +
         scale_color_manual(values = c("#3087BC", "#F7AB61", "#D9464D")) + 
         scale_fill_manual(values = c("#3087BC", "#F7AB61", "#D9464D"))


#### Heatmap ####
#As seen in Figure 2B
# Create logFC annotation for heatmap
fc.anno.dap <- full_join(data.frame(PvsA = p.vs.a$prot, PvsA.lfc = p.vs.a$logFC),
                     data.frame(TvsA = t.vs.a$prot, TvsA.lfc = t.vs.a$logFC),
                     by = c("PvsA" = "TvsA")) %>%
                     replace(is.na(.), 0) %>%
                     column_to_rownames("PvsA")

# Define color scheme
col.fun = colorRamp2(c(-2, 0, 2), c("#5E4FA2", "#FFFFBF", "#9E0142"))
col.fun2 = colorRamp2(c(-1, 0, 1), c("#00429D", "#EEEEEE", "#93003A"))

# Get DAP from the expression matrix
hm.dap = prot.data[match(rownames(fc.anno.dap),rownames(prot.data)),]

# Main Heatmap with Expression Values 
heat.dap = Heatmap(t(scale(t(hm.dap))), name = "Z-score", col = col.fun,
           show_row_names = F,row_dend_gp = gpar(col = "gray"), width = unit(6, "cm"),
           column_dend_gp = gpar(col = "gray"), show_column_names = F, column_split = c(rep("Adult",4), rep("Term", 4), rep("Preterm",4)),
           top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("#3087BC","#F7AB61","#D9464D"), col="white")))
)

#Annotation Heatmap for logFC
anno.lfc = Heatmap(as.matrix(fc.anno.dap), name = "log2FC", col = col.fun2,
           show_row_names = F, show_column_names = F, width = unit(1.5, "cm"),
           column_split = c("PvsA.lfc", "TvsA.lfc"), cluster_column_slices = F, cluster_columns = F,
           column_order = c("PvsA.lfc", "TvsA.lfc"), column_gap = unit(0.3, "mm"),
           top_annotation = HeatmapAnnotation(foo2 = anno_block(show_name=FALSE, gp = gpar(fill = c("#999999", "#CCCCCC"), col="white")))
)

# Add Heatmap Legends
heat.dap.lgd = list(Legend(labels = c("PvsA", "TvsA"), title = "Comparison",
                    legend_gp = gpar(fill = c("#999999", "#CCCCCC"), col="white"))
)

# Print Heatmap
draw(heat.dap+anno.lfc, ht_gap = unit(0.01, "cm"), annotation_legend_list = heat.dap.lgd, merge_legend = TRUE)

#### Venn Diagram ####
# As seen in Figure 2C & Supplemental Figure 2
# The Venn Diagram was drawn using the jvenn application: https://jvenn.toulouse.inrae.fr/app/index.html

#### Reactome Pathway Enrichment of DAPs ####
# Create data frame with Entrez_ID and corresponding logFC value
entrez.p.vs.a <- data.frame("entrez" = p.vs.a$entrez, "logFC" = p.vs.a$logFC)
entrez.t.vs.a <- data.frame("entrez" = t.vs.a$entrez, "logFC" = t.vs.a$logFC)

# Enrichment of DAPs which are higher in Adults compared the Preterms
reactome.p.vs.a.down = enrichPathway(filter(entrez.p.vs.a, logFC < 0)$entrez, readable = TRUE)
# Enrichment of DAPs which are higher in Preterms compared the Adults
reactome.p.vs.a.up = enrichPathway(filter(entrez.p.vs.a, logFC > 0)$entrez, readable = TRUE)
# Enrichment of DAPs which are higher in Adults compared the Terms
reactome.t.vs.a.down = enrichPathway(filter(entrez.t.vs.a, logFC < 0)$entrez, readable = TRUE)
# Enrichment of DAPs which are higher in Terms compared the Adults
reactome.t.vs.a.up = enrichPathway(filter(entrez.t.vs.a, logFC > 0)$entrez, readable = TRUE)

#### Enrichment Dot Plot ####
# As senn in Figure 2C & Supplemental Figure 2
ggplot(reactome.p.vs.a.up, showCategory = 7, 
       aes(Count, fct_reorder(Description, Count))) + 
       geom_segment(aes(xend=0, yend = Description)) +
       geom_point(aes(color=p.adjust, size = 0.3)) +
       scale_color_distiller(palette = "Spectral", direction = 1) +
       theme_minimal() + 
       xlab("DAP Ratio") +
       ylab(NULL) +
       theme(text = element_text(size=6), axis.text.x = element_text(size=6),axis.text.y = element_text(size=6))

ggplot(reactome.p.vs.a.down, showCategory = 7, 
       aes(Count, fct_reorder(Description, Count))) + 
       geom_segment(aes(xend=0, yend = Description)) +
       geom_point(aes(color=p.adjust, size = 0.3)) +
       scale_y_discrete(position = "left", labels = function(y) str_wrap(y, width = 40))+
       theme_minimal() + 
       xlab("DAP Ratio") +
       ylab(NULL) +
       theme(text = element_text(size=6), axis.text.x = element_text(size=6),axis.text.y = element_text(size=6))

ggplot(reactome.t.vs.a.up, showCategory = 7, 
       aes(Count, fct_reorder(Description, Count))) + 
       geom_segment(aes(xend=0, yend = Description)) +
       geom_point(aes(color=p.adjust, size = 0.3)) +
       scale_color_distiller(palette = "Spectral", direction = 1) +
       theme_minimal() + 
       xlab("DAP Ratio") +
       ylab(NULL) +
       theme(text = element_text(size=6), axis.text.x = element_text(size=6),axis.text.y = element_text(size=6))

ggplot(reactome.t.vs.a.down, showCategory = 7, 
       aes(Count, fct_reorder(Description, Count))) + 
       geom_segment(aes(xend=0, yend = Description)) +
       geom_point(aes(color=p.adjust, size = 0.3)) +
       scale_color_distiller(palette = "Spectral", direction = 1) +
       theme_minimal() + 
       xlab("DAP Ratio") +
       ylab(NULL) +
       theme(text = element_text(size=6), axis.text.x = element_text(size=6),axis.text.y = element_text(size=6))

#### Export dataset for integration models ####
write.table(prot.data, "protein_expression_matrix_NEC.txt", sep="\t", quote = F)



