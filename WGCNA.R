###################
#### Libraries ####
###################

library(tidyverse) 
library(WGCNA)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(pheatmap)

#### Functions ####
options(stringsAsFactors = FALSE)
cor <- WGCNA::cor
col.fun = colorRamp2(c(-1, 0, 1), c("#5E4FA2", "#FFFFBF", "#9E0142"))
heat.palette <- colorRampPalette(c("blue", "white", "red"))(n = 255)
breaks <- seq(-1.5, 1.5, length.out = 256)
dataQC <- function(data) {
  # Checking for good samples and genes
  gsg = goodSamplesGenes(data, verbose = 3)
  
  if (!gsg$allOK) {
    if (sum(!gsg$goodGenes)>0) {
      print(paste("Removing genes:", paste(names(data)[!gsg$goodGenes], collapse = ", ")))
    }
    
    if (sum(!gsg$goodSamples)>0) {
      print(paste("Removing samples:", paste(rownames(data)[!gsg$goodSamples], collapse = ", ")))
    }
    
    data = data[gsg$goodSamples, gsg$goodGenes]
  }
  
  # Clustering and plotting the dendrogram
  sampleTree = hclust(dist(data), method = "average")
  
  sizeGrWindow(12,9)
  par(cex = 0.6)
  par(mar = c(0,4,2,0))
  
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", 
       cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
  
  return(data)
}
calcThreshold <- function(data) {
  powers = c(c(1:10), seq(from = 12, to=30, by=2))
  
  sft = pickSoftThreshold(data, powerVector = powers, verbose = 5, networkType = "signed hybrid", corFnc = "cor")
  
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2))
  cex1 = 0.9
  
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  abline(h=0.80,col="red")
  
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
}
calcMEs <- function(mat, mergedColors, trait){
  # Define numbers of genes and samples
  nGenes = ncol(mat);
  nSamples = nrow(mat);
  
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(mat, mergedColors)$eigengenes
  MEs = orderMEs(MEs0)
  moduleTraitCor = cor(MEs, trait, use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
  
  sizeGrWindow(10,6)
  
  # Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  
  par(mar = c(6, 8.5, 3, 3));
  
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(trait),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  
  # Return values for further use
  return(list(MEs = MEs, 
              moduleTraitCor = moduleTraitCor, 
              moduleTraitPvalue = moduleTraitPvalue))
  }
multiData.subset = function(multiData, colIndex = NULL, rowIndex = NULL){
  size = checkSets(multiData);
  if (is.null(colIndex)) colIndex = c(1:size$nGenes);
  if (is.null(rowIndex)) rowIndex = lapply(size$nSamples, function(n) {c(1:n)})
  if (length(rowIndex)!=size$nSets) 
    stop("If given, 'rowIndex' must be a list of the same length as 'multiData'.");
  out = list();
  for (set in 1:size$nSets)
    out[[set]] = list(data = multiData[[set]]$data[rowIndex[[set]], colIndex, drop = FALSE]);
  names(out) = names(multiData);
  out;
}
varNames = function(multiData){
  colnames(multiData[[1]]$data);
}
multiData.mapply = function(FUN, ..., MoreArgs = NULL, mdmaSimplify = TRUE, mda.doCollectGarbage = FALSE){
  dots = list(...);
  if (length(dots)==0) 
    stop("No arguments were specified. Please type ?multiData.mapply to see the help page.");
  dotLengths = sapply(dots, length);
  if (any(dotLengths!=dotLengths[1]))
    stop(spaste("All arguments to vectorize over must have the same length.\n", 
                "Scalar arguments should be put into the 'MoreArgs' argument.\n",
                "Note: lengths of '...' arguments are: ", paste(dotLengths, collapse = ", ")));
  nArgs = length(dots);
  res = list();
  isMultiSet = sapply(dots, isMultiData);
  
  FUN = match.fun(FUN);
  nSets = dotLengths[1];
  for (set in 1:nSets)
  {
    localArgs = list();
    for (arg in 1:nArgs)
      localArgs[[arg]] = if (isMultiSet[arg]) dots[[arg]] [[set]] $ data else dots[[arg]] [[set]];
    names(localArgs) = names(dots);
    res[[set]] = list(data = do.call(FUN, c(localArgs, MoreArgs)));
    if (mda.doCollectGarbage) collectGarbage();
  }
  
  names(res) = names(dots[[1]]);
  
  if (mdmaSimplify)
    return(multiData.simplify(res));
  
  return(res);
}
process.vec.heat <- function(vec, data.rna, meta_condition) {
  merge.data = data.rna[match(vec, rownames(data.rna)),]
  hm.data = cbind(as.data.frame(t(merge.data)), meta_condition)
  
  hm.data = hm.data %>% 
    group_by(meta_condition) %>% 
    summarise_all(mean) %>% 
    t() %>% 
    as.data.frame()
  
  colnames(hm.data) = hm.data[1,]
  hm.data = hm.data[-1,]
  hm.data = as.matrix(sapply(hm.data, as.numeric))
  rownames(hm.data) = vec
  
  return(hm.data)
}
create.heatmap <- function(matrix_data, col_fun) {
  return(Heatmap(t(scale(t(matrix_data))), 
                 name = "Z-score", col = col_fun,
                 row_names_gp = gpar(fontsize = 6),
                 show_row_names = T, cluster_rows = T, 
                 show_column_names = T, cluster_columns = F))
}
if(exists("SoftConnectivity")) rm(SoftConnectivity);
SoftConnectivity=function(datE, power=6,batchsize=1500,MinimumNoSamples=10) {
  no.genes=dim(datE)[[2]]
  no.samples=dim(datE)[[1]]
  sum1=function(x) sum(x,na.rm=T)
  k=rep(NA,no.genes)
  no.batches=as.integer(no.genes/ batchsize)
  if (no.batches>0) {
    for (i in 1:no.batches) {
      print(paste("batch number = ", i))
      index1=c(1:batchsize)+(i-1)* batchsize
      ad1=abs(cor(datE[,index1], datE,use="p"))^power
      ad1[is.na(ad1)]=0
      k[index1]=apply(ad1,1,sum1)
      # If fewer than MinimumNoSamples contain gene expression information for a given
      # gene, then we set its connectivity to 0.
      NoSamplesAvailable=apply(!is.na(datE[,index1]),2,sum)
      k[index1][NoSamplesAvailable< MinimumNoSamples]=NA
    } # end of for (i in 1:no.batches
  } # end of if (no.batches>0)...
  if (no.genes-no.batches*batchsize>0 ) {
    restindex=c((no.batches*batchsize+1):no.genes)
    ad1=abs(cor(datE[,restindex], datE,use="p"))^power
    ad1[is.na(ad1)]=0
    k[restindex]=apply(ad1,1,sum1)
    NoSamplesAvailable=apply(!is.na(datE[,restindex]),2,sum)
    k[restindex][NoSamplesAvailable< MinimumNoSamples]=NA
  } # end of if
  k
} # end of function


#### Import RNAseq dataset & preprocess for WGCNA ####
data.rna <- read.delim("gene_expression_matrix_NEC.txt") # RNAseq data
rna.meta <- read.delim("rna_metadata.txt", stringsAsFactors=TRUE)

# set up the dataset for Adult
rna.adult <- c("A1_INF_M","A1_INF_S","A2_INF_M","A2_INF_S","A3_INF_M","A3_INF_S","A4_INF_M","A4_INF_S","A5_INF_M","A5_INF_S",
               "A1_RSV_M","A1_RSV_S","A2_RSV_M","A2_RSV_S","A3_RSV_M","A3_RSV_S","A4_RSV_M","A4_RSV_S","A5_RSV_M","A5_RSV_S")
data.rna.adult <- data.rna %>% dplyr::select(all_of(rna.adult))
rna.meta$condition <- factor(rna.meta$condition, levels = c("MOCK", "INF", "RSV"))
rna.meta.adult <- rna.meta %>% dplyr::filter(row.names(.) %in% rna.adult) %>% arrange(condition)
data.rna.adult <- t(data.rna.adult[, row.names(rna.meta.adult)])
trait.adult = data.frame(
  MOCK = c(rep(1, 10), rep(0,10)),
  VIRUS = c(rep(0, 10), rep(1,10)),
  IAV = c(rep(0,10), rep(1,5), rep(0,5)),
  RSV = c(rep(0,15), rep(1,5)),
  row.names = rownames(rna.meta.adult))
# Check data quality and remove unwanted genes
data.rna.adult <- dataQC(data.rna.adult)

# set up the dataset for Term
rna.term <- c("T1_INF_M","T1_INF_S","T2_INF_M","T2_INF_S","T3_INF_M","T4_INF_M","T4_INF_S","T5_INF_M","T5_INF_S",
              "T1_RSV_M","T1_RSV_S","T2_RSV_M","T2_RSV_S","T3_RSV_M","T3_RSV_S","T4_RSV_M","T4_RSV_S","T5_RSV_M","T5_RSV_S")
data.rna.term <- data.rna %>% dplyr::select(all_of(rna.term))
rna.meta.term <- rna.meta %>% dplyr::filter(row.names(.) %in% rna.term) %>% arrange(condition)
data.rna.term <- t(data.rna.term[, row.names(rna.meta.term)])
trait.term = data.frame(
  MOCK = c(rep(1, 10), rep(0,9)),
  VIRUS = c(rep(0, 10), rep(1,9)),
  IAV = c(rep(0,10), rep(1,4), rep(0,5)),
  RSV = c(rep(0,14), rep(1,5)),
  row.names = rownames(rna.meta.term))
# Check data quality and remove unwanted genes
data.rna.term <- dataQC(data.rna.term)
#T1_INF_S is a outlier and will be removed
remove <- c("T1_INF_S")
data.rna.term <- data.rna.term[!(rownames(data.rna.term) %in% remove), ]
trait.term <- trait.term[!(rownames(trait.term) %in% remove), ]
rna.meta.term <- rna.meta.term[!(rownames(rna.meta.term) %in% remove), ]

# set up the dataset for Preterm
rna.preterm <- c("P1_INF_M","P1_INF_S","P2_INF_M","P2_INF_S","P3_INF_M","P3_INF_S","P4_INF_M","P4_INF_S","P5_INF_M","P5_INF_S",
                 "P1_RSV_M","P1_RSV_S","P2_RSV_M","P2_RSV_S","P3_RSV_M","P3_RSV_S","P4_RSV_M","P4_RSV_S","P5_RSV_M","P5_RSV_S")
data.rna.preterm <- data.rna %>% dplyr::select(all_of(rna.preterm))
rna.meta.preterm <- rna.meta %>% dplyr::filter(row.names(.) %in% rna.preterm) %>% arrange(condition)
data.rna.preterm <- t(data.rna.preterm[, row.names(rna.meta.preterm)])
trait.preterm = data.frame(
  MOCK = c(rep(1, 10), rep(0,10)),
  VIRUS = c(rep(0, 10), rep(1,10)),
  IAV = c(rep(0,10), rep(1,5), rep(0,5)),
  RSV = c(rep(0,15), rep(1,5)),
  row.names = rownames(rna.meta.preterm))
# Check data quality and remove unwanted genes
data.rna.preterm <- dataQC(data.rna.preterm)

# Create datasets with matching genes
common.columns <- intersect(colnames(data.rna.preterm), colnames(data.rna.term))
common.columns <- intersect(common.columns, colnames(data.rna.adult))

# Subset each dataset to only include the common columns
data.rna.preterm <- data.rna.preterm[, common.columns]
data.rna.term <- data.rna.term[, common.columns]
data.rna.adult <- data.rna.adult[, common.columns]



#### Create ADULT Network ####
# Choose soft-thresholding powers
calcThreshold(data.rna.adult)

# Calculate Network
net.adult <- blockwiseModules(data.rna.adult, power = 14, networkType = "signed hybrid", maxBlockSize = 30000,
                             corType = "pearson", TOMType = "signed", minModuleSize = 50, deepSplit = 2,
                             mergeCutHeight=0.25, numericLabels = TRUE, pamRespectsDendro = T, verbose=3)

table(net.adult$colors)
moduleLabels.adult <- net.adult$colors
mergedColors.adult <- labels2colors(net.adult$colors)
table(mergedColors.adult)

# Generate dendrogram with gene significance values as seen in Figure 4B
iav <- data.frame("IAV" = trait.adult$IAV)
rsv <- data.frame("RSV" = trait.adult$RSV)
GS.iav = as.numeric(cor(data.rna.adult, iav, use = "p"))
GS.rsv = as.numeric(cor(data.rna.adult, rsv, use = "p"))
GS.iavColor = numbers2colors(GS.iav, signed = T)
GS.rsvColor = numbers2colors(GS.rsv, signed = T)
datColors = data.frame(mergedColors.adult,GS.iavColor,GS.rsvColor)[net.adult$blockGenes[[1]],]
plotDendroAndColors(net.adult$dendrograms[[1]], 
                    colors = datColors,
                    groupLabels = c("Module colors","GS.iav","GS.rsv"), 
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = F, guideHang = 0.05)

# Calculate MEs & show Module-trait relationship Heatmap for Supplemental Figure 4
eigengenes.adult <- calcMEs(data.rna.adult, mergedColors.adult, trait.adult)

# Select modules of interest which show correleation in the heatmap
interesting.modules <- c("magenta", "green", "yellow", "salmon", "red", "cyan", "black", "tan")
AM.indices <- c("AM1", "AM2", "AM3", "AM4", "AM5", "AM6","AM7", "AM8")

# Summarize the MEs according to condition, filter interesting modules & plot them as line graph as seen in figure 4C
MEs.adult.sum <- eigengenes.adult[1] %>% as.data.frame() %>% 
                 mutate(group = rna.meta.adult$condition) %>% 
                 pivot_longer(cols = starts_with("ME"), names_to = "modules", values_to = "me") %>%
                 mutate(modules = str_replace(modules, "MEs.ME", ""),
                 modules = factor(modules, levels = interesting.modules)) %>% 
                 group_by(group, modules) %>% 
                 summarize(meanME = mean(me), .groups = 'drop')  %>% 
                 filter(modules %in% interesting.modules) %>% 
                 arrange(modules) %>% as.data.frame()
MEs.adult.sum$subgroup = c(rep("MOCK",6), rep("IAV/RSV", 3), rep("RSV", 15))
MEs.adult.sum$subgroup <- factor(MEs.adult.sum$subgroup, levels = c("MOCK", "IAV/RSV", "RSV"))

ggplot(MEs.adult.sum, aes(x=group, y=meanME, color=modules, group = modules)) +
  geom_line() +
  geom_hline(yintercept=0,col="grey")+
  scale_color_manual(values = interesting.modules) +
  facet_wrap(subgroup~.,ncol=1)+
  theme_minimal()+
  theme(text = element_text(size=6), axis.text.x = element_text(size=6, angle = 45, hjust=1),axis.text.y = element_text(size=6),
        legend.position = "none")

# Filter & summarize gene expression for Heatmap
heat.adult <- data.rna.adult %>% t() %>% as.data.frame() %>%
              mutate(module = mergedColors.adult) %>% na.omit() %>%
              filter(module %in% interesting.modules) %>%
              mutate(module = factor(module, levels = interesting.modules)) %>% arrange(module)
heat.adult.modules <- heat.adult$module
heat.adult.final <- heat.adult[,-21] %>% t() %>% as.data.frame() %>%
                    mutate(group = rna.meta.adult$condition) %>%
                    group_by(group) %>% summarise_all(mean) %>% t() %>% as.data.frame()
colnames(heat.adult.final) = heat.adult.final[1,]
heat.adult.final = heat.adult.final[-1,] %>% mutate_all(function(x) as.numeric(as.character(x)))

anno.adult <- data.frame("Module" = heat.adult.modules, row.names = rownames(heat.adult.final))
anno.adult$Module <- factor(anno.adult$Module, levels = interesting.modules)
color.heat.adult = list(
  Module = c(magenta="magenta", green="green", yellow="yellow", salmon="salmon", red="red", cyan="cyan", black="black", tan="tan")
)

pheatmap(heat.adult.final, scale="row",
         annotation_colors = color.heat.adult,
         color = heat.palette, breaks = breaks,
         annotation_row = anno.adult,
         cluster_rows = F, cluster_cols = F,
         show_rownames=F, show_colnames = T)


# Enrichment of significant modules
enrich.adult <- data.rna.adult %>% t() %>% as.data.frame() %>%
                mutate(module = mergedColors.adult) %>% na.omit() %>%
                mutate(entrez = mapIds(org.Hs.eg.db, rownames(.), 'ENTREZID', 'SYMBOL')) %>% na.omit()

list.enrich.adult <- lapply(seq_along(interesting.modules), function(i) {
    enrich.adult %>%
    filter(module == interesting.modules[i]) %>%
    pull(entrez)
})
names(list.enrich.adult) <- AM.indices
rea.enrich.adult = compareCluster(list.enrich.adult, fun="enrichPathway", readable=TRUE)
#Export enrichment results for Supplementals
write.table(rea.enrich.adult@compareClusterResult, "wgcna_adult_enrichment_reactome.txt", sep="\t", quote = F)

adult.enrich.vec = c("R-HSA-376176", "R-HSA-71291", #AM1
                     "R-HSA-8856825","R-HSA-8856828", #AM2
                     "R-HSA-913531","R-HSA-1169410", "R-HSA-168928","R-HSA-9013694", #AM3
                     "R-HSA-9612973","R-HSA-1428517", #AM4
                     "R-HSA-5628897", "R-HSA-72172", #AM5
                     "R-HSA-6798695", "R-HSA-195721") #AM7
dotplot(filter(rea.enrich.adult, ID %in% adult.enrich.vec), size=1, font.size = 6, label_format = 30)

#### Create TERM Network ####

# Choose soft-thresholding powers
calcThreshold(data.rna.term)

# Calculate Network
net.term <- blockwiseModules(data.rna.term, power = 12, networkType = "signed hybrid", maxBlockSize = 30000,
                             corType = "pearson", TOMType = "signed", minModuleSize = 50, deepSplit = 2,
                             mergeCutHeight=0.25, numericLabels = TRUE, pamRespectsDendro = T, verbose=3)

table(net.term$colors)
moduleLabels.term <- net.term$colors
mergedColors.term <- labels2colors(net.term$colors)
table(mergedColors.term)

# Generate dendrogram with gene significance values as seen in Figure 4B
iav <- data.frame("IAV" = trait.term$IAV)
rsv <- data.frame("RSV" = trait.term$RSV)
GS.iav = as.numeric(cor(data.rna.term, iav, use = "p"))
GS.rsv = as.numeric(cor(data.rna.term, rsv, use = "p"))
GS.iavColor = numbers2colors(GS.iav, signed = T)
GS.rsvColor = numbers2colors(GS.rsv, signed = T)
datColors = data.frame(mergedColors.term,GS.iavColor,GS.rsvColor)[net.term$blockGenes[[1]],]
plotDendroAndColors(net.term$dendrograms[[1]], 
                    colors = datColors,
                    groupLabels = c("Module colors","GS.iav","GS.rsv"), 
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = F, guideHang = 0.05)

# Calculate MEs & show Module-trait relationship Heatmap for Supplemental Figure 4
eigengenes.term <- calcMEs(data.rna.term, mergedColors.term, trait.term)

# Select modules of interest which show correleation in the heatmap
interesting.modules <- c("greenyellow", "grey60", "salmon", "black", "pink", "blue", "lightgreen", "darkturquoise", "turquoise", "red", "darkred", "magenta")
TM.indices <- c("TM1", "TM2", "TM3", "TM4", "TM5", "TM6", "TM7", "TM8", "TM9", "TM10", "TM11", "TM12")

# Summarize the MEs according to condition, filter interesting modules & plot them as line graph as seen in figure 4C
MEs.term.sum <- eigengenes.term[1] %>% as.data.frame() %>% 
                mutate(group = rna.meta.term$condition) %>% 
                pivot_longer(cols = starts_with("ME"), names_to = "modules", values_to = "me") %>%
                mutate(modules = str_replace(modules, "MEs.ME", ""),
                modules = factor(modules, levels = interesting.modules)) %>% 
                group_by(group, modules) %>% 
                summarize(meanME = mean(me), .groups = 'drop')  %>% 
                filter(modules %in% interesting.modules) %>% 
                arrange(modules) %>% as.data.frame()
MEs.term.sum$subgroup = c(rep("MOCK",12), rep("IAV", 9), rep("RSV", 15))
MEs.term.sum$subgroup <- factor(MEs.term.sum$subgroup, levels = c("MOCK", "IAV", "RSV"))

ggplot(MEs.term.sum, aes(x=group, y=meanME, color=modules, group = modules)) +
  geom_line() +
  geom_hline(yintercept=0,col="grey")+
  scale_color_manual(values = interesting.modules) +
  facet_wrap(subgroup~.,ncol=1)+
  theme_minimal()+
  theme(text = element_text(size=6), axis.text.x = element_text(size=6, angle = 45, hjust=1),axis.text.y = element_text(size=6),
        legend.position = "none")

# Filter & summarize gene expression for Heatmap
heat.term <- data.rna.term %>% t() %>% as.data.frame() %>%
             mutate(module = mergedColors.term) %>% na.omit() %>%
             filter(module %in% interesting.modules) %>%
             mutate(module = factor(module, levels = interesting.modules)) %>% arrange(module)
heat.term.modules <- heat.term$module
heat.term.final <- heat.term[,-19] %>% t() %>% as.data.frame() %>%
                   mutate(group = rna.meta.term$condition) %>%
                   group_by(group) %>% summarise_all(mean) %>% t() %>% as.data.frame()
colnames(heat.term.final) = heat.term.final[1,]
heat.term.final = heat.term.final[-1,] %>% mutate_all(function(x) as.numeric(as.character(x)))

anno.term <- data.frame("Module" = heat.term.modules, row.names = rownames(heat.term.final))
anno.term$Module <- factor(anno.term$Module, levels = interesting.modules)
color.heat.term = list(
  Module = c(greenyellow = "greenyellow", grey60 = "grey60", salmon = "salmon", black = "black", pink = "pink", blue = "blue", lightgreen = "lightgreen", darkturquoise="darkturquoise", turquoise = "turquoise", red = "red", darkred = "darkred", magenta = "magenta")
)

pheatmap(heat.term.final, scale="row",
         annotation_colors = color.heat.term,
         color = heat.palette, breaks = breaks,
         annotation_row = anno.term,
         cluster_rows = F, cluster_cols = F,
         show_rownames=F, show_colnames = T)

# Enrichment of significant modules
enrich.term <- data.rna.term %>% t() %>% as.data.frame() %>%
               mutate(module = mergedColors.term) %>% na.omit() %>%
               mutate(entrez = mapIds(org.Hs.eg.db, rownames(.), 'ENTREZID', 'SYMBOL')) %>% na.omit()

list.enrich.term <- lapply(seq_along(interesting.modules), function(i) {
                    enrich.term %>%
                    filter(module == interesting.modules[i]) %>%
                    pull(entrez)})

names(list.enrich.term) <- TM.indices
rea.enrich.term = compareCluster(list.enrich.term, fun="enrichPathway", readable=TRUE)
#Export enrichment results for Supplementals
write.table(rea.enrich.term@compareClusterResult, "wgcna_term_enrichment_reactome.txt", sep="\t", quote = F)

term.enrich.vec = c("R-HSA-1474244", "R-HSA-1474290", #TM1
                    "R-HSA-5668914", #TM2
                    "R-HSA-8951664","R-HSA-169911", #TM5
                    "R-HSA-168255", "R-HSA-168273", #TM6
                    "R-HSA-191273", #TM7
                    "R-HSA-1445148", #TM8
                    "R-HSA-913531", "R-HSA-168898", #TM9
                    "R-HSA-400206", "R-HSA-1989781", #TM10
                    "R-HSA-9031628","R-HSA-198725", #TM11
                    "R-HSA-194138", "R-HSA-354192") #TM12
dotplot(filter(rea.enrich.term, ID %in% term.enrich.vec), size=1, font.size = 6, label_format = 30)

#### Create PRETERM Network ####
# Choose soft-thresholding powers
calcThreshold(data.rna.preterm)

# Calculate Network
net.preterm <- blockwiseModules(data.rna.preterm, power = 9, networkType = "signed hybrid", maxBlockSize = 30000,
                                corType = "pearson", TOMType = "signed", minModuleSize = 50, deepSplit = 2,
                                mergeCutHeight=0.25, numericLabels = TRUE, pamRespectsDendro = T, verbose=3)

table(net.preterm$colors)
moduleLabels.preterm <- net.preterm$colors
mergedColors.preterm <- labels2colors(net.preterm$colors)
table(mergedColors.preterm)

# Generate dendrogram with gene significance values as seen in Figure 4B
iav <- data.frame("IAV" = trait.preterm$IAV)
rsv <- data.frame("RSV" = trait.preterm$RSV)
GS.iav = as.numeric(cor(data.rna.preterm, iav, use = "p"))
GS.rsv = as.numeric(cor(data.rna.preterm, rsv, use = "p"))
GS.iavColor = numbers2colors(GS.iav, signed = T)
GS.rsvColor = numbers2colors(GS.rsv, signed = T)
datColors = data.frame(mergedColors.preterm,GS.iavColor,GS.rsvColor)[net.preterm$blockGenes[[1]],]
plotDendroAndColors(net.preterm$dendrograms[[1]], 
                    colors = datColors,
                    groupLabels = c("Module colors","GS.iav","GS.rsv"), 
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = F, guideHang = 0.05)

# Calculate MEs & show Module-trait relationship Heatmap for Supplemental Figure 4
eigengenes.preterm <- calcMEs(data.rna.preterm, mergedColors.preterm, trait.preterm)

# Select modules of interest which show correleation in the heatmap
interesting.modules <- c("pink", "magenta", "yellow", "darkgrey", "greenyellow", "brown", "tan", "lightgreen")
PM.indices <- c("PM1", "PM2", "PM3", "PM4", "PM5", "PM6", "PM7", "PM8")

# Summarize the MEs according to condition, filter interesting modules & plot them as line graph as seen in figure 4C
MEs.preterm.sum <- eigengenes.preterm[1] %>% as.data.frame() %>% 
                   mutate(group = rna.meta.preterm$condition) %>% 
                   pivot_longer(cols = starts_with("ME"), names_to = "modules", values_to = "me") %>%
                   mutate(modules = str_replace(modules, "MEs.ME", ""),
                   modules = factor(modules, levels = interesting.modules)) %>% 
                   group_by(group, modules) %>% 
                   summarize(meanME = mean(me), .groups = 'drop')  %>% 
                   filter(modules %in% interesting.modules) %>% 
                   arrange(modules) %>% as.data.frame()
MEs.preterm.sum$subgroup = c(rep("MOCK",9), rep("IAV", 12), rep("RSV", 3))
MEs.preterm.sum$subgroup <- factor(MEs.preterm.sum$subgroup, levels = c("MOCK", "IAV", "RSV"))

ggplot(MEs.preterm.sum, aes(x=group, y=meanME, color=modules, group = modules)) +
  geom_line() +
  geom_hline(yintercept=0,col="grey")+
  scale_color_manual(values = interesting.modules) +
  facet_wrap(subgroup~.,ncol=1)+
  theme_minimal()+
  theme(text = element_text(size=6), axis.text.x = element_text(size=6, angle = 45, hjust=1),axis.text.y = element_text(size=6),
        legend.position = "none")

# Filter & summarize gene expression for Heatmap
heat.preterm <- data.rna.preterm %>% t() %>% as.data.frame() %>%
                mutate(module = mergedColors.preterm) %>% na.omit() %>%
                filter(module %in% interesting.modules) %>%
                mutate(module = factor(module, levels = interesting.modules)) %>% arrange(module)
heat.preterm.modules <- heat.preterm$module
heat.preterm.final <- heat.preterm[,-21] %>% t() %>% as.data.frame() %>%
                      mutate(group = rna.meta.preterm$condition) %>%
                      group_by(group) %>% summarise_all(mean) %>% t() %>% as.data.frame()
colnames(heat.preterm.final) = heat.preterm.final[1,]
heat.preterm.final = heat.preterm.final[-1,] %>% mutate_all(function(x) as.numeric(as.character(x)))

anno.preterm <- data.frame("Module" = heat.preterm.modules, row.names = rownames(heat.preterm.final))
anno.preterm$Module <- factor(anno.preterm$Module, levels = interesting.modules)
color.heat.preterm = list(
  Module = c(pink="pink", magenta="magenta", yellow="yellow", darkgrey="darkgrey", greenyellow="greenyellow", brown="brown", tan="tan", lightgreen="lightgreen")
)

pheatmap(heat.preterm.final, scale="row",
         annotation_colors = color.heat.preterm,
         color = heat.palette, breaks = breaks,
         annotation_row = anno.preterm,
         cluster_rows = F, cluster_cols = F,
         show_rownames=F, show_colnames = T)

# Enrichment of significant modules
enrich.preterm <- data.rna.preterm %>% t() %>% as.data.frame() %>%
                  mutate(module = mergedColors.preterm) %>% na.omit() %>%
                  mutate(entrez = mapIds(org.Hs.eg.db, rownames(.), 'ENTREZID', 'SYMBOL')) %>% na.omit()

list.enrich.preterm <- lapply(seq_along(interesting.modules), function(i) {
                       enrich.preterm %>%
                       filter(module == interesting.modules[i]) %>%
                       pull(entrez)})

names(list.enrich.preterm) <- PM.indices
rea.enrich.preterm = compareCluster(list.enrich.preterm, fun="enrichPathway", readable=TRUE)
write.table(rea.enrich.preterm@compareClusterResult, "rea_enrich_wgcna_preterm.txt", quote = F, sep="\t")

preterm.enrich.vec = c("R-HSA-379726", #PM2
                       "R-HSA-2426168", "R-HSA-1989781","R-HSA-400206", #PM4
                       "R-HSA-71387", "R-HSA-71336",#PM5
                       "R-HSA-913531", "R-HSA-168898","R-HSA-168928", #PM6
                       "R-HSA-73894", "R-HSA-72163", "R-HSA-5676590" ) #PM7
dotplot(filter(rea.enrich.preterm, ID %in% preterm.enrich.vec), size=1, font.size = 6, label_format = 30)


#### Consensus Module Analysis ####

# Assemble datasets
nSets = 3
setLabels = c("Preterm - IAV", "Term - IAV", "Adult - IAV", "Preterm - RSV", "Term - RSV", "Adult - RSV")
shortLabels = c("Preterm", "Term", "Adult")
multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = as.data.frame(data.rna.preterm))
multiExpr[[2]] = list(data = as.data.frame(data.rna.term))
multiExpr[[3]] = list(data = as.data.frame(data.rna.adult))
exprSize = checkSets(multiExpr)
gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK

# Define data set dimensions
nGenes = exprSize$nGenes;
nSamples = exprSize$nSamples;

# Assemble traits
allTraits = rbind(trait.preterm, trait.term, trait.adult)
Traits = vector(mode="list", length = nSets);
for (set in 1:nSets){
  setSamples = rownames(multiExpr[[set]]$data);
  traitRows = match(setSamples, rownames(allTraits));
  Traits[[set]] = list(data = allTraits[traitRows,]);
}
varNames(Traits)
gsg = goodSamplesGenesMS(Traits);
multiTraitsQC = multiData.subset(Traits, gsg$goodGenes);
nTraits = checkSets(multiTraitsQC)$nGenes;
traitNames = varNames(multiTraitsQC);

# Build consensus network with previously calculated powers
net = blockwiseConsensusModules(multiExpr, power = c(9,12,14), maxBlockSize = 30000,
                                networkType = "signed hybrid", corType = "pearson",
                                TOMDenom = "mean", minModuleSize = 50, deepSplit = 2,
                                mergeCutHeight=0.25, consensusQuantile = 0.25,
                                numericLabels = TRUE, pamRespectsDendro = F,
                                verbose=3)
table(net$colors)
consMEs = net$multiMEs
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
consTree = net$dendrograms[[1]]

# Calculate gene significance and plot dendrogram as seen in figure 6A
MEs = multiSetMEs(multiExpr, universalColors = moduleLabels, excludeGrey= TRUE);
nModules = ncol(MEs[[1]]$data)
GS.iav = p.GS.iav = matrix(0, nGenes, nSets);
MESignif.iav = p.MESignif.iav = matrix(0, nModules, nSets);

for (set in 1:nSets)
{
  bcp = corAndPvalue(multiExpr[[set]]$data, Traits[[set]]$data$IAV)
  GS.iav[, set] = bcp$cor
  p.GS.iav[, set] = bcp$p;
  bcp = corAndPvalue(MEs[[set]]$data, Traits[[set]]$data$IAV);
  MESignif.iav[, set] = bcp$cor;
  p.MESignif.iav[, set] = bcp$p;
}
GS.rsv = p.GS.rsv = matrix(0, nGenes, nSets);
MESignif.rsv = p.MESignif.rsv = matrix(0, nModules, nSets);
for (set in 1:nSets)
{
  bcp = corAndPvalue(multiExpr[[set]]$data, Traits[[set]]$data$RSV)
  GS.rsv[, set] = bcp$cor
  p.GS.rsv[, set] = bcp$p;
  bcp = corAndPvalue(MEs[[set]]$data, Traits[[set]]$data$RSV);
  MESignif.rsv[, set] = bcp$cor;
  p.MESignif.rsv[, set] = bcp$p;
}
GS = cbind(GS.iav, GS.rsv)

sizeGrWindow(13,9);
plotDendroAndColors(consTree,
                    cbind(numbers2colors(GS,signed = TRUE), labels2colors(moduleLabels)),
                    colorHeight = 0.5,
                    dendroLabels = FALSE,
                    main = paste("NEC Consensus modules"),
                    addGuide = TRUE,
                    hang = 0.01, guideHang = 0.04,
                    marAll = c(0.2, 8, 1.5, 2.2))

# Metaanalysis calculating Z and p value as seen in Figure 6B
MEsignif0 = multiData.mapply(corAndPvalue, multiTraitsQC, MEs, mdmaSimplify = FALSE);
MEsignif = pMEsignif = array(0, dim = c(nSets, nTraits, nModules))
for (set in 1:nSets)
{
  MEsignif[set, , ] = MEsignif0[[set]]$data$cor;
  pMEsignif[set, , ] = MEsignif0[[set]]$data$p;
}

MEMAs = list()
colName = "Z.RootDoFWeights"
for (t in 1:nTraits)
{
  printFlush(paste("Working on", traitNames[t]));
  MEMAs[[t]] = metaAnalysis(MEs, multiData.subset(multiTraitsQC, t), useRankPvalue = FALSE,
                            corFnc = cor)
}

labLevels = sort(unique(moduleLabels[moduleLabels!=0]));
modSizes = table(moduleLabels[moduleLabels!=0]);
xLabels = spaste("ME", labels2colors(labLevels));
mat = MEsignif
dim(mat) = c(nSets*nTraits, nModules);
textMat = spaste(round(MEsignif, 2), "\n", signif(pMEsignif, 1));
dim(textMat) = dim(mat);
max = max(abs(MEsignif), na.rm = TRUE)

for (t in 1:nTraits)
{
  index = ( (t-1)*nSets+1):(t*nSets);
  Z = MEMAs[[t]] [ , match(colName, colnames(MEMAs[[t]]))];
  Z1 = Z/max(abs(Z)) * max;
  mat1 = rbind(Z1, mat[index, ]);
  pColName = sub("Z", "p", colName);
  p = MEMAs[[t]] [ , match(pColName, colnames(MEMAs[[t]]))];
  text1 = rbind(spaste(signif(Z, 2), "\n", signif(p, 1)), textMat[index, ]);
  yLabels1 = c("meta-analysis Z and p", shortLabels);
  par(mar = c(7, 13, 2.2, 0.3));
  labeledHeatmap(t(mat1),
                 xLabels = yLabels1,
                 yLabels = xLabels,
                 ySymbols = xLabels,
                 textMat = t(text1),
                 cex.text = 0.5,
                 main = spaste("Module eigengene significance for ", varNames(Traits)[t]),
                 colors = blueWhiteRed(100),
                 zlim = c(-max, max), setStdMargins = FALSE);
}

# Enrichment of significant consensus modules
probes <- colnames(data.rna.preterm)
geneInfo <- data.frame(geneSymbol = probes, moduleColor = labels2colors(moduleLabels))
entrez <- mapIds(org.Hs.eg.db, probes, 'ENTREZID', 'SYMBOL')
entrez <- as.data.frame(entrez)
geneInfo$entrez <- entrez$entrez
geneInfo <- na.omit(geneInfo)
cons.interesting.module <- c("red","tan","lightcyan","grey60", #mock
                             "brown", "lightgreen", #IAV
                             "pink", "purple", "greenyellow", "salmon", "midnightblue") 

geneInfo <- data.rna.preterm %>% 
            colnames() %>% 
            data.frame(geneSymbol = ., moduleColor = labels2colors(moduleLabels)) %>% 
            cbind(entrez = mapIds(org.Hs.eg.db, .$geneSymbol, 'ENTREZID', 'SYMBOL')) %>%
            filter(moduleColor %in% cons.interesting.module) %>%
            mutate(moduleColor = factor(moduleColor, levels = cons.interesting.module)) %>%
            na.omit()

rea.cons <- compareCluster(entrez~moduleColor, data=geneInfo, fun="enrichPathway", readable=TRUE)
rea.cons.filt <- rea.cons@compareClusterResult
rea.cons.filt$Cluster <- factor(rea.cons.filt$Cluster, levels = cons.interesting.module)
write.table(rea.cons.filt, "rea_cons_filt.txt", sep="\t", quote = F)

rea.vect <- c("R-HSA-200425", #Mock Red
              "R-HSA-72766", "R-HSA-71291",#Mock Tan
              "R-HSA-913531","R-HSA-168164","R-HSA-168188","R-HSA-166166","R-HSA-168181","R-HSA-168638", #brown
              "R-HSA-191273","R-HSA-1989781","R-HSA-8957322","R-HSA-6798695","R-HSA-8951664", #lightgreen
              "R-HSA-2990846", "R-HSA-3108232", #pink
              "R-HSA-195721", "R-HSA-9634600", #salmon
              "R-HSA-72086", "R-HSA-6781827") #midnightblue


rea.vect <- c("R-HSA-72766","R-HSA-9010553","R-HSA-71291","R-HSA-376176","R-HSA-15869", #darkred
              "R-HSA-913531","R-HSA-168164","R-HSA-168188","R-HSA-166166","R-HSA-168181","R-HSA-168638", #blue
              "R-HSA-5693537","R-HSA-5693532","R-HSA-69620", #salmon
              "R-HSA-191273","R-HSA-1989781","R-HSA-8957322","R-HSA-6798695","R-HSA-8951664", #royalblue
              "R-HSA-5617833","R-HSA-1852241","R-HSA-5620912", #magenta
              "R-HSA-5696398","R-HSA-72086") #darkturquoise
rea.cons.filt <- rea.cons.filt %>% filter(ID %in% rea.vect) %>% filter(Cluster != "purple")

ggplot(rea.cons.filt, aes(x=Count, fct_reorder(Description, Count))) + 
       geom_segment(aes(xend=0, yend = Description)) +
       geom_point(aes(color=p.adjust, size = 0.3)) +
       scale_color_distiller(palette = "Spectral", direction = 1) +
       facet_wrap(.~Cluster, scales = "free_y", ncol = 1)+
       theme(text = element_text(size=6), axis.text.x = element_text(size=6),axis.text.y = element_text(size=6))


# Create Heatmpa of interesting enrichment of consensus modules
vectors = list(
  tlr = c("BIRC3","CASP8","MAP2K4","MEF2A","ALPK1","NFKB2","LGMN","NFKBIA","N4BP1","RIPK2","GSDMD","NFKBIB",
          "NOD1","MAP3K8","NFKB1","UBE2D3","UNC93B1","BIRC2","PPP2CA","CREB1","TICAM1",
          "RIPK3","IRAK2","TANK","RIPK1","TLR2","PPP2R1B","NLRC5","UBC","TAB3","TLR3","HSP90B1","FADD","TNIP2",
          "MAP2K1","FOS","MYD88","TRAF6","TBK1","USP18","SOCS1","IRF7","TLR7","IKBKE"),
  chol = c("MSMO1","IDI1","FDFT1","SC5D","HMGCS1","ACAT2","EBP","FDPS"),
  neddy = c("PSMD7","PSMD14","UCHL3","UBE2D2","PSMB7","UBA3","PSMC2"),
  capp = c("PAPOLA","HNRNPC","NUP50","CDC5L","POLDIP3","SNW1","PNN","CRNKL1","WBP4","HNRNPA2B1","UPF3B","NSRP1","PPIL4","PPIG","MFAP1","EIF4A3","SDE2","SLU7","CPSF2","ZRSR2","DDX23","SNRNP35","GTF2F2"),
  ppar = c("ME1","FDFT1","CCNC","HMGCS1","PPARG")
)

# Apply function to each vector
processed.data.heat = lapply(vectors, process.vec.heat, data.rna, rna.meta$condition)

# Create heatmaps
heatmaps = lapply(processed.data.heat, create.heatmap, col_fun=col.fun)
draw(heatmaps$tlr %v% heatmaps$chol %v% heatmaps$ppar %v% heatmaps$neddy %v% heatmaps$capp)


#### Differential Connectivity Analysis ####

# Calculate soft connectivity
kp=SoftConnectivity(data.rna.preterm, 9)
kt=SoftConnectivity(data.rna.term, 12)
ka=SoftConnectivity(data.rna.adult,14)
Kp=kp/max(kp)
Kt=kt/max(kt)
Ka=ka/max(ka)

# Calculate connectivity between the networks
diff.con.avsp <- as.data.frame(Ka-Kp)
colnames(diff.con.avsp) <- "diff"
diff.con.avsp$group <- rep("Group", dim(diff.con.avsp)[1])
summary(diff.con.avsp$diff > (0.4))

diff.con.avst <- as.data.frame(Ka-Kt)
colnames(diff.con.avst) <- "diff"
diff.con.avst$group <- rep("Group", dim(diff.con.avst)[1])
summary(diff.con.avst$diff < (-0.4))


diff.con.pvst <- as.data.frame(Kp-Kt)
colnames(diff.con.pvst) <- "diff"
diff.con.pvst$group <- rep("Group", dim(diff.con.pvst)[1])
summary(diff.con.pvst$diff > (0.4))


# Plot connectivity as seen in Figure 5A
ggplot(data = diff.con.avsp, aes(x = diff, y = group)) +
       geom_point(aes(color = ifelse(diff > 0.4, "Above 0.4", ifelse(diff < -0.4, "Below -0.4", "Inside [-0.4, 0.4]"))),
       size = 0.01, position = position_jitter(width = 0.0000001, height = 0.1)) +
       scale_x_continuous(breaks = c(-1, -0.4, 0, 0.4, 1), limits = c(-1, 1)) +
       scale_color_manual(values = c("Above 0.4" = "#5E4FA2", "Below -0.4" = "#D9464D", "Inside [-0.4, 0.4]" = "black")) +
       geom_vline(xintercept = c(-0.4, 0.4)) +
       theme_minimal() +
       theme(legend.position = "none")

ggplot(data = diff.con.avst, aes(x = diff, y = group)) +
       geom_point(aes(color = ifelse(diff > 0.4, "Above 0.4", ifelse(diff < -0.4, "Below -0.4", "Inside [-0.4, 0.4]"))),
       size = 0.01, position = position_jitter(width = 0.0000001, height = 0.1)) +
       scale_x_continuous(breaks = c(-1, -0.4, 0, 0.4, 1), limits = c(-1, 1)) +
       scale_color_manual(values = c("Above 0.4" = "#5E4FA2", "Below -0.4" = "#F7AB61", "Inside [-0.4, 0.4]" = "black")) +
       geom_vline(xintercept = c(-0.4, 0.4)) +
       theme_minimal() +
       theme(legend.position = "none")

ggplot(data = diff.con.pvst, aes(x = diff, y = group)) +
       geom_point(aes(color = ifelse(diff > 0.4, "Above 0.4", ifelse(diff < -0.4, "Below -0.4", "Inside [-0.4, 0.4]"))),
       size = 0.01, position = position_jitter(width = 0.0000001, height = 0.1)) +
       scale_x_continuous(breaks = c(-1, -0.4, 0, 0.4, 1), limits = c(-1, 1)) +
       scale_color_manual(values = c("Above 0.4" = "#D9464D", "Below -0.4" = "#F7AB61", "Inside [-0.4, 0.4]" = "black")) +
       geom_vline(xintercept = c(-0.4, 0.4)) +
       theme_minimal() +
       theme(legend.position = "none")


diff.con.preterm <- data.frame(genes = colnames(data.rna.preterm[, diff.con.avsp$diff < -0.4])) %>% 
                    cbind(entrez = mapIds(org.Hs.eg.db, .$genes, 'ENTREZID', 'SYMBOL')) %>% na.omit()
enrich.diff.con.preterm <- enrichGO(diff.con.preterm$entrez, OrgDb='org.Hs.eg.db', ont="BP")

diff.con.adult <- data.frame(genes = colnames(data.rna.adult[, diff.con.avsp$diff > 0.4])) %>% 
                  cbind(entrez = mapIds(org.Hs.eg.db, .$genes, 'ENTREZID', 'SYMBOL')) %>% na.omit()
enrich.diff.con.adult <- enrichGO(diff.con.adult$entrez, OrgDb='org.Hs.eg.db', ont="BP")

diff.con.term <- data.frame(genes = colnames(data.rna.term[, diff.con.avst$diff < -0.4])) %>% 
  cbind(entrez = mapIds(org.Hs.eg.db, .$genes, 'ENTREZID', 'SYMBOL')) %>% na.omit()
enrich.diff.con.term <- enrichGO(diff.con.term$entrez, OrgDb='org.Hs.eg.db', ont="BP")

diff.con.adult2 <- data.frame(genes = colnames(data.rna.adult[, diff.con.avst$diff > 0.4])) %>% 
  cbind(entrez = mapIds(org.Hs.eg.db, .$genes, 'ENTREZID', 'SYMBOL')) %>% na.omit()
enrich.diff.con.adult2 <- enrichGO(diff.con.adult2$entrez, OrgDb='org.Hs.eg.db', ont="BP")

write.table(diff.con.preterm, "diff.con.preterm.txt", quote = FALSE, sep="\t")

preterm.go.vec <- c("GO:0051092", "GO:0050900", "GO:0032757", "GO:0032722", "GO:0032612", "GO:0030101")
ggplot(filter(enrich.diff.con.preterm, ID %in% preterm.go.vec), showCategory = 15, 
       aes(Count, fct_reorder(Description, Count))) + 
       geom_segment(aes(xend=0, yend = Description)) +
       geom_point(aes(color="#D9464D", size = 0.3)) +
       scale_y_discrete(position = "left", labels = function(y) str_wrap(y, width = 40))+
       theme_minimal() + 
       xlab("Count") +
       ylab(NULL) +
       theme(text = element_text(size=7), axis.text.x = element_text(size=7),
             axis.text.y = element_text(size=7),legend.position = "none")

term.go.vec <- c("GO:0051607","GO:1903900","GO:0140888","GO:0002221","GO:0046596","GO:0033108")
ggplot(filter(enrich.diff.con.term, ID %in% term.go.vec), showCategory = 15, 
       aes(Count, fct_reorder(Description, Count))) + 
       geom_segment(aes(xend=0, yend = Description)) +
       geom_point(aes(color="#F7AB61", size = 0.3)) +
       scale_y_discrete(position = "left", labels = function(y) str_wrap(y, width = 40))+
       theme_minimal() + 
       xlab("Count") +
       ylab(NULL) +
       theme(text = element_text(size=7), axis.text.x = element_text(size=7),
             axis.text.y = element_text(size=7),legend.position = "none")

ggplot(enrich.diff.con.adult2, showCategory = 6, 
       aes(Count, fct_reorder(Description, Count))) + 
       geom_segment(aes(xend=0, yend = Description)) +
       geom_point(aes(color="#5E4FA2", size = 0.3)) +
       scale_y_discrete(position = "left", labels = function(y) str_wrap(y, width = 40))+
       theme_minimal() + 
       xlab("Count") +
       ylab(NULL) +
       theme(text = element_text(size=7), axis.text.x = element_text(size=7),
             axis.text.y = element_text(size=7),legend.position = "none")


