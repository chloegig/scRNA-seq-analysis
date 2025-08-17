# Lung Cancer VISTA single cell RNA seq analysis
# GSE176091 from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176091 
# developed by /Chloe Giglio
# last updated 8/17/2025


#
##
####
#######
########### Preprocessing 
#######
####
##
#


#
##
### Run when starting up
##
#

install.packages(c('tidyverse', 'pheatmap'))
install.packages(c('RColorBrewer','scales', 'cowplot','patchwork','grid','gridExtra','harmony','ggplot2',
                   'dplyr','NMF'))
install.packages(c('ggalluvial', 'data.table'))

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("dittoSeq")
install.packages(c('ggdendro','tidyr', 'reshape'))
install.packages('ggplot2')
BiocManager::install("enrichplot")
BiocManager::install("pathview")
BiocManager::install("LRBaseDbi")
BiocManager::install("AnnotationHub")
BiocManager::install("grid")
BiocManager::install("ComplexHeatmap")
BiocManager::install("BiocNeighbors")
install.packages('devtools')
devtools::install_github("sqjin/CellChat")
BiocManager::install("Seurat")


# Downloading Packages and Setting Working Directory
```{r}
library(tidyverse)
library(Seurat)
library(pheatmap)
library(RColorBrewer)
library(scales)
library(cowplot)
library(patchwork)
library(grid)
library(gridExtra)
library(harmony)
library(ggplot2) 
library("dittoSeq")
library(clusterProfiler)
library ("enrichplot")
library ("pathview")
library("LRBaseDbi")
library("AnnotationHub")
library("ggdendro")
library("grid")
library("tidyr")
library ("reshape")
library(dplyr)
library ("ComplexHeatmap")
library("BiocNeighbors")
library(CellChat)
library(NMF)
library(ggalluvial)
library(data.table)
library (ComplexHeatmap)
library (ggpubr)
options(stringsAsFactors = FALSE)


# set the working directory
setwd("C:/Documents/LC") 
set.seed(1383)


#
##
### Creating Seurat object
##
#

# Create initial Seurat objects
#CD45_Control_1
# set the folder to pull the files from
data_dir <- "CD45_Control_1"
# check that the folder has the correct files
# should output [1] "barcodes.tsv.gz" "features.tsv.gz" "matrix.mtx.gz"
list.files(data_dir)
# convert the three files into a Seurat Object 
CD45_Control_1 <- Read10X(data.dir = data_dir)
CD45_Control_1 <- CreateSeuratObject(counts=CD45_Control_1, project= "CD45_Control_1")

###CD45_Control_2
data_dir <- "CD45_Control_2"
list.files(data_dir)
CD45_Control_2 <- Read10X(data.dir = data_dir)
CD45_Control_2 <- CreateSeuratObject(counts=CD45_Control_2, project= "CD45_Control_2")

###CD45_Adjuvant_1
data_dir <- "CD45_Adjuvant_1"
list.files(data_dir)
CD45_Adjuvant_1 <- Read10X(data.dir = data_dir)
CD45_Adjuvant_1 <- CreateSeuratObject(counts=CD45_Adjuvant_1, project= "CD45_Adjuvant_1")

###CD45_Adjuvant_2
data_dir <- "CD45_Adjuvant_2"
list.files(data_dir)
CD45_Adjuvant_2 <- Read10X(data.dir = data_dir)
CD45_Adjuvant_2 <- CreateSeuratObject(counts=CD45_Adjuvant_2, project= "CD45_Adjuvant_2")

###CD45_CA170_1
data_dir <- "CD45_CA170_1"
list.files(data_dir)
CD45_CA170_1 <- Read10X(data.dir = data_dir)
CD45_CA170_1 <- CreateSeuratObject(counts=CD45_CA170_1, project= "CD45_CA170_1")

###CD45_CA170_2
data_dir <- "CD45_CA170_2"
list.files(data_dir)
CD45_CA170_2 <- Read10X(data.dir = data_dir)
CD45_CA170_2 <- CreateSeuratObject(counts=CD45_CA170_2, project= "CD45_CA170_2")

###CD45_Kvax_1
data_dir <- "CD45_Kvax_1"
list.files(data_dir)
CD45_Kvax_1 <- Read10X(data.dir = data_dir)
CD45_Kvax_1 <- CreateSeuratObject(counts=CD45_Kvax_1, project= "CD45_Kvax_1")

###CD45_Kvax_2
data_dir <- "CD45_Kvax_2"
list.files(data_dir)
CD45_Kvax_2 <- Read10X(data.dir = data_dir)
CD45_Kvax_2 <- CreateSeuratObject(counts=CD45_Kvax_2, project= "CD45_Kvax_2")

###CD45_CA170_KVax_1
data_dir <- "CD45_CA170_KVax_1"
list.files(data_dir)
CD45_CA170_KVax_1 <- Read10X(data.dir = data_dir)
CD45_CA170_KVax_1 <- CreateSeuratObject(counts=CD45_CA170_KVax_1, project= "CD45_CA170_KVax_1")

###CD45_CA170_KVax_2
data_dir <- "CD45_CA170_KVax_2"
list.files(data_dir)
CD45_CA170_KVax_2 <- Read10X(data.dir = data_dir)
CD45_CA170_KVax_2 <- CreateSeuratObject(counts=CD45_CA170_KVax_2, project= "CD45_CA170_KVax_2")


# merge the different treatment objects into one object to analyze them together
# add.cell.ids labels each individual condition
lung_CD45 <- merge(x= CD45_Control_1, y=c(CD45_Control_2, CD45_Adjuvant_1, CD45_Adjuvant_2, 
                              CD45_CA170_1, CD45_CA170_2, CD45_Kvax_1, CD45_Kvax_2, 
                              CD45_CA170_KVax_1, CD45_CA170_KVax_2), add.cell.ids=c("CD45_Control_1",
                                           "CD45_Control_2", "CD45_Adjuvant_1", "CD45_Adjuvant_2", 
                                           "CD45_CA170_1", "CD45_CA170_2", "CD45_Kvax_1", 
                                           "CD45_Kvax_2", "CD45_CA170_KVax_1", "CD45_CA170_KVax_2"), project = "lung_CD45")

#
##
### Quality Control (QC)
##
#

# add the mitochondrial DNA (mt) stats
lung_CD45[["percent.mt"]] <- PercentageFeatureSet(object = lung_CD45, pattern = "^mt-")

# plot the number of genes detected in each cell (nFeatures), the total number of molecules
# detected within a cell (nCounts), and percent.mt
# Low nFeature_RNA for a cell indicates that it may be dead/dying or an empty droplet. 
# High nCount_RNA and/or nFeature_RNA indicates that the "cell" may in fact be a doublet (or a multiplet).
VlnPlot(lung_CD45, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3, pt.size=0)

# Filter cells based on the nFeatures and mitochondrial content by subsetting out the cells that 
# go beyond the threshold. These numbers will vary dataset to dataset.
# The typical percent.mt cut off is 10%. 
lung_CD45 <- subset(lung_CD45, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
# revisualize the QC metrics
VlnPlot(lung_CD45, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3, pt.size=0)

# normalize the gene expression
lung_CD45 <- NormalizeData(lung_CD45)
# find the variable features in the object
lung_CD45 <- FindVariableFeatures(lung_CD45, selection.method = "vst", nfeatures=2000)
# produce the variable feature plot where the most extreme points are the most variable
VariableFeaturePlot(lung_CD45)

# name rows
all_genes45 <- rownames(lung_CD45)

# set the global max size to accomodate all the memory
options(future.globals.maxSize = 1 * 1024^6)

# scale data. Can take a few minutes to run
lung_CD45 <- ScaleData(lung_CD45, features = all_genes45)

#
##
### Run Principle Component Analysis (PCA) and Uniform Manifold Approximation and Projection (UMAP)
##
#

# Run PCA
lung_CD45 <- RunPCA(lung_CD45, features = VariableFeatures(object = lung_CD45))
# this plots the PCA according to principle component 1 and principle component 2
DimPlot(lung_CD45, reduction = "pca")
# Visualization of PCA between treatment groups. In particular, look for whether any samples are very underrepresented
DimPlot(lung_CD45, reduction = "pca", split.by = "orig.ident")
# choose a number of principal components to move forward with 
ElbowPlot(lung_CD45)

# Cluster cells for UMAP
# change dims to reflect the number of principal components. In this case, I used 15 PCs. 
lung_CD45 <- FindNeighbors(lung_CD45, dims = 1:15)
# change the resolution based on clustering strategy
# For normal clustering, use resolution = 0.5
# For over clustering, use resolution = 1.5
lung_CD45 <- FindClusters(lung_CD45, resolution = 1.5)

# Run UMAP
lung_CD45 <- RunUMAP(lung_CD45, dims = 1:15)
# Visualize UMAP
DimPlot(lung_CD45, label = FALSE, reduction="umap", split.by = "orig.ident")


#
##
### Save and read Seurat object
##
#


# This point marks all the pre-processing for this analysis. 
# Save the RDS file to access this seurat object for analysis. 
saveRDS(lung_CD45, file = "scrna_seq_LC45.rds")

# Read file
lung_CD45 <- readRDS("scrna_seq_LC45.rds")


# combine idents (treatment groups) of lung_CD45. This combines CD45_Control_1 and CD45_Control_2
# into one ident called Control, etc.
lung_CD45 <- SetIdent(lung_CD45, value = 'orig.ident')
lung_CD45<- RenameIdents(lung_CD45, 'CD45_Control_1' = 'Control',
                         'CD45_Control_2' = 'Control', 'CD45_Adjuvant_1' = 'Adjuvant', 
                         'CD45_Adjuvant_2' = 'Adjuvant', 'CD45_CA170_1' = 'CA170',
                         'CD45_CA170_2' = 'CA170',  'CD45_Kvax_1' = 'Kvax',
                         'CD45_Kvax_2' = 'Kvax', 'CD45_CA170_KVax_1' = 'Combo', 
                         'CD45_CA170_KVax_2' = 'Combo')
lung_CD45<- RenameIdents(lung_CD45, 'Control' = 'Control',
                         'Adjuvant' = 'Adjuvant', 
                         'CA170' = 'CA170',
                         'Kvax' = 'Kvax',
                         'CA170_Kvax' = 'Combo')
# Stash cell identity genotypes)
lung_CD45[['orig.ident']] <- Idents(lung_CD45)
#set ident to seurat clusters
lung_CD45 <- SetIdent(lung_CD45, value = 'seurat_clusters')


#
##
####
#######
########### Sub clustering and identifying cell types in the CD45+ object
#######
####
##
#

# Find the markers that define each cluster and save data in a csv file
lung_CD45_markers <- FindAllMarkers(lung_CD45, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# arrange the data nicely
lung_CD45_markers %>% group_by(cluster) %>% top_n(n=25, wt = avg_log2FC)
# save the markers in a comma-separated values table that can be accessed in excel
write.csv(lung_CD45_markers, file = "lung_CD45_markers.csv")


# dotplot for CD45 cluster identification based on key gene markers
plot<-DotPlot (object = lung_CD45, features=c("Cd3d","Cd3g", "S100a9", "Bank1", "Cd79a", "Cd79b", 
                                              "Klrd1", "Lyz2", "Chil3", "Lpl", "Wfdc17", 
                                              "Siglech"), dot.scale = 7, cols = c("lightgrey", "blue"), scale=T)
plot + theme(axis.text.x = element_text(angle = 45, hjust=1))


# rename clusters of CD45 according to cell type
# list out the cluster names in order and with the same spelling 
# IMPORTANT: the length of this list MUST match the number of total clusters 
#lung_CD45 <- SetIdent(lung_CD45, value = 'RNA_snn_res.1.5')
new.cluster.ids.CD45id <- c("B cell", "B cell", "T cell", "B cell", "T cell", "Neutrophil", "T cell", "T cell",
                            "Myeloid", "NK cell", "T cell", "Myeloid", "T cell", "Neutrophil", "Myeloid",
                            "Myeloid", "T cell", "NKT cell", "Myeloid", "Myeloid", "B cell", "T cell", 
                            "T cell", "B cell", "T cell", "Myeloid", "Neutrophil", "Myeloid", "Myeloid",
                            "Dendritic cell", "Neutrophil", "T cell", "T cell", "Myeloid", "Myeloid", "T cell",
                            "Neutrophil", "Myeloid", "B cell")
names(new.cluster.ids.CD45id) <- levels(lung_CD45)
lung_CD45  <- RenameIdents(lung_CD45, new.cluster.ids.CD45id)
#lung_CD45[['RNA_snn_res.1.5']] <- Idents(lung_CD45)
#lung_CD45 <- SetIdent(lung_CD45, value = 'cell_type')

colors <- c("B cell" = "coral2", "Dendritic cell" = "magenta", "Myeloid" ="springgreen4", 
            "Neutrophil" = "olivedrab3", "NK cell" ="cyan3", "NKT cell" = "mediumpurple1", "T cell" = "darkgoldenrod3")

# visualize cell clusters separated by treatment group
DimPlot(lung_CD45, label = FALSE, split.by= "orig.ident", cols = colors)

#
##
### proportions 
##
#
#overall
# combine cell types into separate objects
Bcell <- subset(lung_CD45, ident = "B cell")
Tcell <- subset(lung_CD45, ident = c("T cell"))
Neutrophil <- subset(lung_CD45, ident = c("Neutrophil"))
Myeloid <- subset (lung_CD45, ident = c("Myeloid"))
NKT <- subset (lung_CD45, ident = "NKT cell")
Dendritic <- subset (lung_CD45, ident = "Dendritic cell")
colors <- c("B cell" = "coral2", "Dendritic cell" = "magenta", "Myeloid" ="springgreen4", 
            "Neutrophil" = "olivedrab3", "NK cell" ="cyan3", "NKT cell" = "mediumpurple1", "T cell" = "darkgoldenrod2")
dittoBarPlot(object = lung_CD45, var = "RNA_snn_res.1.5",var.labels.rename = new.cluster.ids.CD45id,
             group.by = "orig.ident", data.out = TRUE, x.reorder = c(4,1,2,5,3), color.panel = colors)
metadata <- lung_CD45@meta.data
cell_counts <- table(metadata$RNA_snn_res.1.5, metadata$orig.ident)
print(cell_counts)


#B cells
Bcell_subcluster <- readRDS("Bcell_subcluster.rds")
Bcell_subcluster <- SetIdent(Bcell_subcluster, value = 'orig.ident')
#Bcell_subcluster<- RenameIdents(Bcell_subcluster, 'Control' ='Control', 'Adjuvant' = 'Adjuvant',
                                'CA170' = 'CA170', 'Kvax' = 'Kvax', 'CA170_KVax' = 'Combo')
# Stash cell identity genotypes)
Bcell_subcluster[['orig.ident']] <- Idents(Bcell_subcluster)
#set ident to seurat clusters
Bcell_subcluster <- SetIdent(Bcell_subcluster, value = 'seurat_clusters')
new.cluster.id.Bcell_subclusterid <- c("memory", "memory", "follicular", "follicular", "follicular", "follicular", "plasma cell", "memory", "plasma cell")
names(new.cluster.id.Bcell_subclusterid) <- levels(Bcell_subcluster) 
Bcell_subcluster  <- RenameIdents(Bcell_subcluster, new.cluster.id.Bcell_subclusterid)

#organize orig.ident by subtype instead of treatment group
#Bcell_subcluster[['RNA_snn_res.0.5']] <- Idents(Bcell_subcluster)
#Bcell_subcluster <- SetIdent(Bcell_subcluster, value = 'cell_type')
# combine cell types into separate objects
follicular <- subset(Bcell_subcluster, ident = "follicular")
memory <- subset(Bcell_subcluster, ident = "memory")
plasma <- subset(Bcell_subcluster, ident = "plasma cell")
# assign colors to cell types 
cols <- c("memory" = "orange", "follicular" = "deepskyblue2", "plasma cell" = "forestgreen")
# visualize B cell types separated by treatment group
DimPlot(Bcell_subcluster, label.size = 4, raster = FALSE, split.by="orig.ident", cols = cols)
dittoBarPlot(object = Bcell_subcluster, group.by = "orig.ident", var = "RNA_snn_res.1.5",
              data.out = TRUE)
metadata <- Bcell_subcluster@meta.data
cell_counts <- table(metadata$RNA_snn_res.0.5, metadata$orig.ident)
print(cell_counts)


#T cells
Tcell_subcluster <- readRDS("Tcell_subcluster.rds")
Tcell_subcluster <- SetIdent(Tcell_subcluster, value = 'orig.ident')
Tcell_subcluster<- RenameIdents(Tcell_subcluster, 'CD45_Control_1' = 'Control',
                         'CD45_Control_2' = 'Control', 'CD45_Adjuvant_1' = 'Adjuvant', 
                         'CD45_Adjuvant_2' = 'Adjuvant', 'CD45_CA170_1' = 'CA170',
                         'CD45_CA170_2' = 'CA170',  'CD45_Kvax_1' = 'Kvax',
                         'CD45_Kvax_2' = 'Kvax', 'CD45_CA170_KVax_1' = 'CA170_Kvax',
                         'CD45_CA170_KVax_2' = 'CA170_Kvax')
# Stash cell identity genotypes)
Tcell_subcluster[['orig.ident']] <- Idents(Tcell_subcluster)
#set ident to seurat clusters
Tcell_subcluster <- SetIdent(Tcell_subcluster, value = 'seurat_clusters')
# renaming T cell clusters 
new.cluster.ids.CD45Tid <- c("Helper T", "CD8+/CD4+ Naive", "Helper T", "CD8+/CD4+ Naive", 
                             "Helper T", "CD8+ Naive", "CD8+/CD4+ Mem eff", "CD8+ CTL", "CD4+ CTL")
names(new.cluster.ids.CD45Tid) <- levels(Tcell_subcluster)
Tcell_subcluster  <- RenameIdents(Tcell_subcluster, new.cluster.ids.CD45Tid)
#Tcell_subcluster[['RNA_snn_res.0.5']] <- Idents(Tcell_subcluster)
#Tcell_subcluster <- SetIdent(Tcell_subcluster, value = 'cell_type')
# assign colors to each cell type
cols <- c("Helper T" = "orange", "CD8+/CD4+ Naive" = "skyblue", "CD8+ Naive" = "seagreen", "CD8+/CD4+ Mem eff" = "yellow", "CD8+ CTL" = "dodgerblue4", "CD4+ CTL" = "orangered2")
# visualize B cell types separated by treatment group
DimPlot(object = Tcell_subcluster,  label.size = 4, raster = FALSE,  split.by = "orig.ident", cols = cols)
# combine cell types into separate objects
HelperT <- subset(Tcell_subcluster, ident = "Helper T")
NaiveT <- subset(Tcell_subcluster, ident = c("CD8+/CD4+ Naive", "CD8+ Naive"))
Effector<- subset(Tcell_subcluster, ident = c("CD8+/CD4+ Mem eff", "CD8+ CTL", "CD4+ CTL"))
CTL <- subset (Tcell_subcluster, ident = c("CD8+ CTL", "CD4+ CTL"))
CD8CTL <- subset (Tcell_subcluster, ident = "CD8+ CTL")
CD4CTL <- subset (Tcell_subcluster, ident = "CD4+ CTL")
# bargraph depicting the percent composition of each B cell type
dittoBarPlot(object = Tcell_subcluster, var = "RNA_snn_res.0.5", group.by = "orig.ident", 
             var.labels.rename = new.cluster.ids.CD45Tid, data.out = TRUE, x.reorder = c(4,1,2,5,3))
metadata <- Tcell_subcluster@meta.data
cell_counts <- table(metadata$RNA_snn_res.0.5, metadata$orig.ident)
print(cell_counts)

#
##
####
#######
########### CellChat analysis on CD45+ object
#######
####
##
#

# normalize data matrix
data.input <- lung_CD45[["RNA"]]$data
# add labels for the treatment groups
labels <- Idents(lung_CD45)
# create a dataframe of the cell labels
meta <- data.frame(labels = labels, row.names = names(labels)) 
#create CellChat object
cellChat_full <- createCellChat(object = lung_CD45, group.by = "ident", assay = "RNA") 

# add metadata
cellChat_full <- addMeta(cellChat_full, meta = meta)
# set "labels" as default cell identity
cellChat_full <- setIdent(cellChat_full, ident.use = "labels")
# show factor levels of the cell labels
levels(cellChat_full@idents) 
# number of cells in each cell group
groupSize <- as.numeric(table(cellChat_full@idents)) 

#subsetting out the genes that are not contained in the CD45+ object  
gene_data <- cellChat_full@data
genes_to_remove <- c("H2-BI", "H2-Ea-ps")
subset_gene_data <- gene_data[!rownames(gene_data) %in% genes_to_remove, ]
cellChat_full@data <- subset_gene_data

# Set the ligand-receptor interaction database to the mouse database
CellChatDB <- CellChatDB.mouse 
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB)
# Ensure other related slots are also updated if necessary. For example, the CellChatDB interaction data
filtered_interactions <- CellChatDB$interaction$ligand.receptor[!(CellChatDB$interaction$ligand.receptor[, "ligand"] %in% genes_to_remove | CellChatDB$interaction$ligand.receptor[, "receptor"] %in% genes_to_remove), ]
CellChatDB$interaction$ligand.receptor <- filtered_interactions
cellChat_full@DB <- CellChatDB

#Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cellChat_full <- subsetData(cellChat_full)
future::plan("multisession", workers = 4) # do parallel

# set the global max size to accomodate all the memory
options(future.globals.maxSize = 1 * 1024^6)
cellChat_full <- identifyOverExpressedGenes(cellChat_full)
cellChat_full <- identifyOverExpressedInteractions(cellChat_full)

#Compute the communication probability and infer cellular communication network
cellChat_full <- computeCommunProb(cellChat_full, type = "triMean")
cellChat_full <- filterCommunication(cellChat_full, min.cells = 10)

#Infer the cell-cell communication at a signaling pathway level
cellChat_full <- computeCommunProbPathway(cellChat_full)

#Calculate the aggregated cell-cell communication network
cellChat_full <- aggregateNet(cellChat_full)
groupSize <- as.numeric(table(cellChat_full@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellChat_full@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# View the available pathways
significant_pathways <- cellChat_full@netP$pathways
print(significant_pathways)

# choose the most highly expressed pathways
pathways.show <- c("MHC-I","SELPLG","CCL","CD52" ,"CD45","MHC-II" ,"MIF","CD22","SPP1" ,"APP", "BST2") 

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway 
netAnalysis_contribution(cellChat_full, signaling = pathways.show)

# Save the CellChat object
saveRDS(cellChat_fullT, file = "cellchat_fullT.rds")


# heatmaps of the top receptors and ligands
# top 20 receptors
genes_of_interest <- c("Sell", "Cd8b1", "Siglecg", "Cd8a", "Cd8a", "Ccr1", "Mrc1", "Cd22", "Ptprc", "Cd74",
                       "Ccr1", "Pira2", "Ccr5", "Klrd1", "Klrc2", "Ccr2", "Cxcr4", "Cd44", "Cd4", "Itga4")

# top 20 ligands
genes_of_interest <-c("Selplg", "H2-D1", "H2-K1", "Cd52", "Ccl5", "Ptprc", "Cd22", "H2-T23", "H2-Q7", 
                      "App", "Ccl6", "Bst2", "Mif", "Spp1", "H2-Aa", "H2-Ab1", "H2-Eb1", "H2-T22", 
                      "H2-DMa", "H2-Q6")
genes_of_interest <- c("H2-D1", "H2-K1", "H2-T23", "H2-Q7","H2-Aa", "H2-Ab1", "H2-Eb1", "H2-T22", 
                       "H2-DMa", "H2-Q6")

expression_data <- FetchData(lung_CD45, vars = genes_of_interest)

#get metadata for treatment groups from object
treatment_groups <- unique(lung_CD45@meta.data$orig.ident)
print(treatment_groups)
# Initialize a matrix to store the mean expression values
expression_by_group <- matrix(nrow = length(genes_of_interest), ncol = length(treatment_groups))
rownames(expression_by_group) <- genes_of_interest
colnames(expression_by_group) <- treatment_groups

# Calculate the mean expression for each gene in each treatment group
for (group in treatment_groups) {
  cells_in_group <- WhichCells(lung_CD45, expression = orig.ident == group)
  expression_data_group <- GetAssayData(lung_CD45, slot = "data")[genes_of_interest, cells_in_group, drop=FALSE]
  expression_by_group[, group] <- rowMeans(expression_data_group)
}

# Create the heatmap
pheatmap(expression_by_group, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         scale = "row",
         color = colorRampPalette(c("blue", "white", "red"))(50), name = "Expression", 
         angle_col = "0", 
         fontsize_col = 12)


#
##
####
#######
########### CD45+ B cell analysis
#######
####
##
#

#
##
### sub clustering and identifying cell types in B cell object
##
#

# sub setting CD45 object for B cells 
Bcell_subcluster <- subset(lung_CD45, idents="B cell")
# sub clustering B cell object for more specific cell types
Bcell_subcluster <- FindNeighbors(Bcell_subcluster, dims = 1:6)
Bcell_subcluster <- FindClusters(Bcell_subcluster, resolution = 0.5)
# make UMAP B cells
Bcell_subcluster <- RunUMAP(Bcell_subcluster, dims = 1:6)
DimPlot(Bcell_subcluster, label = TRUE, reduction="umap")
# save sub cluster
saveRDS(Bcell_subcluster, file = "Bcell_subcluster.rds")
Bcell_subcluster <- readRDS("Bcell_subcluster.rds")
# find markers
CD45_Bcell_markers <- FindAllMarkers(Bcell_subcluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# arranges the data nicely
CD45_Bcell_markers %>% group_by(cluster) %>% top_n(n=25, wt = avg_log2FC)
# saves the markers in a comma-separated values table that can be accessed in excel
write.csv(CD45_Bcell_markers, file = "CD45_Bcell_markers.csv")
#save Bcell45_markers
saveRDS(CD45_Bcell_markers, file = "CD45_Bcell_markers.rds")
CD45_Bcell_markers <- readRDS("CD45_Bcell_markers.rds")

# identifying sub clusters with dotplot and key genes
plot<-DotPlot (object = Bcell_subcluster, dot.scale = 7, cols = c("lightgrey", "blue"), scale=T, 
               features=c("Ms4a1", "Cd19", "Sell", "Cd79a", "Cd79b", "Cd40", "Cr2", "Lmo2", "Bank1", 
                          "Il4ra", "Fcer2a", "Bach2","Mzb1", "Tnfrsf17" , "Sdc1", "Jchain", "Ebf1", 
                          "Cd24a", "Eif4ebp1","Tnfrsf13c", "Ltb", "Ctsh", "Ly6a", "Bcl2", "Fcrl1"))
plot + theme(axis.text.x = element_text(angle = 45, hjust=1))

# renaming B cell clusters
new.cluster.id.Bcell_subclusterid <- c("memory", "memory", "follicular", "follicular", "follicular", "follicular", "plasma cell", "memory", "plasma cell")
names(new.cluster.id.Bcell_subclusterid) <- levels(Bcell_subcluster) 
Bcell_subcluster  <- RenameIdents(Bcell_subcluster, new.cluster.id.Bcell_subclusterid)
# assign colors to cell types 
cols <- c("memory" = "orange", "follicular" = "skyblue", "plasma cell" = "seagreen")
# visualize B cell types separated by treatment group
DimPlot(Bcell_subcluster, label.size = 4, raster = FALSE, split.by="orig.ident", cols = cols)

# combine cell types into separate objects
follicular <- subset(Bcell_subcluster, ident = "follicular")
memory <- subset(Bcell_subcluster, ident = "memory")
plasma <- subset(Bcell_subcluster, ident = "plasma cell")

# save Bcell_subcluster
saveRDS(Bcell_subcluster, file = "Bcell_subcluster.rds")
Bcell_subcluster <- readRDS("Bcell_subcluster.rds")
Bcell_subcluster <- SetIdent(Bcell_subcluster, value = 'orig.ident')
Bcell_subcluster<- RenameIdents(Bcell_subcluster, 'Control' ='Control', 'Adjuvant' = 'Adjuvant',
                                'CA170' = 'CA170', 'Kvax' = 'Kvax', 'CA170_KVax' = 'Combo')
Bcell_subcluster[['orig.ident']] <- Idents(Bcell_subcluster)

#
##
### Heatmap for cell subtype markers
##
#

#memory markers
genes_of_interest <- c("Cd80","Cd19" ,"Cd27", "Cd38", "Cd40", "Pax5", "Spib")
#plasma markers
genes_of_interest <- c("Cd19", "Cd27", "Cd38", "Cd93", "Cxcr4", "Ly6k", "Irf4", "Xbp1")
#follicular markers
genes_of_interest<- c("Cd19", "Cd22", "Cd38", "Cd93", "Cxcr5")

expression_data <- FetchData(Bcell_subcluster, vars = genes_of_interest)

#get metadata for treatment groups from object
treatment_groups <- unique(Bcell_subcluster@meta.data$orig.ident)
print(treatment_groups)
# Initialize a matrix to store the mean expression values
expression_by_group <- matrix(nrow = length(genes_of_interest), ncol = length(treatment_groups))
rownames(expression_by_group) <- genes_of_interest
colnames(expression_by_group) <- treatment_groups

# Calculate the mean expression for each gene in each treatment group
for (group in treatment_groups) {
  cells_in_group <- WhichCells(Bcell_subcluster, expression = orig.ident == group)
  expression_data_group <- GetAssayData(Bcell_subcluster, slot = "data")[genes_of_interest, cells_in_group, drop=FALSE]
  expression_by_group[, group] <- rowMeans(expression_data_group)
}
expression_by_group <- as.matrix(expression_by_group)

# Create the heatmap
pheatmap(expression_by_group, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         scale = "row",
         color = colorRampPalette(c("blue", "white", "red"))(50), name = "Expression",
         angle_col = "0", fontsize_col = 12)

### troubleshooting if index out of range 
treatment_groups <- unique(Bcell_subcluster@meta.data$orig.ident)
print(treatment_groups)

expression_data <- GetAssayData(Bcell_subcluster, slot = "data")
missing_genes <- setdiff(genes_of_interest, rownames(expression_data))
if (length(missing_genes) > 0) {
  stop("The following genes are missing from the expression data: ", paste(missing_genes, collapse = ", "))
}


#
##
### Pathway analysis
##
#


# bargraph depicting the percent composition of each B cell type
dittoBarPlot(object = Bcell_subcluster, var = "RNA_snn_res.0.5", group.by = "orig.ident", 
             var.labels.rename = new.cluster.id.Bcell_subclusterid, data.out = TRUE, x.reorder = c(4,1,2,5,3))

# Most gene signatures were obtained from the Gene Ontology (GO) database. Refer to my paper references for
# the past literature from which I obtained the gene signatures not found in the GO database. 

## B cell proliferation
# create a list of genes associated with B cell proliferation
genes <- c("Abl1", "Ada", "Ahr", "Atad5", "Bax", "Bcl2", "Bcl6", "Bmi1", "Bst1", "Card11", "Cd19", 
           "Cd22", "Cd27", "Cd38", "Cd40", "Cd40lg", "Cd70", "Cd74", "Cd79a", "Cd81", "Cd180", "Cd320",
           "Cdkn1a", "Cfb", "Chrnb2", "Clcf1", "Cr2", "Ephb2", "Fosl2", "Gapt", "Gm13271", "Gm13272", 
           "Gm13275", "Gm13276", "Gm13277", "Gm13283", "Gpr183", "Hspd1", "Ifna1", "Ifna2", "Ifna4", "Ifna5", 
           "Ifna6", "Ifna7", "Ifna9", "Ifna11", "Ifna12", "Ifna13", "Ifna14", "Ifna15", "Ifna16", "Ifnab", 
           "Ifnb1", "Ifne", "Ifnk", "Ifnz", "Ikzf3", "Il2", "Il3", "Il4", "Il5",
           "Il7", "Il7r", "Il9", "Il9r", "Il13", "Il21", "Irs2", "Lef1", "Mef2c", "Mif", "Mzb1", "Nckap1l",
           "Nfatc2", "Nfkbiz", "Ntn1", "Plcl2", "Prkcd", "Prlr", "Ptprc", "Rag2", "Rasgrp1", "Sash3", "Shb", 
           "Siglecg", "Slc39a10", "Tcf3", "Tfrc", "Ticam1", "Tirap", "Tlr4", "Tlr9", "Tnfrsf4", "Tnfrsf13c", 
           "Tnfsf13b", "Vav3","Wnt3a")
# calculate average expression of each gene and average it
memory <- AddModuleScore(memory, features=list(genes), name= "Bcellproliferation")
# assign consistent colors associated with each treatment group
cols <- c("Control" = "lightblue", "Adjuvant" = "darkorange", "CA170" = "forestgreen","Kvax"= "purple","Combo" = "gold2")
cols <- c("memory" = "orange", "follicular" = "skyblue", "plasma cell" = "seagreen")
# visualize the average gene expression in each treatment group
VlnPlot(memory,features= "Bcellproliferation1", pt.size = 0, group.by = "orig.ident", col=cols, y.max=0.5) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95)+
  geom_boxplot(width=0.1, fill="white") +
  geom_signif(
    comparisons = list(c("Control", "CA170")),
    annotations = "****",
    textsize = 5,
    map_signif_level = TRUE,
    y_position = 0.35,
    tip_length=0.02
  ) +
  geom_signif(
    comparisons = list(c("Control", "Combo")),
    annotations = "****",
    textsize = 5,
    map_signif_level = TRUE,
    y_position = 0.4,
    tip_length=0.02
  ) + 
  geom_signif(
    comparisons = list(c("CA170", "Combo")),
    annotations = "****",
    textsize = 5,
    map_signif_level = TRUE,
    y_position = 0.35,
    tip_length=0.02
  ) + 
  geom_signif(
    comparisons = list(c("Kvax", "Combo")),
    annotations = "****",
    textsize = 5,
    map_signif_level = TRUE,
    y_position = 0.45,
    tip_length=0.02
  )
df <- FetchData(Bcell_subcluster, vars = c("Bcellproliferation1", "orig.ident"))
medians <- aggregate(Bcellproliferation1 ~ orig.ident, data = df, FUN = mean)
print(medians)

## Follicular Differentiation
genes <- c("Irf8", "Plcg2", "Spi1", "Mki67", "Bcl6")
Bcell_subcluster <- AddModuleScore(Bcell_subcluster, features=list(genes),name= "FollicularDifferentiation")
cols <- c("Control" = "lightblue", "Adjuvant" = "darkorange", "CA170" = "forestgreen","Kvax"= "purple","Combo" = "gold")
VlnPlot(Bcell_subcluster,features= "FollicularDifferentiation1", pt.size = 0, group.by = "orig.ident", col=cols, y.max = 1.9)+
  stat_summary(fun.y = median, geom='point', size = 13, colour = "black", shape = 95)+
  geom_boxplot(width=0.1, fill="white") + 
  geom_signif(
    comparisons = list(c("Control", "Kvax")),
    annotations = "**",
    textsize = 4,
    map_signif_level = TRUE,
    y_position = 1.37,
    tip_length=0.02
  ) + 
  geom_signif(
    comparisons = list(c("Control", "CA170")),
    annotations = "****",
    textsize = 4.5,
    map_signif_level = TRUE,
    y_position = 1.55,
    tip_length=0.02
  ) + 
  geom_signif(
    comparisons = list(c("Control", "Combo")),
    annotations = "****",
    textsize = 4.5,
    map_signif_level = TRUE,
    y_position = 1.73,
    tip_length=0.02
  ) + 
  geom_signif(
    comparisons = list(c("CA170", "Combo")),
    annotations = "****",
    textsize = 4.5,
    map_signif_level = TRUE,
    y_position = 1.55,
    tip_length=0.02
  ) + 
  geom_signif(
    comparisons = list(c("Kvax", "Combo")),
    annotations = "****",
    textsize = 4.5,
    map_signif_level = TRUE,
    y_position = 1.37,
    tip_length=0.02
  )
df <- FetchData(Bcell_subcluster, vars = c("FollicularDifferentiation1", "orig.ident"))
medians <- aggregate(FollicularDifferentiation1 ~ orig.ident, data = df, FUN = mean)
print(medians)

## Memory Differentiation
genes <- c("St3gal1","Cd46","Pck1","Tsc1", "Cd27","Tnfsf4", "Aicda", "Prdm1", "Bcl6", "Cd27", "Cd38", "Cd80", "Cd86","Cd19", "Pax5", "Aicda", "Bcl6", "Mki67", "Cd40", "Icosl", "Cr2","Cxcr4", "Cxcr5", "Tnfsf13b")
Bcell_subcluster <- AddModuleScore(Bcell_subcluster, features=list(genes),name= "MemoryDifferentiation")
cols <- c("Control" = "lightblue", "Adjuvant" = "darkorange", "CA170" = "forestgreen","Kvax"= "purple","Combo" = "gold")
VlnPlot(Bcell_subcluster,features= "MemoryDifferentiation1", pt.size = 0, group.by = "orig.ident", col=cols, y.max = 0.8) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95)+
  geom_boxplot(width=0.1, fill="white")  + 
  geom_signif(
    comparisons = list(c("Control", "CA170")),
    annotations = "****",
    textsize = 5,
    map_signif_level = TRUE,
    y_position = 0.48,
    tip_length=0.02
  ) + 
  geom_signif(
    comparisons = list(c("Control", "Kvax")),
    annotations = "****",
    textsize = 5,
    map_signif_level = TRUE,
    y_position = 0.58,
    tip_length=0.02
  ) + 
  geom_signif(
    comparisons = list(c("Control", "Combo")),
    annotations = "****",
    textsize = 5,
    map_signif_level = TRUE,
    y_position = 0.68,
    tip_length=0.02
  )

df <- FetchData(Bcell_subcluster, vars = c("MemoryDifferentiation1", "orig.ident"))
medians <- aggregate(MemoryDifferentiation1 ~ orig.ident, data = df, FUN = median)
print(medians)

## Memory proliferation 
genes <- c("Cd80","Cd19" ,"Cd27", "Cd38", "Cd40", "Pax5", "Spib", "St3gal1","Cd46","Pck1","Tsc1", "Cd27","Tnfsf4", "Aicda", "Prdm1", "Bcl6", "Cd27", "Cd38", "Cd80", "Cd86","Cd19", "Pax5", "Aicda", "Bcl6", "Mki67", "Cd40", "Icosl", "Cr2","Cxcr4", "Cxcr5", "Tnfsf13b")
memory <- AddModuleScore(memory, features=list(genes),name= "memMemoryDifferentiation")
cols <- c("Control" = "lightblue", "Adjuvant" = "darkorange", "CA170" = "forestgreen","Kvax"= "purple","Combo" = "gold")
VlnPlot(memory,features= "memMemoryDifferentiation1", pt.size = 0, group.by = "orig.ident", col=cols, y.max=0.55) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95)+
  geom_boxplot(width=0.1, fill="white") + 
  geom_signif(
    comparisons = list(c("Control", "CA170")),
    annotations = "****",
    textsize = 4.5,
    map_signif_level = TRUE,
    y_position = 0.4,
    tip_length=0.02
  )+ 
  geom_signif(
    comparisons = list(c("Control", "Kvax")),
    annotations = "****",
    textsize = 4.5,
    map_signif_level = TRUE,
    y_position = 0.47,
    tip_length=0.02
  )
df <- FetchData(memory, vars = c("memMemoryDifferentiation1", "orig.ident"))
medians <- aggregate(memMemoryDifferentiation1 ~ orig.ident, data = df, FUN = median)
print(medians)

## PlasmaPosDifferentiation 
genes <- c ("Enpp1", "Il2", "Il10", "Itm2a", "Lgals1", "Lgals8", "Nfkbiz", "Nkx2-3", "Xbp1", "Xbp1", "Irf4", "Sdc1","Cd27", "Cd38", "Cd40", "Cd19", "Pax5", "Bcl6", "Aicda", "Prdm1")
Bcell_subcluster <- AddModuleScore(Bcell_subcluster, features=list(genes), name= "plasmaPosDifferentiation")
cols <- c("Control" = "lightblue", "Adjuvant" = "darkorange", "CA170" = "forestgreen","Kvax"= "purple",  "Combo" = "gold")
VlnPlot(Bcell_subcluster,features= "plasmaPosDifferentiation1", pt.size = 0, group.by = "orig.ident", col=cols, y.max = 0.75)+
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95)+
  geom_boxplot(width=0.1, fill="white") + 
  geom_signif(
    comparisons = list(c("Control", "CA170")),
    annotations = "****",
    textsize = 4.5,
    map_signif_level = TRUE,
    y_position = 0.5,
    tip_length=0.02
  )+ 
  geom_signif(
    comparisons = list(c("Control", "Kvax")),
    annotations = "****",
    textsize = 4.5,
    map_signif_level = TRUE,
    y_position = 0.58,
    tip_length=0.02
  )+ 
  geom_signif(
    comparisons = list(c("Control", "Combo")),
    annotations = "****",
    textsize = 4.5,
    map_signif_level = TRUE,
    y_position = 0.66,
    tip_length=0.02
  )+ 
  geom_signif(
    comparisons = list(c("CA170", "Combo")),
    annotations = "****",
    textsize = 4.5,
    map_signif_level = TRUE,
    y_position = 0.5,
    tip_length=0.02
  )+ 
  geom_signif(
    comparisons = list(c("Kvax", "Combo")),
    annotations = "****",
    textsize = 4.5,
    map_signif_level = TRUE,
    y_position = 0.58,
    tip_length=0.02
  )

df <- FetchData(plasma, vars = c("plasmaPosDifferentiation1", "orig.ident"))
medians <- aggregate(plasmaPosDifferentiation1 ~ orig.ident, data = df, FUN = median)
print(medians)


#
##
### Statistical testing
##
#


# Creating the csv files for statistical testing. The csv file contains the gene expression level for
# each gene across all treatment groups.

# proliferation
data <- data.table(x1 = Bcell_subcluster$RNA_snn_res.0.5,
                   x2=Bcell_subcluster$Bcellproliferation1)
write.csv(data, "Bcellsubcluster_Bproliferation1_for_stats.csv")

# follicular differentiation
data <- data.table(x1 = Bcell_subcluster$orig.ident,
                   x2=Bcell_subcluster$FollicularDifferentiation1)
write.csv(data, "Bcellsubcluster_FollicularDifferentiation1_for_stats.csv")

# memory differentiation
data <- data.table(x1 = Bcell_subcluster$orig.ident,
                   x2=Bcell_subcluster$MemoryDifferentiation1)
write.csv(data, "Bcellsubcluster_MemoryDifferentiation1_for_stats.csv")

# memory proliferation
data <- data.table(x1 = memory$orig.ident,
                   x2=memory$memMemoryDifferentiation1)
write.csv(data, "memory_memMemoryDifferentiation1_for_stats.csv")


# plasma differentiation
data <- data.table(x1 = Bcell_subcluster$orig.ident,
                   x2=Bcell_subcluster$PlasmaPosDifferentiation1)
write.csv(data, "Bcellsubcluster_PlasmaPosDifferentiation1_for_stats.csv")


#
##
####
#######
########### CD45+ T cell analysis
#######
####
##
#

#
##
### sub clustering and identifying cell types in T cell object
##
#

Tcell_subcluster <- readRDS("Tcell_subcluster.rds")
Tcell_subcluster <- SetIdent(Tcell_subcluster, value = 'orig.ident')
Tcell_subcluster<- RenameIdents(Tcell_subcluster, 'CD45_Control_1' = 'Control',
                                'CD45_Control_2' = 'Control', 'CD45_Adjuvant_1' = 'Adjuvant', 
                                'CD45_Adjuvant_2' = 'Adjuvant', 'CD45_CA170_1' = 'CA170',
                                'CD45_CA170_2' = 'CA170',  'CD45_Kvax_1' = 'Kvax',
                                'CD45_Kvax_2' = 'Kvax', 'CD45_CA170_KVax_1' = 'CA170_Kvax',
                                'CD45_CA170_KVax_2' = 'CA170_Kvax')
# Stash cell identity genotypes)
Tcell_subcluster[['orig.ident']] <- Idents(Tcell_subcluster)
#set ident to seurat clusters
Tcell_subcluster <- SetIdent(Tcell_subcluster, value = 'seurat_clusters')
# renaming T cell clusters 
new.cluster.ids.CD45Tid <- c("Helper T", "CD8+/CD4+ Naive", "Helper T", "CD8+/CD4+ Naive", 
                             "Helper T", "CD8+ Naive", "CD8+/CD4+ Mem eff", "CD8+ CTL", "CD4+ CTL")
names(new.cluster.ids.CD45Tid) <- levels(Tcell_subcluster)
Tcell_subcluster  <- RenameIdents(Tcell_subcluster, new.cluster.ids.CD45Tid)
#Tcell_subcluster[['RNA_snn_res.0.5']] <- Idents(Tcell_subcluster)
#Tcell_subcluster <- SetIdent(Tcell_subcluster, value = 'cell_type')
# assign colors to each cell type
cols <- c("Helper T" = "orange", "CD8+/CD4+ Naive" = "skyblue", "CD8+ Naive" = "seagreen", "CD8+/CD4+ Mem eff" = "yellow", "CD8+ CTL" = "dodgerblue4", "CD4+ CTL" = "orangered2")
# visualize B cell types separated by treatment group
DimPlot(object = Tcell_subcluster,  label.size = 4, raster = FALSE,  split.by = "orig.ident", cols = cols)

# find markers
CD45_Tcell_markers <- FindAllMarkers(Tcell_subcluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# arranges the data nicely
CD45_Tcell_markers %>% group_by(cluster) %>% top_n(n=25, wt = avg_log2FC)
# saves the markers in a comma-separated values table that can be accessed in excel
write.csv(CD45_Tcell_markers, file = "CD45_Tcell_markers.csv")
#save Tcell45_markers
saveRDS(CD45_Tcell_markers, file = "CD45_Tcell_markers.rds")
CD45_Tcell_markers <- readRDS("CD45_Tcell_markers.rds")

# identifying sub clusters with dotplot and key genes
plot<- DotPlot (object = Tcell_subcluster, dot.scale = 10, cols = c("lightgrey", "blue"), scale=T, 
         features=c("Cd4", "Cd8b1", "Cd8a", "Sell", "Dapl1", "S1pr1", "Ccr7", "Lef1","Icos", "Cxcr5", "Cd44", 
                    "Cxcr6", "Lgals3", "Itgb1", "Ccr2", "Itgae", "Foxp3", "Il2ra", "Ikzf2", 
                    "Ctla4", "Hif1a", "Nkg7", "Ccl5", "Rora",
                    "Gzmb","Tcf7", "Itga4", "Gzma", "Tnf", "Klrg1"))
plot + theme(axis.text.x = element_text(angle = 45, hjust=1))

# renaming T cell clusters 
new.cluster.ids.CD45Tid <- c("Helper T", "CD8+/CD4+ Naive", "Helper T", "CD8+/CD4+ Naive", 
                           "Helper T", "CD8+ Naive", "CD8+/CD4+ Mem eff", "CD8+ CTL", "CD4+ CTL")
names(new.cluster.ids.CD45Tid) <- levels(Tcell_subcluster)
Tcell_subcluster  <- RenameIdents(Tcell_subcluster, new.cluster.ids.CD45Tid)
# assign colors to each cell type
cols <- c("Helper T" = "orange", "CD8+/CD4+ Naive" = "skyblue", "CD8+ Naive" = "seagreen", "CD8+/CD4+ Mem eff" = "yellow", "CD8+ CTL" = "dodgerblue4", "CD4+ CTL" = "orangered2")
# visualize B cell types separated by treatment group
DimPlot(object = Tcell_subcluster,  label.size = 4, raster = FALSE,  split.by = "orig.ident", cols = cols)

# combine cell types into separate objects
HelperT <- subset(Tcell_subcluster, ident = "Helper T")
NaiveT <- subset(Tcell_subcluster, ident = c("CD8+/CD4+ Naive", "CD8+ Naive"))
Effector<- subset(Tcell_subcluster, ident = c("CD8+/CD4+ Mem eff", "CD8+ CTL", "CD4+ CTL"))
CTL <- subset (Tcell_subcluster, ident = c("CD8+ CTL", "CD4+ CTL"))
CD8CTL <- subset (Tcell_subcluster, ident = "CD8+ CTL")
CD4CTL <- subset (Tcell_subcluster, ident = "CD4+ CTL")


#
##
### heatmap for subtypes based off markers
##
#
#Naive markers
genes_of_interest <- c("Cd44", "Ccr7", "Cd4", "Lef1", "Sell", "Il7r", "S1pr1", "Dapl1")
#Helper markers
genes_of_interest <- c("Cd4", "Cxcr5", "Icos", "Bcl6", "Gata3")
#Mem eff markers
genes_of_interest<- c("Cd44","Id3", "Stat3", "Slamf6", "Cd69", "Cxcr3", "Lgals1", "Lgals3", "Cxcr6", "Ccr2", "Bhlhe40", "Tbx21","Bcl6", "Id2", "Rora", "Ifng", "Il21", "Il4")
#CTL markers
genes_of_interest<- c("Cd69", "Cd44", "Gzmb", "Gzma", "Prf1")

expression_data <- FetchData(Tcell_subcluster, vars = genes_of_interest)

#get metadata for treatment groups from object
treatment_groups <- unique(Tcell_subcluster@meta.data$orig.ident)
print(treatment_groups)
# Initialize a matrix to store the mean expression values
expression_by_group <- matrix(nrow = length(genes_of_interest), ncol = length(treatment_groups))
rownames(expression_by_group) <- genes_of_interest
colnames(expression_by_group) <- treatment_groups

# Calculate the mean expression for each gene in each treatment group
for (group in treatment_groups) {
  cells_in_group <- WhichCells(Tcell_subcluster, expression = orig.ident == group)
  expression_data_group <- GetAssayData(Tcell_subcluster, slot = "data")[genes_of_interest, cells_in_group, drop=FALSE]
  expression_by_group[, group] <- rowMeans(expression_data_group)
}
expression_by_group <- as.matrix(expression_by_group)

# Create the heatmap
pheatmap(expression_by_group, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         scale = "row",
         color = colorRampPalette(c("blue", "white", "red"))(50), name = "Expression")

### troubleshooting if index out of range 
treatment_groups <- unique(Tcell_subcluster@meta.data$orig.ident)
print(treatment_groups)

expression_data <- GetAssayData(Tcell_subcluster, slot = "data")
missing_genes <- setdiff(genes_of_interest, rownames(expression_data))
if (length(missing_genes) > 0) {
  stop("The following genes are missing from the expression data: ", paste(missing_genes, collapse = ", "))
}


#
##
### pathway analysis
##
#


# bargraph depicting the percent composition of each B cell type
dittoBarPlot(object = Tcell_subcluster, var = "RNA_snn_res.0.5", group.by = "orig.ident", 
             var.labels.rename = new.cluster.ids.CD45Tid, data.out = TRUE, x.reorder =c(4,1,2,5,3))

# Most gene signatures were obtained from the Gene Ontology (GO) database. Refer to my paper references for
# the past literature from which I obtained the gene signatures not found in the GO database. 

## helper T differentiation

# Create a named vector for significance labels
genes <- c("Anxa1", "Brd2", "Brd4", "Ccl19", "Ccr2", "Ccr7", "Ep300", "Hlx", "Il4ra", "Il6", "Il18", "Il23a", "Irf1", "Malt1", "Nfkbid", "Nfkbiz", "Nlrp3", "Prkcz", "Rara", "Ripk2", "Shb", "Socs5", "Tnfsf4")
HelperT <- AddModuleScore(HelperT, features=list(genes), name= "Helper_Differentiation")
my_colors <- c("Control" = "lightblue", "Adjuvant" = "darkorange", "CA170" = "forestgreen", "Kvax"= "purple", "Combo" = "gold")
VlnPlot(HelperT,features= "Helper_Differentiation1", pt.size = 0, group.by ="orig.ident", cols = my_colors, y.max = 0.75) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95)+
  geom_boxplot(width=0.1, fill="white")  +
  geom_signif(
    comparisons = list(c("Control", "Kvax")),
    annotations = "****",
    textsize = 5.5,
    map_signif_level = TRUE,
    y_position = 0.52,
    tip_length=0.03
  )+
  geom_signif(
    comparisons = list(c("Control", "CA170")),
    annotations = "****",
    textsize = 5.5,
    map_signif_level = TRUE,
    y_position = 0.42,
    tip_length=0.03
  )+
  geom_signif(
    comparisons = list(c("Control", "Combo")),
    annotations = "****",
    textsize = 5.5,
    map_signif_level = TRUE,
    y_position = 0.66,
    tip_length=0.03
  )+
  geom_signif(
    comparisons = list(c("CA170", "Combo")),
    annotations = "****",
    textsize = 5.5,
    map_signif_level = TRUE,
    y_position = 0.56,
    tip_length=0.03
  )+
  geom_signif(
    comparisons = list(c("Kvax", "Combo")),
    annotations = "****",
    textsize = 5.5,
    map_signif_level = TRUE,
    y_position = 0.45,
    tip_length=0.03
  )
# get median values
df <- FetchData(HelperT, vars = c("Helper_Differentiation1", "orig.ident"))
# Compute median for each group
medians <- aggregate(Helper_Differentiation1 ~ orig.ident, data = df, FUN = median)
# Print the median values
print(medians)

## helper T immune response 
genes <- c("Arid5a", "Ascl2", "Batf", "Brd2", "Brd4", "Card9", "Ccl20", "Ccr6", "Clec4n",  
           "Entpd7", "Ep300", "Ephb2", "Foxp3", "Il2", "Il4", "Il6", "Il6ra", "Il17ra", "Il23a", 
           "Il27ra", "Irf4", "Jak2", "Lgals9", "Loxl3", "Ly9", "Malt1",  
           "Nfkbid", "Nfkbiz", "Nlrp3", "Nlrp10", "Notch1", "Otud5", "Pf4", "Prkcq", "Rc3h1",
           "Rc3h2", "Rora", "Rorc", "Slamf6", "Smad7", "Stat3", "Tbx21", "Tgfb1", "Tnfsf18", "Traf3ip2", 
           "Tyk2", "Zbtb7b", "Zc3h12a", "Anxa1", "Arid5a", "Ascl2", "Bcl3", "Ccl19", "Ccr2", "Ccr7", 
           "Cracr2a", "Gadd45g", "Havcr2", "Hlx", "Hmgb1", "Hras", "Il1b", "Il1r1", "Il1rl1", "Il4", 
           "Il4ra", "Il12a", "Il12b", "Il12rb1", "Il18", "Il18bp", "Il18r1", "Il18rap", "Il23a", "Il23r", 
           "Il27", "Il27ra", "Il33", "Irf1", "Jak3", "Lef1", "Mtor", "Nfkbiz", "Nlrp10", "Pla2g4a", "Relb",
           "Ripk2", "Sema4a", "Slamf1", "Slc11a1", "Socs5", "Spn", "Stat4", "Stat6", "Tbx21", "Tmem98", 
           "Tnfsf4", "Traf6", "Vegfa", "Xcl1","Cd74", "Cd81", "Clcf1", "Dennd1b", "Gata3", "Ido1", "Il4", 
           "Il4ra", "Il6", "Il18", "Il33", "Nlrp3", "Nod2", "Prkcz", "Rara", "Rsad2", "Tnfsf4", "Xcl1")
HelperT <-AddModuleScore(HelperT, features=list(genes), name= "HelperT_immuneresponse")
cols <- c("Control" = "lightblue", "Adjuvant" = "darkorange", "CA170" = "forestgreen", "Kvax"= "purple", "Combo" = "gold")
VlnPlot(HelperT,features= "HelperT_immuneresponse1", pt.size = 0, group.by ="orig.ident", col=cols, y.max = 0.45)+
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95)+
  geom_boxplot(width=0.1, fill="white")+
  geom_signif(
    comparisons = list(c("Control", "Kvax")),
    annotations = "****",
    textsize = 5.5,
    map_signif_level = TRUE,
    y_position = 0.31,
    tip_length=0.03
  )+
  geom_signif(
    comparisons = list(c("Control", "CA170")),
    annotations = "****",
    textsize = 5.5,
    map_signif_level = TRUE,
    y_position = 0.23,
    tip_length=0.03
  )+
  geom_signif(
    comparisons = list(c("Control", "Combo")),
    annotations = "****",
    textsize = 5.5,
    map_signif_level = TRUE,
    y_position = 0.39,
    tip_length=0.03
  )+
  geom_signif(
    comparisons = list(c("CA170", "Combo")),
    annotations = "****",
    textsize = 5.5,
    map_signif_level = TRUE,
    y_position = 0.25,
    tip_length=0.03
  )+
  geom_signif(
    comparisons = list(c("Kvax", "Combo")),
    annotations = "****",
    textsize = 5.5,
    map_signif_level = TRUE,
    y_position = 0.33,
    tip_length=0.03
  )
df <- FetchData(HelperT, vars = c("HelperT_immuneresponse1", "orig.ident"))
medians <- aggregate(HelperT_immuneresponse1 ~ orig.ident, data = df, FUN = median)
print(medians)

## Helper T cytokine production
genes <- c("Cd81", "Dennd1b", "Gata3", "Il4", "Il6", "Nlrp3", "Prkcz", "Rsad2", "Arid5a", 
           "Il1b", "Il1r1", "Il12a", "Il12b", "Il18", "Il18r1", "Il18rap", "Slamf1", "Tbx21", "Xcl1")
HelperT <-AddModuleScore(HelperT, features=list(genes), name= "HelperT_cytokineproduction")
cols <- c("Control" = "lightblue", "Adjuvant" = "darkorange", "CA170" = "forestgreen", "Kvax"= "purple", "Combo" = "gold")
VlnPlot(HelperT,features= "HelperT_cytokineproduction1", pt.size = 0, group.by ="orig.ident", col=cols, y.max = 1)+
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95)+
  geom_boxplot(width=0.1, fill="white")+
  geom_signif(
    comparisons = list(c("Control", "CA170")),
    annotations = "****",
    textsize = 5.5,
    map_signif_level = TRUE,
    y_position = 0.50,
    tip_length=0.03
  )+
  geom_signif(
    comparisons = list(c("Control", "Kvax")),
    annotations = "****",
    textsize = 5.5,
    map_signif_level = TRUE,
    y_position = 0.61,
    tip_length=0.03
  )+
  geom_signif(
    comparisons = list(c("Control", "Combo")),
    annotations = "****",
    textsize = 5.5,
    map_signif_level = TRUE,
    y_position = 0.76,
    tip_length=0.03
  )+
  geom_signif(
    comparisons = list(c("CA170", "Combo")),
    annotations = "****",
    textsize = 5.5,
    map_signif_level = TRUE,
    y_position = 0.9,
    tip_length=0.03
  )+
  geom_signif(
    comparisons = list(c("Kvax", "Combo")),
    annotations = "****",
    textsize = 5.5,
    map_signif_level = TRUE,
    y_position = 0.65,
    tip_length=0.03
  )
df <- FetchData(HelperT, vars = c("HelperT_cytokineproduction1", "orig.ident"))
medians <- aggregate(HelperT_cytokineproduction1 ~ orig.ident, data = df, FUN = median)
print(medians)

## T cell extravasation 
genes <-c("Ccl2", "Ccl5", "Ccr2", "Cd99l2", "Crk", "Crkl", "F11r", "Fadd", "Icam1", "Il27ra", "Itgal", "Med23", "Ripk3")
Tcell_subcluster <- AddModuleScore(Tcell_subcluster, features=list(genes), name= "T_cell_extravasation")
cols <- c("Control" = "lightblue", "Adjuvant" = "darkorange", "CA170" = "forestgreen", "Kvax"= "purple", "Combo" = "gold")
VlnPlot(Tcell_subcluster,features= "T_cell_extravasation1", pt.size = 0, group.by ="orig.ident", col = cols, y.max = 1.3)+
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95)+
  geom_boxplot(width=0.1, fill="white")+
  geom_signif(
    comparisons = list(c("Control", "Kvax")),
    annotations = "****",
    textsize = 5.5,
    map_signif_level = TRUE,
    y_position = 0.8,
    tip_length=0.03
  )+
  geom_signif(
    comparisons = list(c("Control", "Combo")),
    annotations = "****",
    textsize = 5.5,
    map_signif_level = TRUE,
    y_position = 0.95,
    tip_length=0.03
  )+
  geom_signif(
    comparisons = list(c("CA170", "Combo")),
    annotations = "****",
    textsize = 5.5,
    map_signif_level = TRUE,
    y_position = 1.15,
    tip_length=0.03
  )+
  geom_signif(
    comparisons = list(c("Kvax", "Combo")),
    annotations = "****",
    textsize = 5.5,
    map_signif_level = TRUE,
    y_position = 1.01,
    tip_length=0.03
  )
df <- FetchData(Tcell_subcluster, vars = c("T_cell_extravasation1", "orig.ident"))
medians <- aggregate(T_cell_extravasation1 ~ orig.ident, data = df, FUN = mean)
print(medians)

## CD8 alpha beta extravasation
genes <- c("Ccl5", "Ccr2", "Ripk3", "Ccl3", "Icam1", "Vcam1")
CTL <- AddModuleScore(CTL, features=list(genes), name= "CD8_alphabeta_T_extravasation")
cols <- c("Control" = "lightblue", "Adjuvant" = "darkorange", "CA170" = "forestgreen", "Kvax"= "purple", "Combo" = "gold")
VlnPlot(CTL,features= "CD8_alphabeta_T_extravasation1", pt.size = 0, group.by ="orig.ident", col = cols, y.max = 1.7)+
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95)+
  geom_boxplot(width=0.1, fill="white") +
  geom_signif(
    comparisons = list(c("Control", "Kvax")),
    annotations = "****",
    textsize = 5.5,
    map_signif_level = TRUE,
    y_position = 1.25,
    tip_length=0.02
  ) +
  geom_signif(
    comparisons = list(c("Control", "CA170")),
    annotations = "*",
    textsize = 5.5,
    map_signif_level = TRUE,
    y_position = 1.05,
    tip_length=0.02
  ) +
  geom_signif(
    comparisons = list(c("Control", "Combo")),
    annotations = "****",
    textsize = 5.5,
    map_signif_level = TRUE,
    y_position = 1.45,
    tip_length=0.02
  ) 
df <- FetchData(CTL, vars = c("CD8_alphabeta_T_extravasation1", "orig.ident"))
medians <- aggregate(CD8_alphabeta_T_extravasation1 ~ orig.ident, data = df, FUN = mean)
print(medians)

## T cell chemotaxis
genes <- c("Adam10", "Adam17", "Ccl3", "Ccl5", "Ccl21a", "Ccl26", "Ccl27a", "Ccr2", "Ccr7", "Cxcl10", "Cxcl11", "Cxcl12", "Cxcl13", "Cxcl16", "Cxcr3", "Gpr183", "Lgals9", "Oxsr1", "Plec", "Slc12a2", "Stk39", "Tmem102", "Tnfsf14", "Wnk1", "Wnt5a", "Xcl1")
Tcell_subcluster <- AddModuleScore(Tcell_subcluster, features=list(genes), name= "Chemotaxis")
cols <- c("Control" = "lightblue", "Adjuvant" = "darkorange", "CA170" = "forestgreen", "Kvax"= "purple", "Combo" = "gold")
VlnPlot(Tcell_subcluster,features= "Chemotaxis1", pt.size = 0, group.by ="orig.ident", col = cols, y.max= 1)+
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95)+
  geom_boxplot(width=0.1, fill="white")  +
  geom_signif(
    comparisons = list(c("Control", "CA170")),
    annotations = "****",
    textsize = 5.5,
    map_signif_level = TRUE,
    y_position = 0.65,
    tip_length=0.02
  ) +
  geom_signif(
    comparisons = list(c("Control", "Kvax")),
    annotations = "****",
    textsize = 5.5,
    map_signif_level = TRUE,
    y_position = 0.78,
    tip_length=0.02
  ) +
  geom_signif(
    comparisons = list(c("Control", "Combo")),
    annotations = "****",
    textsize = 5.5,
    map_signif_level = TRUE,
    y_position = 0.9,
    tip_length=0.02
  ) +
  geom_signif(
    comparisons = list(c("CA170", "Combo")),
    annotations = "****",
    textsize = 5.5,
    map_signif_level = TRUE,
    y_position = 0.6,
    tip_length=0.02
  ) +
  geom_signif(
    comparisons = list(c("Kvax", "Combo")),
    annotations = "****",
    textsize = 5.5,
    map_signif_level = TRUE,
    y_position =0.74,
    tip_length=0.02
  )
df <- FetchData(Tcell_subcluster, vars = c("Chemotaxis1", "orig.ident"))
means <- aggregate(Chemotaxis1 ~ orig.ident, data = df, FUN = mean)
print(means)

## CTL chemotaxis
genes <-c("Ccr5", "Cxcr3", "Cxcr6", "Ccl5", "Cxcl10", "Cxcl9", "Ccl3", "Ccl4", "Itgal", "Icam1", "Icam2", "Sele", "Sell", "S1pr1", "Cd44")
CTL <- AddModuleScore(CTL, features=list(genes), name= "CTL_Chemotaxis")
cols <- c("Control" = "lightblue", "Adjuvant" = "darkorange", "CA170" = "forestgreen", "Kvax"= "purple", "Combo" = "gold")
VlnPlot(CTL,features= "CTL_Chemotaxis1", pt.size = 0, group.by ="orig.ident", col =cols, y.max = 1.2)+
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95)+
  geom_boxplot(width=0.1, fill="white") + 
  geom_signif(
    comparisons = list(c("Control", "Kvax")),
    annotations = "****",
    textsize = 5.5,
    map_signif_level = TRUE,
    y_position = 1,
    tip_length=0.02)  +
  geom_signif(
    comparisons = list(c("Control", "Combo")),
    annotations = "****",
    textsize = 5.5,
    map_signif_level = TRUE,
    y_position = 0.82,
    tip_length=0.02
  ) 
df <- FetchData(CTL, vars = c("CTL_Chemotaxis1", "orig.ident"))
medians <- aggregate(CTL_Chemotaxis1 ~ orig.ident, data = df, FUN = median)
print(medians)

## CTL activity
genes <- c("Ceacam1", "Nckap1l", "Rab27a", "Stx11","Il2", "Il4", "Il7", "Il10", "Il15", "Ifng", "Cd28",
           "Mr1", "Muc4", "Slc22a13","Prf1", "Gzmb", "Gzma", "Ifng", "Fasl", "Tnf", "Nkg7", "Cd8a", 
           "Tbx21", "Eomes", "Stat1", "Ccl3", "Ccl4")
CTL <-AddModuleScore(CTL, features=list(genes), name= "CTL_degranulation_activation_tumorcytotoxicity")
cols <- c("Control" = "lightblue", "Adjuvant" = "darkorange", "CA170" = "forestgreen", "Kvax"= "purple", "Combo" = "gold")
VlnPlot(CTL,features= "CTL_degranulation_activation_tumorcytotoxicity1", pt.size = 0, group.by ="orig.ident", col=cols, y.max = 0.7)+
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95)+
  geom_boxplot(width=0.1, fill="white")+
  geom_signif(
    comparisons = list(c("Control", "Combo")),
    annotations = "**",
    textsize = 5.5,
    map_signif_level = TRUE,
    y_position = 0.62,
    tip_length=0.02
  ) 
df <- FetchData(CTL, vars = c("CTL_degranulation_activation_tumorcytotoxicity1", "orig.ident"))
medians <- aggregate(CTL_degranulation_activation_tumorcytotoxicity1 ~ orig.ident, data = df, FUN = mean)
print(medians)

#
##
### statistical testing 
##
#


# helper T immune response
data <- data.table(x1 = HelperT$orig.ident,
                   x2=HelperT$HelperT_immuneresponse1)
write.csv(data, "HelperT_immuneresponse1_for_stats.csv")

# helper T cytokine production
data <- data.table(x1 = HelperT$orig.ident,
                   x2=HelperT$HelperT_cytokineproduction1)
write.csv(data, "HelperT_cytokineproduction1_for_stats.csv")

# helper T activation
data <- data.table(x1 = HelperT$orig.ident,
                   x2=HelperT$Helper_Activation1)
write.csv(data, "Helper_Helper_Activation1_for_stats.csv")

# T cell extravastion
data <- data.table(x1 = Tcell_subcluster$orig.ident,
                   x2=Tcell_subcluster$T_cell_extravasation1)
write.csv(data, "Tcellsubcluster_Tcellextravasation1_for_stats.csv")

# CD8 T cell extravasation
data <- data.table(x1 = CTL$orig.ident,
                   x2 = CTL$CD8_alphabeta_T_extravasation1)
write.csv(data, "CTL_CD8_alphabeta_T_extravasation1_for_stats.csv")

# T cell chemotaxis
data <- data.table(x1 = Tcell_subcluster$orig.ident,
                   x2=Tcell_subcluster$Chemotaxis1)
write.csv(data, "Tcellsubcluster_Chemotaxis1_for_stats.csv")

# CTL chemotaxis
data <- data.table(x1 = CTL$orig.ident,
                   x2=CTL$CTL_Chemotaxis1)
write.csv(data, "CTL_CTL_Chemotaxis1_for_stats.csv")

# CTL activity
data <- data.table(x1 = CTL$orig.ident,
                   x2=CTL$CTL_degranulation_activation_tumorcytotoxicity_cytokines1)
write.csv(data, "CTL_degranulation_activation_tumorcytotoxicity_cytokines_stats.csv")


#
##
### heatmap
##
#


## Helper T cell activation + cytokine production
genes_of_interest <- c("Cd81", "Dennd1b", "Gata3", "Il4", "Il6", "Nlrp3", "Prkcz", "Rsad2", "Arid5a", 
              "Il1b", "Il1r1", "Il12a", "Il18", "Il18r1", "Il18rap", "Slamf1", "Tbx21",
              "Cebpb", "Tcirg1", "Xcl1", "Prkcq", "Tnfsf4", "Ifng", "Il12rb1", "Il17a", "Il22", 
              "Stat3")
## immune cell adhesion and migration 
genes_of_interest <- c ("Cd44", "Cxcr4", "Ccr1", "Ccr5", "Ccr7")

expression_data <- FetchData(HelperT, vars = genes_of_interest)

#get metadata for treatment groups from object
treatment_groups <- unique(HelperT@meta.data$orig.ident)
print(treatment_groups)
# Initialize a matrix to store the mean expression values
expression_by_group <- matrix(nrow = length(genes_of_interest), ncol = length(treatment_groups))
rownames(expression_by_group) <- genes_of_interest
colnames(expression_by_group) <- treatment_groups

# Calculate the mean expression for each gene in each treatment group
for (group in treatment_groups) {
  cells_in_group <- WhichCells(HelperT, expression = orig.ident == group)
  expression_data_group <- GetAssayData(HelperT, slot = "data")[genes_of_interest, cells_in_group, drop=FALSE]
  expression_by_group[, group] <- rowMeans(expression_data_group)
}

expression_by_group <- as.matrix(expression_by_group)

# Create the heatmap
pheatmap(expression_by_group, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         scale = "row",
         color = colorRampPalette(c("blue", "white", "red"))(50), name = "Expression")

### troubleshooting if index out of range 
treatment_groups <- unique(HelperT@meta.data$orig.ident)
print(treatment_groups)

expression_data <- GetAssayData(HelperT, slot = "data")
missing_genes <- setdiff(genes_of_interest, rownames(expression_data))
if (length(missing_genes) > 0) {
  stop("The following genes are missing from the expression data: ", paste(missing_genes, collapse = ", "))
}



#
##
####
#######
########### saving R file as PDF
#######
####
##
#

knitr::stitch('scrna_seq_LC - Copy.r')


