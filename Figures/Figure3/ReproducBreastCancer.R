### 
library(ggplot2); library(gridExtra)
library(tidyverse)
library(reshape2)


touch <- readRDS('touchstone_metrics.rds')
public <- readRDS('public_metrics.rds')
## updating manually suzuki sample
public[public$dataset %in% 'CosMx_BreastCancer_Suzuki', ]$TxPerCell <- 596.68
public[public$dataset %in% 'CosMx_BreastCancer_Suzuki', ]$TxPerNuc <- 360.692
  
  
  
## BR reproducibility - Touchstone vs Public dataset

breast.public <-  public[public$TissueType %in% 'BreastCancer',]
colnames(breast.public)[1] <- 'sample_id'
breast.touch <- touch[touch$TissueType %in% 'BR', ]


## PCA using all metrics from BR public and touch
which_metrics <- colnames(touch)[c(6,7,11:15,17,21,24:25)]

breast.touch$Origin <- 'Touchstone'
breast.public$Origin <- 'Public'
both <- rbind(breast.touch[, c('Origin','platform', 'TissueType', which_metrics)], breast.public[, c('Origin','platform', 'TissueType', which_metrics)])
both$Origin <- factor(both$Origin, levels = c("Touchstone", "Public"))


pca1 <- prcomp(scale(both[, which_metrics]))
pca <- as.data.frame(pca1$x)
#rownames(pca) <- make.names(public$dataset, unique = T)
pca$Platform <- both$platform
pca$TissueType <- both$TissueType
pca$Origin <- both$Origin

#pca$Institution <- touch$Institution
p <- ggplot(pca, aes(PC1, PC2, label = rownames(pca), col = Origin, shape = Platform)) + geom_point(size =8 ) +
   xlab(paste0('PC1 ', summary(pca3)$importance[2, 1] * 100, '%') ) +
  ylab(paste0('PC2 ', summary(pca3)$importance[2, 2] * 100, '%') )+ ggtitle('Touchstone and Public - Breast Cancer datasets') +
  scale_color_manual(values = c('#87ceeb','#f08080' )) + theme_bw() +  theme(text=element_text(size=20))

#pdf('plots/BreastCancerReproduc_AllMetrics_PCA_PublicTouchstone.pdf', 10,5)
png('plots/BreastCancerReproduc_AllMetrics_PCA_PublicTouchstone.png', width=11,height=8, units = 'in', res= 300)
p
dev.off()


## Plot metrics both datasets by platform
metrics <-colnames(both)[-c(1:3)]

plots <- list()
for(i in metrics) {
  p <- ggplot(both, aes_string(y = i, x = "Origin", fill = "Origin"))  + geom_boxplot()  + theme_bw() +
    scale_fill_manual(values = c('#87ceeb','#f08080' )) + theme(text=element_text(size=20))
  plots[[i]] <- p
}
legend <- get_legend(plots[[1]])
plots <- lapply(plots, function(p) p + theme(legend.position = "none"))
library(cowplot)
# Combine all plots into one figure with one legend
combined_plot <- plot_grid(plotlist = plots, ncol = 4)
p2 <- plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, 0.1))
#pdf("plots/BreastCancerReproduc_AllMetrics_PublicTouchstone.pdf", 18, 10)
png("plots/BreastCancerReproduc_AllMetrics_PublicTouchstone.png", width=18,height = 12, res = 300, units ='in')
p2
dev.off()

##### Biology related reproducibility

br.objs <- list()
for(i in breast.touch$sample_id) {
  br.objs[[i]] <- readSpatial(sample_id = i, path = paste0('/mnt/scratch1/Touchstone_data/new_data/', i), platform = breast.touch[breast.touch$sample_id %in% i, ]$platform,seurat = T)
}

br.objs[['Xenium_BreastCancer_Suzuki']] <- readRDS('/mnt/scratch1/Touchstone_data/public_data/Xenium_BreastCancer_Suzuki/seurat_obj.rds')
br.objs[['Xenium_BreastCancer_Suzuki']][['RNA']] <- br.objs[['Xenium_BreastCancer_Suzuki']][['Nanostring']]
br.objs[['Xenium_BreastCancer_Lobular']] <- readSpatial(sample_id = 'Xenium_BreastCancer_Lobular', path = '/mnt/scratch1/Touchstone_data/public_data/Xenium_BreastCancer_Lobular/', 
                                                        platform = 'Xenium', seurat = T)
br.objs[['Xenium_BreastCancer_DuctalCarcinoma']] <- readSpatial(sample_id = 'Xenium_BreastCancer_DuctalCarcinoma', path = '/mnt/scratch1/Touchstone_data/public_data/Xenium_BreastCancer_DuctalCarcinoma/', 
                                                        platform = 'Xenium', seurat = T)

br.objs[['CosMx_BreastCancer_Suzuki']] <- readRDS('/mnt/scratch1/Touchstone_data/public_data/CosMx_BreastCancer_Suzuki/seurat_obj.rds')
br.objs[['CosMx_BreastCancer_Suzuki']][['RNA']] <- br.objs[['CosMx_BreastCancer_Suzuki']][['Nanostring']]


#sys_probes <- grep("SystemControl", rownames(br.objs[['CosMx_BreastCancer_Suzuki']]), value=T)
neg_probes <- grep("NegP", rownames(br.objs[['CosMx_BreastCancer_Suzuki']]), value=T)
br.objs[['CosMx_BreastCancer_Suzuki']][["ControlProbe"]] <- CreateAssayObject(
  counts = br.objs[['CosMx_BreastCancer_Suzuki']][["Nanostring"]]$counts[neg_probes,]
)
## Make "Nanostring" assay
tx_probes <- rownames(br.objs[['CosMx_BreastCancer_Suzuki']])[!rownames(br.objs[['CosMx_BreastCancer_Suzuki']]) %in%neg_probes]
br.objs[['CosMx_BreastCancer_Suzuki']][["RNA"]] <- CreateAssayObject(
  counts = br.objs[['CosMx_BreastCancer_Suzuki']][["Nanostring"]]$counts[tx_probes,]
)
DefaultAssay(br.objs[['CosMx_BreastCancer_Suzuki']]) <- "RNA"
br.objs[['CosMx_BreastCancer_Suzuki']][["Nanostring"]] <- NULL

##### Cell Metadata
print("Getting additional cell metadata")
cell_meta <- data.table::fread(
  file.path(path, list.files(path, pattern="*metadata_file*"))
)
cell_meta$cell_ID <- paste0(cell_meta$cell_ID, '_', cell_meta$fov)
cell_meta <- cell_meta[cell_meta$cell_ID %in% colnames(br.objs[['CosMx_BreastCancer_Suzuki']]), ]
#It's excessive, but we'll add all metadata into the object
br.objs[['CosMx_BreastCancer_Suzuki']]@meta.data <- cbind(br.objs[['CosMx_BreastCancer_Suzuki']]@meta.data, cell_meta)
br.objs[['CosMx_BreastCancer_Suzuki']]@meta.data$fov <- factor(paste0("FOV", br.objs[['CosMx_BreastCancer_Suzuki']]@meta.data$fov))
#seu_obj@meta.data$cell_area <- seu_obj$Area.um2
br.objs[['CosMx_BreastCancer_Suzuki']]@meta.data$transcript_counts <- br.objs[['CosMx_BreastCancer_Suzuki']]$nCount_RNA
br.objs[['CosMx_BreastCancer_Suzuki']]@meta.data$negprobe_counts <- br.objs[['CosMx_BreastCancer_Suzuki']]$nCount_ControlProbe


## Get autocorrelation


## Get autocorrelation
library(SingleCellExperiment); library(SpatialExperiment); library(SpatialFeatureExperiment)
library(BiocParallel)
morans <- lapply(br.objs, getMorans)

## necessary to work getMorans
coords <- GetTissueCoordinates(br.objs[['CosMx_BreastCancer_Suzuki']])
coords <- as.matrix(coords[, 1:2])
colnames(coords) <- c("Tissue_1", "Tissue_2")
rownames(coords) <- colnames(br.objs[['CosMx_BreastCancer_Suzuki']])
br.objs[['CosMx_BreastCancer_Suzuki']][["tissue"]] <- CreateDimReducObject(coords, key="Tissue_", assay="RNA")
br.objs[['CosMx_BreastCancer_Suzuki']]$sample_id <- 'CosMx_BreastCancer_Suzuki'
br.objs[['CosMx_BreastCancer_Suzuki']]$platform <- 'CosMx'

tmp <- getMorans(br.objs[['CosMx_BreastCancer_Suzuki']])
morans[['CosMx_BreastCancer_Suzuki']] <- tmp

# Step 3: Use Reduce and intersect to find common values in the specified column
common_values <- Reduce(intersect, lapply(morans, function(df) df[['gene']]))

# Step 4: Filter data frames based on common values
common_morans <- lapply(morans, function(df) df[df[['gene']] %in% common_values, ])
common_morans <- do.call(rbind, common_morans)
common_morans$Origin <- as.factor(common_morans$sample_id)
levels(common_morans$Origin) <- c('Public',rep('Touchstone', 8), rep('Public', 3))
common_morans$Origin <- as.character(common_morans$Origin)
common_morans$Origin <- factor(common_morans$Origin, levels = c("Public", "Touchstone"))


e <- ggplot(common_morans, aes(x = Origin, y = value, fill = Origin))


# Combine with box plot to add median and quartiles
# Change fill color by groups, remove legend
p3 <- e + geom_violin() + 
  geom_boxplot(width = 0.1) +
  theme(legend.position = "none") + 
  xlab("Moran's I") + ylab("Sample origin") +
  theme_bw() + 
  theme(text = element_text(size = 20)) + 
  scale_fill_manual(values = c('#f08080', '#87ceeb'))  + coord_flip()

#pdf('plots/BreastCancerReproduc_MoransI_PublicTouchstone.pdf', 10, 8)
png('plots/BreastCancerReproduc_MoransI_PublicTouchstone.png', width = 11,height =  8, units ='in', res =300)
p3
dev.off()


#### annotation using agnostic reference 

brca_wu <- readRDS("/media/ResearchHome/plummgrp/home/common/Touchstone_data/brca_wu.rds")
brca_wu$celltype_major <- as.character(brca_wu$celltype_major)

#br.objs[1:11] <- lapply(br.objs[1:11], function(x) { x <- annotateData(x, ref = brca_wu, celltype_meta = 'celltype_major')})
br.objs <- lapply(br.objs, function(x) { x <- annotateData(x, ref = brca_wu, celltype_meta = 'celltype_major')})



prop.cellt <- lapply(br.objs , function(x) {
  x <- as.data.frame( prop.table(table(x$celltype_pred)))
})

# Use lapply to add a new column with the name of the data frame
prop.cellt <- lapply(names(prop.cellt), function(name) {
  df <- prop.cellt[[name]]
  df$sample_id <- name
  return(df)
})

prop.cellt <- Reduce(rbind, prop.cellt)
prop.cellt$Origin <- as.factor(prop.cellt$sample_id)
levels(prop.cellt$Origin) <- c('Public',rep('Touchstone', 8), rep('Public', 3))

## sd per celltype
sd_prop <- prop.cellt%>%
  group_by(Origin, Var1) %>%
  summarise(sd_Freq = sd(Freq, na.rm = TRUE))

sd_public<- sd_prop[sd_prop$Origin == 'Public',]
sd_public[order(sd_public$sd_Freq, decreasing = T), ]

sd_touch<- sd_prop[sd_prop$Origin == 'Touchstone',]
sd_touch[order(sd_touch$sd_Freq, decreasing = T), ]


#### BAR PLOT CELL PROPORTIONS FOR ALL SAMPLES ###

color_map <- brewer.pal(9, 'Set1')
names(color_map) <- levels(as.factor(prop.cellt$Var1))


p5 <- ggplot(prop.cellt, aes(Freq, sample_id, fill = Var1)) + geom_bar(stat='identity') + #facet_grid(rows=  prop.cellt$Origin,scales="free") + 
  scale_fill_manual(values = brewer.pal(9, 'Set1')) + theme_bw() + theme(text = element_text(size = 20)) +
  xlab('Cell type proportion') + ylab('') 
pdf('ReproducBreastCancer_CellTProp_BarPlots.pdf', 15, 6 )
p5
dev.off()

## Representative ImageDimPlot - Public and Touchstone ##
png('plots/PublicBreastCancer_CellType_ImagePlot.png', width = 10, height = 5, units = 'in', res = 150)
p4<- ImageDimPlot(br.objs[['Xenium_BreastCancer_DuctalCarcinoma']], group.by = 'celltype_pred', cols = color_map, size = 0.4, dark.background = F) & NoLegend()
p4
dev.off()

png('plots/TouchstoneBreastCancer_CellType_ImagePlot.png', width = 10, height = 5, units = 'in', res = 150)
p4.1 <- ImageDimPlot(br.objs[['TOUCHSTONE_UOA_XR_FFPE_BR_2_SYD_C_R1']], group.by = 'celltype_pred', cols = color_map, size = 0.4, dark.background = F) & NoLegend()
p4.1
dev.off()

### 


## All plots in one
p1 + p2 + p3 + p4 + p4.1 + p5

# Plot all ImageDimPlot
plots <- list()
for( i in names(br.objs)) {
  plots[[i]] <- ImageDimPlot(br.objs[[i]], group.by = 'celltype_pred', cols = color_map, size = 0.4, dark.background = F) & NoLegend()
}
do.call(grid.arrange, c(plots, ncol = 5)) 





