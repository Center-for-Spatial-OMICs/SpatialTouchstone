## Plotting Touchstone metrics and public datasets metrics
library(ggplot2); library(RColorBrewer); library(gridExtra)

public <- readRDS('public_metrics.rds')
public[122:127,]$expMat <- gsub('cell_feature_matrix.h5', '', public[122:127,]$expMat)

touch <- readRDS('touchstone_metrics.rds')

# updated metrics (changed the way capturing noise from neg probes) - May 2024
touch$specificityFDR <- mapply(getGlobalFDR, tx_file = touch$tx_file, 
                               platform = touch$platform, cellSegMeta =  touch$cell_meta, 
                               MoreArgs = list(features = NULL))

touch$SigNoiseRatio <- mapply(getMeanSignalRatio, expMat = touch$expMat, 
                               platform = touch$platform, 
                               MoreArgs = list(features = NULL))

touch$DynamicRange <- mapply(getMaxRatio, expMat = touch$expMat, 
                              platform = touch$platform, 
                              MoreArgs = list(features = NULL))

## adding PA XR samples
basePath = '/mnt/scratch1/Touchstone_data/new_data/'
df_samples <- data.frame(sample_id = dir(basePath, pattern = 'TOUCHSTONE_WCM_XR_FFPE_PA'),
                         platform = "Xenium", expMat = NA, tx_file = NA, cell_meta = NA)

df_samples$expMat <- paste0(basePath,df_samples$sample_id, '/cell_feature_matrix/')
df_samples$tx_file <- paste0(basePath,df_samples$sample_id, '/transcripts.csv.gz')
df_samples$cell_meta <- paste0(basePath,df_samples$sample_id, '/cells.csv.gz')

results <- getAllMetrics(df_samples)
results$TxPerCellNorm <- results$TxPerCell / results$PanelSize
results$TissueType <- 'PA'
results$Institution <- 'WCM'
results$ComplexityNorm <- results$Complexity / results$PanelSize
results$Complexity_CommonGenes <- NA
results$Complexity_CommonGenesNorm <- NA
results$TxPerNucNorm <- results$TxPerNuc / results$PanelSize
results$Disease <- 'N'
results$Origin <- 'Touchstone'

#merge touchstone with new samples -- 
results <- results[, colnames(results) %in% colnames(touch)]
touch <- rbind(touch, results)
#saveRDS(touch, 'touchstone_metrics.rds')

###




# Add disease type
library(stringr)
touch$Disease <- str_extract(touch$sample_id, pattern = '_N_|_C_') %>% gsub("\\_", "", .)

public$specificityFDR <- mapply(getGlobalFDR, tx_file = public$t , platform = public$platform)
public$SigNoiseRatio <- mapply(getMeanSignalRatio, expMat = public$expMat, 
                              platform = public$platform, 
                              MoreArgs = list(features = NULL))

public$DynamicRange <- mapply(getMaxRatio, expMat = public$expMat, 
                             platform = public$platform, 
                             MoreArgs = list(features = NULL))

public$Complexity <- mapply(getComplexity, expMat = public$expMat , platform = public$platform)
public$DynamicRange <- mapply(getMaxRatio, expMat = public$expMat , platform = public$platform)
public$PanelSize <- mapply(getPanelSize, expMat = public$expMat , platform = public$platform)
public$ComplexityNorm <- public$Complexity / public$PanelSize
public$TxPerCellNorm <- public$TxPerCell / public$PanelSize
public$TxPerNucNorm<- public$TxPerNuc / public$PanelSize


# Add tissueType from sample name
public$dataset
public$TissueType <- gsub('CosMx_|Xenium_','',public$dataset)
public$TissueType <- gsub('\\_.*','',public$TissueType)
public$TissueType <- as.factor(public$TissueType)
levels(public$TissueType)[c(2:3,8,11)] <- c('Brain', 'Brain', 'Lung', 'Pancreas')
public$TissueType <- as.character(public$TissueType)
saveRDS(public, 'public_metrics.rds')
#write.csv(public[, -c(3:5)], file = 'public_datasets_metrics.csv')

## Plot metrics public datasets
metrics <- colnames(public)[c(7:17, 19:20, 22)]
plots <- list()
for(i in metrics) {
  p <- ggplot(public, aes_string(y = i, x = "platform", fill = "TissueType"))  + geom_boxplot() + scale_fill_manual(values = brewer.pal(8, 'Paired'))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  plots[[i]] <- p
}

do.call(grid.arrange, c(plots, ncol = 5)) 


## Plot metrics touchstone datasets
plots <- list()
for(i in metrics) {
  p <- ggplot(touch, aes_string(y = i, x = "Disease", fill = "TissueType"))  + geom_boxplot() + scale_fill_manual(values = brewer.pal(6, 'Paired'))
  plots[[i]] <- p
}

do.call(grid.arrange, c(plots, ncol = 5)) 



# PCA of metrics #
library(plotly); library(RColorBrewer); library(htmlwidgets); library(tidyr)


#which_metrics <- colnames(touch)[c(7,9,11:15,17,21,24:25)]
which_metrics <- colnames(touch)[c(6,7,11:15,17,21,24:25)]


pca2 <- prcomp(scale(touch[, which_metrics]))
pca <- as.data.frame(pca2$x)
rownames(pca) <- touch$sample_id
pca$Platform <- touch$platform
pca$TissueType <- touch$TissueType
pca$Institution <- touch$Institution
p <- ggplot(pca, aes(PC1, PC2, label = rownames(pca), colour = TissueType, shape = Platform)) + geom_point(size =5) +
  scale_colour_manual(values = brewer.pal(6, 'Paired')) + xlab(paste0('PC1 ', summary(pca2)$importance[2, 1] * 100, '%') ) +
  ylab(paste0('PC2 ', summary(pca2)$importance[2, 2] * 100, '%') )+ ggtitle('PCA of all metrics - color by Tissue type')
p2 <- ggplot(pca, aes(PC1, PC2, label = rownames(pca), colour = Institution, shape = Platform)) + geom_point(size =5) +
  scale_colour_manual(values = brewer.pal(4, 'Set2'))+ xlab(paste0('PC1 ', summary(pca2)$importance[2, 1] * 100, '%') ) +
  ylab(paste0('PC2 ', summary(pca2)$importance[2, 2] * 100, '%') ) + ggtitle('PCA of all metrics - color by Institution')

pp <- subplot(ggplotly(p) , ggplotly(p2) , shareX = T, shareY = T, titleX = T,
        titleY = T)
saveWidget(pp, "PCA_TouchstoneMetricsScaled.html")
#--------#


pca3 <- prcomp(scale(public[, which_metrics]))
pca <- as.data.frame(pca3$x)
rownames(pca) <- make.names(public$dataset, unique = T)
pca$Platform <- public$platform
pca$TissueType <- public$TissueType
#pca$Institution <- touch$Institution
p <- ggplot(pca, aes(PC1, PC2, label = rownames(pca), colour = TissueType, shape = Platform)) + geom_point(size =5) +
  scale_colour_manual(values = brewer.pal(8, 'Dark2')) + xlab(paste0('PC1 ', summary(pca3)$importance[2, 1] * 100, '%') ) +
  ylab(paste0('PC2 ', summary(pca3)$importance[2, 2] * 100, '%') )+ ggtitle('PCA metrics - Public Datasets - color by Tissue type')

pp <- ggplotly(p)
saveWidget(pp, "PCA_PublicDatasetsMetricsScaled.html")

#--------#

## Touchstone and Public together! ##
library(Polychrome)
touch$Origin <- 'Touchstone'
public$Origin <- 'Public'
both <- rbind(touch[, c('Origin','platform', 'TissueType', which_metrics)], public[, c('Origin','platform', 'TissueType', which_metrics)])
both$Origin <- factor(both$Origin, levels = c("Touchstone", "Public"))


pca4 <- prcomp(scale(both[, which_metrics]))
pca <- as.data.frame(pca4$x)
#rownames(pca) <- make.names(public$dataset, unique = T)
pca$Platform <- both$platform
pca$TissueType <- both$TissueType
pca$Origin <- both$Origin

#pca$Institution <- touch$Institution
p <- ggplot(pca, aes(PC1, PC2, label = rownames(pca), colour = TissueType, shape = Platform)) + geom_point(size =5) +
  scale_colour_manual(values = as.character(kelly.colors(14))) + xlab(paste0('PC1 ', summary(pca3)$importance[2, 1] * 100, '%') ) +
  ylab(paste0('PC2 ', summary(pca3)$importance[2, 2] * 100, '%') )+ ggtitle('PCA metrics - Touchstone and Public datasets') +
  theme_bw() +  theme(text=element_text(size=20))

p2 <- ggplot(pca, aes(PC1, PC2, label = rownames(pca), colour = Origin, shape = Platform)) + geom_point(size =5) +
   xlab(paste0('PC1 ', summary(pca3)$importance[2, 1] * 100, '%') ) + scale_colour_manual(values = c('#87ceeb','#f08080' )) +
  ylab(paste0('PC2 ', summary(pca3)$importance[2, 2] * 100, '%') )+ ggtitle('PCA metrics - Touchstone and Public datasets') +
  theme_bw() +  theme(text=element_text(size=20))
  

#pdf('plots/AllSamples-AllMetrics-PCA_TissueType.pdf', 10,8)
png('plots/AllSamples-AllMetrics-PCA_TissueType.png', width = 10,height = 8, res = 300, units = 'in')
p
dev.off()

#pdf('plots/AllSamples-AllMetrics-PCA_Origin.pdf', 10,8)
png('plots/AllSamples-AllMetrics-PCA_Origin.png', width =  10 , height =  8, res = 300, units ='in')
p2
dev.off()





pp <- subplot(ggplotly(p) , ggplotly(p2) , shareX = T, shareY = T, titleX = T,
              titleY = T)
saveWidget(pp, "PCA_Touchstone_Public_MetricsScaled.html")




## reproducibility 
table(touch$Institution, touch$TissueType)

## Plot every metric vs all others by tissue type + regression line
ggplot(touch, aes(Complexity_CommonGenesNorm, TxPerCellNorm, col = platform)) + geom_point(size = 4) +  geom_smooth(method = "lm", se = F) + facet_wrap(~ TissueType)


metrics <- colnames(touch)[c(6:17, 21, 24:25)]
plots <- list()

save_plots_pdf <- function(plots, start_index, end_index) {
  pdf(file_name, width = 25, height = 15)  # Set the size of the PDF
  do.call(grid.arrange, plots[start_index:end_index]) # Print each plot into the PDF
  dev.off()  # Close the PDF device
}

make_pair_plots <- function(df, batch_size , fileName = 'plots/ScatterPlot_', metrics) {
  for (i in metrics) {
    for (j in metrics) {
      if (i != j) {  # avoids plotting a metric against itself
        plot_title <- paste("Scatter plot of", i, "vs", j)
        p <- ggplot(df, aes_string(x = i, y = j, col = "platform")) +
          geom_point(size =3) +
          geom_smooth(method = "lm", se = FALSE) +
          facet_wrap(~TissueType) +
          ggtitle(plot_title) + scale_color_manual(values = c(brewer.pal(2, 'Dark2')))
        # Save each plot with a unique name based on the metric names
        plots[[paste(i, j, sep = "_vs_")]] <- p
      }
    }
  }
  
  
  fnames <- gsub("_.*$", "", names(plots)[seq(1,156, by = 12)])
  
  
  batch_size <- batch_size # Define how many plots per PDF
  num_plots <- length(plots)
  num_files <- ceiling(num_plots / batch_size)  # Calculate the number of files needed
  
  for (file_num in 1:num_files) {
    start_index <- (file_num - 1) * batch_size + 1
    end_index <- min(file_num * batch_size, num_plots)
    file_name <- paste0(fileName, fnames[file_num], "_vs_All.pdf")  # Create a filename for each batch
    
    save_plots_pdf(plots, start_index, end_index, file_name)
  }
}

make_pair_plots(public, batch_size = 14, metrics = metrics, fileName = 'plots/ScatterPlot_PublicDatasets_')

## Reproducibility 
## looking at only PR samples
repro <- touch[touch$TissueType %in% 'PR',]
repro <- repro[-1,] # remove sample from EAP
repro <- repro[, -c(3:5)]
repro_scaled <- repro
repro_scaled[, -c(1:2,15:16,23:24)] <- scale(repro_scaled[, -c(1:2,15:16,23:24)] )


## Plot metrics public datasets
metrics <- colnames(public)[c(6:20, 22)]
plots <- list()
for(i in metrics) {
  p <- ggplot(repro, aes_string(y = i, x = "Institution", col = "Institution"))  + geom_point(size = 3) + facet_wrap(~ platform)
  plots[[i]] <- p
}

do.call(grid.arrange, c(plots, ncol = 5)) 


#calculate standard deviation and variance 
xenium_repro <- lapply(repro_scaled[repro_scaled$platform %in% 'Xenium', metrics], sd) %>% unlist 
xenium_repro <- data.frame(Metric = names(xenium_repro), SD = xenium_repro)

p1 <- ggplot(xenium_repro, aes(SD, reorder(Metric, SD))) + geom_bar(stat='identity') + xlab('Standard deviation - between Institutions') +
  ylab('Calculated metric')

xenium_repro_var <- lapply(repro[repro$platform %in% 'Xenium', metrics], var) %>% unlist 
xenium_repro_var <- data.frame(Metric = names(xenium_repro_var), Var = xenium_repro_var)

p2 <- ggplot(xenium_repro_var, aes(Var, reorder(Metric, Var))) + geom_bar(stat='identity') + xlab('Variance - between Institutions') +
  ylab('Calculated metric') 

plots[[17]] <- p1
plots[[18]] <- p2

do.call(grid.arrange, c(plots, ncol = 5)) 

### Show correlation of each gene per pair of same sample across institution
xenium_repro <- repro$sample_id[repro$platform == 'Xenium']

# get path of expression matrix files
paths <- touch$expMat[touch$sample_id %in% xenium_repro]
paths <- paths[-4] ## removing XRCP sample

exp_list <- list()
for( i in paths) {
  exp_list[i] <- getMeanExpression(expMat = i, platform = 'Xenium')
}
names(exp_list) <- xenium_repro

exp_xen_repro <- as.data.frame(exp_list)

colnames(exp_xen_repro) <- c('Xenium_STJ', 'Xenium_UOA_R1', 'Xenium_UOA_R2', 'Xenium_WCM')

#exp_xen_repro <- exp_xen_repro[, -4]

library(ggplot2)
library(gridExtra) 

# Function to create scatter plots with Pearson correlation and regression line
plot_pairs <- function(df, color_line = 'blue') {
  column_names <- colnames(df)
  plot_list <- list()  # To store ggplot objects
  
  # Generate all combinations of column pairs
  combinations <- combn(column_names, 2, simplify = FALSE)
  
  # Loop through each combination and create a scatter plot
  for (comb in combinations) {
    # Calculate Pearson correlation
    correlation <- cor(df[[comb[1]]], df[[comb[2]]], method = "pearson", use = "complete.obs")
    
    # Create the scatter plot
    plot <- ggplot(df, aes_string(x = comb[1], y = comb[2])) +
      geom_point(alpha = 0.6) +  # Adjust point transparency
      geom_smooth(method = "lm", color = color_line, se = FALSE, size = 0.5) +  # Thinner line
      labs(title = paste(comb[1], "vs", comb[2],
                         "\nPearson Correlation: ", sprintf("%.2f", correlation)),
           x = comb[1], y = comb[2]) +
      theme_minimal()
    
    # Store the plot with correlation in the title
    plot_list[[paste(comb[1], comb[2], sep = "_vs_")]] <- plot
  }
  
  return(plot_list)
}

plots <- plot_pairs(exp_xen_repro)

do.call(grid.arrange, c(plots, ncol = 3))

## For CosMx
cosmx_repro <- repro$sample_id[repro$platform == 'CosMx']

paths <- touch$expMat[touch$sample_id %in% cosmx_repro]

exp_list <- list()
for( i in paths) {
  exp_list[i] <- getMeanExpression(expMat = i, platform = 'CosMx')
}
names(exp_list) <- cosmx_repro

exp_csmx_repro <- as.data.frame(exp_list)

colnames(exp_csmx_repro) <- c('CosMx_STJ', 'CosMx_UOA_R1', 'CosMx_UOA_R2', 'CosMx_WCM')

plots2 <- plot_pairs(exp_csmx_repro, color_line = 'orange')
do.call(grid.arrange, c(plots2, ncol = 3))

do.call(grid.arrange, c(plots, plots2, ncol = 3))









