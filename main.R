#Libraries needed for the script
library(Seurat)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(rstudioapi)

#Global variables to use
project_path = "VLab_RNA_Seq"
github_path = getwd()
data_path = 'data'
setwd(file.path(github_path,project_path))

#Get data from file location
mtx_dirs = dir(data_path, recursive = F, full.names = F)

#Create Seurat objects
for (mtx in mtx_dirs){
  cts = ReadMtx(mtx = file.path(data_path,mtx,paste0(mtx,'_matrix.mtx.gz')),
                features = file.path(data_path,mtx,paste0(mtx,'_features.tsv.gz')),
                cells = file.path(data_path,mtx,paste0(mtx,'_barcodes.tsv.gz')))
  
  assign(mtx, CreateSeuratObject(counts=cts))
}

#Merge datasets 
merged_seurat = merge(GSM6736410, c(GSM6736411, GSM6736412, GSM6736413),
                      add.cells.ids = ls()[4:7],
                      project = 'GSM')
saveRDS(merged_seurat, file = "GSM_merged_seurat.rds")

View(merged_seurat@meta.data$nCount_RNA)
