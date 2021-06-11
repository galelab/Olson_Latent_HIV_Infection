#Description: Analysis for CD4 T Cells Single cell (GSE176386)
#Author: Leanne S. Whitmmore 
#

#######################
# Libraries
#######################
library(dplyr)
library(Seurat) #v3.2.3
library(ggplot2)
library(stringr)
library(data.table)
library(VennDiagram)
library(monocle3)
library(garnett)
library(org.Hs.eg.db)
library(pheatmap)
source("SC_DE_Fig.r")

##########################
# Functions
##########################
generate_folder <- function(foldername) {
    workDir <- getwd()
    subDir <- foldername
    results_path <- file.path(workDir, subDir)
    if (file.exists(subDir)) {
    } else {
        dir.create(results_path)
    }
    return(results_path)
}

theme_Publication <- function(base_size = 14, base_family = "arial") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size = base_size, base_family = base_family)
    + theme(
            plot.title = element_text(
                face = "bold",
                size = rel(1.2), hjust = 0.5
            ),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold", size = rel(1)),
            axis.title.y = element_text(angle = 90, vjust = 2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour = "#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size = unit(0.2, "cm"),
            legend.margin = unit(0.1, "cm"),
            legend.title = element_text(size = 10, face = "bold"),
            plot.margin = unit(c(10, 5, 5, 5), "mm"),
            strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
            strip.text = element_text(face = "bold")
        ))
}

scale_fill_Publication <- function(...) {
    library(scales)
    discrete_scale("fill", "Publication", manual_pal(values = c("#386cb0", "#fdb462", "#7fc97f", "#ef3b2c", "#662506", "#a6cee3", "#fb9a99", "#984ea3", "#ffff33")), ...)
}

scale_colour_Publication <- function(...) {
    library(scales)
    discrete_scale("colour", "Publication", manual_pal(values = c("#386cb0", "#fdb462", "#7fc97f", "#ef3b2c", "#662506", "#a6cee3", "#fb9a99", "#984ea3", "#ffff33")), ...)
}

run_normalization <- function(immune.combined, SD = TRUE, SCTRANSFORM = FALSE) {
    print("STATUS: normalizing data")
    all.genes <- rownames(immune.combined)
    if (SD == TRUE) {
        immune.combined <- ScaleData(immune.combined,
            vars.to.regress = "percent.mt",
            features = all.genes, verbose = FALSE
        )
    } else if (SCTRANSFORM == TRUE) {
        immune.combined <- SCTransform(immune.combined,
            vars.to.regress = "percent.mt",
            verbose = FALSE
        )
    }
    return(immune.combined)
}

generate_monocle_cds <- function(sample_path) {
    SAMPLE <- load_cellranger_data(sample_path, umi_cutoff = 200)

    SAMPLE.counts <- SAMPLE@assays@data@listData$counts
    gene_meta_data <- rowData(SAMPLE)
    cell_meta_data <- colData(SAMPLE)

    cds <- new_cell_data_set(SAMPLE.counts,
        cell_metadata = cell_meta_data,
        gene_metadata = gene_meta_data
    )
    return(cds)
}

feature_reduction <- function(immune.combined, result_folder) {
    print("STATUS: performing PCA and UMAP")

    # results_path <- generate_folder(result_folder)
    # unlink(paste0(result_folder, "/*"))

    immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
    immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)

    return(immune.combined)
}


initial_SC_Seurat <- function(data_dir, sampleID, mock=FALSE) {
    pbmc.data <- Read10X(data.dir = data_dir)

    pbmc <- CreateSeuratObject(counts = pbmc.data,
                               project = sampleID,
                               min.cells = 0,
                               min.features = 200)
    pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
    VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    ggsave(paste0(results_path, sampleID, "_VnPlotMt.png"), dpi = 300)

    pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
    VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    ggsave(paste0(results_path, sampleID, "_VnPlotMt_Filtered.png"), dpi = 300)
    z <- pbmc["hiv-positive", ]
    x1 <- z@meta.data[z@meta.data$nFeature_RNA != 0, ]
    HIVdata <- c()
    for (cell in colnames(pbmc)) {
        if (cell %in% rownames(x1)) {
            HIVdata <- c(HIVdata, "HIV")
        } else {
            HIVdata <- c(HIVdata, "noHIV")
        }
    }
    pbmc$HIV <- HIVdata

    if (isTRUE(mock)) {
        print("STATUS: removing HIV infected cells from analysis")
        cells <- pbmc@meta.data[pbmc@meta.data$HIV == "HIV", ]
        print(paste0("STATUS: number of HIV cells being removed ", length(rownames(cells))))
        if (length(cells) > 0) {
            nohivcells <- setdiff(rownames(pbmc@meta.data), rownames(cells))
            pbmc <- subset(pbmc, cells = nohivcells)
        }
    }

    cds <- generate_monocle_cds(str_remove(data_dir, "outs/filtered_feature_bc_matrix"))

    cds <- cds[, colnames(pbmc)]


    pbmc <- NormalizeData(pbmc,
        normalization.method = "LogNormalize",
        scale.factor = 10000
    )


    # Find features (genes) that are highly variable from cell-to-cell
    pbmc <- FindVariableFeatures(pbmc,
        selection.method = "vst",
        nfeatures = 2000
    )

    # Identify the 10 most highly variable genes
    top10 <- head(VariableFeatures(pbmc), 10)

    # plot variable features with and without labels
    plot1 <- VariableFeaturePlot(pbmc)
    LabelPoints(plot = plot1, points = top10, repel = TRUE)
    ggsave(paste0(results_path, sampleID, "_HighlyVariableGenes.png"), dpi = 300)

    pbmc <- ScaleData(pbmc)
    pbmc <- RunPCA(pbmc, features=VariableFeatures(pbmc))
    pbmc <- RunUMAP(pbmc, dims = 1:10)

    VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
    ggsave(paste0(results_path, sampleID, "_PCA_Loadings.png"), dpi = 300)

    z <- pbmc["hiv-positive", ]
    x1 <- z@meta.data[z@meta.data$nFeature_RNA != 0, ]
    if (length(rownames(x1))>0) {
        print(paste0("STATUS: Generating HIV figure for ", sampleID))
        pbmchiv <- subset(pbmc, cells = rownames(x1))

        VlnPlot(pbmchiv, features = c("hiv-positive"))
        ggsave(paste0(results_path, sampleID, "_hiv_VnPlotMt_onlyhivcells.png"), dpi = 300)
    }

    return(list("seurat"=pbmc, "monocle"=cds))
}

get_DE_between_conditions <- function(ident_1, ident_2, compare,
                                      immune.combined,
                                      result_folder, fontsize = 9, height=8, foldchange = 0.26,
                                      pval = 0.05, percent_cells = NULL, dotscale=3.5) {
    contrasting <- colorRampPalette(rev(c("firebrick4", "firebrick1", "white", "steelblue1", "steelblue4")))(100)
    print("STATUS: getting DEs...")
    ### Parameters for FindMarkers
    #### test.use: mast (default is wilcox)
    #### min.pct: 0.1(default: 0.1) filter out genes (features) that are detected at less than 10 percent frequency in cells in ident_1 or ident_2
    #### logfc: 0 (default is .25) logfc must be higher than 0 (I set this at 0 because I filter this out at a later stage - this was primarily done to
    ####                          understand how the tool (Seurat) works
    #### min.cells.feature: 3 (default is 3) minimum number of cells expressing the feature in at least one of the two groups (similar to min.pct)
    #### min.cells.group: 3 (defualt is 3) minimum number of cells in the group
    #### max.cells.per.ident: Inf (default Inf-means no down sampling) Down sample each identity class (cluster of cells) to a max number
    #### min.dif.pct: Inf (default Inf) Only test genes that show minimum difference in the fraction of detection between the two identities (cluster)
    if (file.exists(paste0(result_folder, "DEgenes_full_",ident_1, "_", ident_2, "_", compare, ".txt"))) {
        message("STATUS: DE has already been done filtering from full file ")
        DEgenes <- read.table(paste0(result_folder, "DEgenes_full_", ident_1, "_", ident_2, "_", compare, ".txt"), header = TRUE, row.names = 1)
    } else {
        DEgenes <- FindMarkers(immune.combined,
            ident.1 = ident_1,
            ident.2 = ident_2,
            test.use = "MAST",
            logfc.threshold = 0
        )
        write.table(data.frame(DEgenes),
            paste0(
                result_folder,
                "DEgenes_full_",
                ident_1,
                "_", ident_2, "_", compare,
                ".txt"
            ),sep = "\t", quote = FALSE)
    }

    DEgenes[["gene_name"]] <- rownames(DEgenes)
    DE_sig_final <- DEgenes %>%
        filter(avg_logFC >= foldchange | avg_logFC <= -foldchange) %>%
        dplyr::select(gene_name, p_val, avg_logFC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
    DE_sig_final <- DE_sig_final %>%
        filter(p_val_adj <= pval) %>%
        dplyr::select(gene_name, p_val, avg_logFC, pct.1, pct.2, p_val_adj)
    
    if (!is.null(percent_cells)) {
 
        DE_sig_final <- DE_sig_final %>%
            filter(pct.1 >= percent_cells | pct.2 >= percent_cells) %>%
            dplyr::select(gene_name, p_val, avg_logFC, pct.1, pct.2, p_val_adj)
    }
    rownames(DE_sig_final) <- DE_sig_final$gene_name

    DE_sig_final$gene_name <- NULL
    if (is.null(percent_cells)) {
        write.table(data.frame(DE_sig_final),
            paste0(
                result_folder,
                "DEgenes_sig_",
                ident_1, "_", ident_2, "_", compare,
                ".txt"
            ),
            sep = "\t", quote = FALSE
        )
    } else {
        write.table(data.frame(DE_sig_final),
            paste0(
                result_folder,
                "DEgenes_sig_",
                ident_1, "_", ident_2, "_", compare,
                "_",  percent_cells, ".txt"
            ),
            sep = "\t", quote = FALSE
        )
    }
    if (length(rownames(DE_sig_final))>0) {
        DotPlot(immune.combined, 
            idents = c(ident_1, ident_2),
            features = rownames(DE_sig_final),
            assay="RNA",
            col.min=0, col.max=2.5,
            cols = c("black", "#eb3306"),
            dot.scale = dotscale) + coord_flip() +
            theme(axis.text.y = element_text(size = fontsize))
        if (is.null(percent_cells)) {
            ggsave(paste0(
                result_folder, "DEgenes_sig_", ident_1, "_",
                ident_2, "_", compare, ".png"
            ), height = 8, width = 6, units = "in", dpi = 350)
        } else {
            ggsave(paste0(
                result_folder, "DEgenes_sig_", ident_1, "_",
                ident_2, "_", compare,"_", percent_cells, ".png"
            ), height = 8, width = 6, units = "in", dpi = 350)
        }
    }
    if (length(rownames(DE_sig_final))>1){
        message("STATUS: Making BarPlot")
        pl <- SC_DE_barplot(DE_sig_final, str_remove(ident_1, "CD4-T-cell_"), str_remove(ident_2, "CD4-T-cell_"),fontsize=fontsize)
        if (is.null(percent_cells)) {
            ggsave(paste0(
                result_folder, "DEgenes_sig_", ident_1, "_",
                ident_2, "_", compare, "barPlot.png"
            ),
            width = 5, height = height, units = "in", dpi = 300
            )
        } else {
            ggsave(paste0(
                result_folder, "DEgenes_sig_", ident_1, "_",
                ident_2, "_", compare, "_", percent_cells, "barPlot.png"
            ),
            width = 5, height = height, units = "in", dpi = 300
            )
        }
    }
    return(DE_sig_final)
}


get_percentages <- function(tabletemp_orig) {
    tabletemp <- data.frame(table(tabletemp_orig))
    tabletemp$Freq <- tabletemp$Freq / length(tabletemp_orig)
    return(tabletemp)
}

get_virus_percentages <- function(tabletemp_orig, sample) {
    hiv_final <- data.frame()
    for (cell in unique(tabletemp_orig)) {
        x <- which(tabletemp_orig == cell)
        total_per <- length(x) / length(tabletemp_orig)
        z <- immune.combined_sd$HIV[names(x)]
        hiv <- as.data.frame(table(z))
        hiv$Freq <- (hiv$Freq / length(x)) * total_per
        hiv$cell <- rep(cell, length(hiv$Freq))
        hiv$sample <- rep(sample, length(hiv$Freq))
        hiv_final <- rbind(hiv_final, hiv)
    }
    return(hiv_final)
}

get_virus_percentages_total <- function(tabletemp_orig, sample) {
    hiv_final <- data.frame()
    for (cell in unique(tabletemp_orig)) {
        x <- which(tabletemp_orig == cell)
        total_per <- length(x) / length(tabletemp_orig)
        z <- immune.combined_sd$HIV[names(x)]
        hiv <- as.data.frame(table(z))
        hiv$Freq <- (hiv$Freq / length(x))
        hiv$cell <- rep(cell, length(hiv$Freq))
        hiv$sample <- rep(sample, length(hiv$Freq))
        hiv_final <- rbind(hiv_final, hiv)
    }
    return(hiv_final)
}

#############################################
# Main Program - Process individual samples 
#############################################

## ----Paths to results folder and raw data
SAMPLE_PATH <- "/vol08/ngs/P51/HIV/HIV01_scLatency/WhitmoreAnalysis/7th_run/"
results_path <- "3.SC_Analysis/"
generate_folder(results_path)

## ----PROCESS INDIVIDUAL SAMPLES
G001 <- initial_SC_Seurat(paste0(SAMPLE_PATH, "/G001_HIV01_scLatency_CD4-T-cell_M_8hr_D1_Lib1/outs/filtered_feature_bc_matrix"),
    sampleID = "CD4-T-cell_M_8hr_D1", mock=TRUE
)
G002 <- initial_SC_Seurat(paste0(SAMPLE_PATH, "/G002_HIV01_scLatency_CD4-T-cell_M_IFNb_8hr_D1_Lib1/outs/filtered_feature_bc_matrix"),
    sampleID = "CD4-T-cell_M_IFNb_8hr_D1", mock=TRUE
)
G003 <- initial_SC_Seurat(paste0(SAMPLE_PATH, "/G003_HIV01_scLatency_CD4-T-cell_HIV_8hr_D1_Lib1/outs/filtered_feature_bc_matrix"),
    sampleID = "CD4-T-cell_HIV_8hr_D1"
)
G004 <- initial_SC_Seurat(paste0(SAMPLE_PATH, "/G004_HIV01_scLatency_CD4-T-cell_HIV_IFNb_8hr_D1_Lib1/outs/filtered_feature_bc_matrix"),
    sampleID = "CD4-T-cell_HIV_IFNb_8hr_D1"
)
G005 <- initial_SC_Seurat(paste0(SAMPLE_PATH, "/G005_HIV01_scLatency_CD4-T-cell_M_8hr_D2_Lib1/outs/filtered_feature_bc_matrix"),
    sampleID = "CD4-T-cell_M_8hr_D2", mock=TRUE
)
G006 <- initial_SC_Seurat(paste0(SAMPLE_PATH, "/G006_HIV01_scLatency_CD4-T-cell_M_IFNb_8hr_D2_Lib1/outs/filtered_feature_bc_matrix"),
    sampleID = "CD4-T-cell_M_IFNb_8hr_D2", mock=TRUE
)
G007 <- initial_SC_Seurat(paste0(SAMPLE_PATH, "/G007_HIV01_scLatency_CD4-T-cell_HIV_8hr_D2_Lib1/outs/filtered_feature_bc_matrix"),
    sampleID = "CD4-T-cell_HIV_8hr_D2"
)
G008 <- initial_SC_Seurat(paste0(SAMPLE_PATH, "/G008_HIV01_scLatency_CD4-T-cell_HIV_IFNb_8hr_D2_Lib1/outs/filtered_feature_bc_matrix"),
    sampleID = "CD4-T-cell_HIV_IFNb_8hr_D2"
)

## -----COMBINED ALL SAMPLES TOGETHER FOR ANALSYIS
sample.list <- list(
    "CD4-T-cell_M_8hr_D1" = G001$seurat, "CD4-T-cell_M_IFNb_8hr_D1" = G002$seurat,
    "CD4-T-cell_HIV_8hr_D1" = G003$seurat, "CD4-T-cell_HIV_IFNb_8hr_D1" = G004$seurat,
    "CD4-T-cell_M_8hr_D2" = G005$seurat, "CD4-T-cell_M_IFNb_8hr_D2" = G006$seurat,
    "CD4-T-cell_HIV_8hr_D2" = G007$seurat, "CD4-T-cell_HIV_IFNb_8hr_D2" = G008$seurat
)
saveRDS(sample.list, "3.sample.list.rds")

## -----INTEGRATE AND NORMALIZE ALL SAMPLES TOGETHER
generate_folder("3.SC_integrated_analysis/")
immune.anchors <- FindIntegrationAnchors(object.list = sample.list, dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)
saveRDS(immune.combined, "3.Immune.combinedobject.rds")
# immune.combined <- readRDS("3.Immune.combinedobject.rds")

## -----Normalize all samples
DefaultAssay(immune.combined) <- "integrated"
immune.combined_sd <- run_normalization(immune.combined)

## -----Normalize RNA assay (TAKES A LONG TIME)
DefaultAssay(immune.combined_sd) <- "RNA"
immune.combined_sd <- run_normalization(immune.combined_sd)
saveRDS(immune.combined_sd, "3.Immune.combinedobject_norm.rds")
# immune.combined_sd <- readRDS("3.Immune.combinedobject_norm.rds")


## -----Feature reduction
DefaultAssay(immune.combined_sd) <- "integrated"
immune.combined_sd <- feature_reduction(immune.combined_sd, "./3.SC_integrated_analysis/")
saveRDS(immune.combined_sd, "3.Immune.combinedobject_norm_fr.rds")
immune.combined_sd <- readRDS("3.Immune.combinedobject_norm_fr.rds")

###############################
# Classify cell type
############################### 
generate_folder("3.SC_Celltype_classification/")
batch.name <- c("G001", "G002", "G003", "G004", "G005", "G006", "G007", "G008")
cds.list <- list(
    G001$monocle, G002$monocle, G003$monocle, G004$monocle,
    G005$monocle, G006$monocle, G007$monocle, G008$monocle
)
names(cds.list) <- batch.name
big_cds <- combine_cds(cds.list)

## ------Align/Normalize CDS object
big_cds <- preprocess_cds(big_cds, num_dim = 100, preprocess_method = 100)
big_cds <- align_cds(big_cds, alignment_group = "sample")

big_cds <- reduce_dimension(big_cds,
    preprocess_method = "PCA",
    reduction_method = "UMAP"
)
saveRDS(big_cds, "3.SC_Celltype_classification/2.big_cds.rds")
# big_cds <- readRDS("3.SC_Celltype_classification/2.big_cds.rds")

## ------Build training data set using markers
marker_file_path <- "hs_markers_Tcells.txt"
npc_classifier <- train_cell_classifier(
    cds = big_cds,
    marker_file = marker_file_path,
    db = org.Hs.eg.db,
    cds_gene_id_type = "ENSEMBL",
    num_unknown = 500,
    cores=16,
    marker_file_gene_id_type = "SYMBOL"
)
saveRDS(npc_classifier, "3.SC_Celltype_classification/3.npc_classifier.rds")
npc_classifier <- readRDS("3.SC_Celltype_classification/3.npc_classifier.rds")

## ------Classify cells
big_cds1 <- classify_cells(big_cds, npc_classifier,
    db = org.Hs.eg.db,
    cluster_extend = TRUE,
    cds_gene_id_type = "ENSEMBL"
)
saveRDS(big_cds1, "3.SC_Celltype_classification/3.big_cds_clusterext.rds")
big_cds1 <- readRDS("3.SC_Celltype_classification/3.big_cds_clusterext.rds")
marker_check <- check_markers(big_cds1, marker_file_path,
    db = org.Hs.eg.db,
    cds_gene_id_type = "ENSEMBL",
    marker_file_gene_id_type = "SYMBOL"
)
plot_markers(marker_check)
ggsave("3.SC_Celltype_classification/3.Marker_Fig.png", dpi = 500)

## ------add cell type classifications to integrated seurat object
g_cluster_ext <- big_cds1$cluster_ext_type
immune.combined_sd[["garnett_cluster_extend"]] <- g_cluster_ext
Idents(immune.combined_sd) <- immune.combined_sd$garnett_cluster_extend
DimPlot(immune.combined_sd, reduction = "umap")
ggsave("3.SC_Celltype_classification/3.garnett_clusters.png", width = 6.5, height = 4, dpi = 500)

## ------find clusters in integrated seurat object
immune.combined_sd <- FindNeighbors(immune.combined_sd, dims = 1:10)
immune.combined_sd <- FindClusters(immune.combined_sd) # Default values ued
Idents(immune.combined_sd) <- immune.combined_sd$seurat_clusters
DimPlot(immune.combined_sd, reduction = "umap")
ggsave("3.SC_Celltype_classification/3.seruat_clusters.png", width = 6.5, height = 4, dpi = 500)

## ------for cells with the unknown classification assign them the cell type that is most prevelant in a seurat cluster
immune.combined_sd[["garnett_cluster_extend_lw"]] <- immune.combined_sd$garnett_cluster_extend
for (cluster in unique(immune.combined_sd$seurat_clusters)) {
    celltype_count <- c()
    clust_list <- which(immune.combined_sd$seurat_clusters == cluster)
    tempdf <- immune.combined_sd@meta.data[clust_list, ]
    celltypes <- unique(tempdf$garnett_cluster_extend)
    for (celltype in celltypes) {
        if (celltype != "Unknown") {
            celltype_count[celltype] <- length(which(tempdf$garnett_cluster_extend == celltype))
        }
    }
    maxvalue <- max(celltype_count)
    celltypefinal <- names(which(celltype_count == maxvalue))
    print (paste0("STATUS: max cell count in cluster ", cluster, " is ", celltypefinal))
    tempdfunknown <- immune.combined_sd@meta.data[immune.combined_sd@meta.data$garnett_cluster_extend == "Unknown", ]
    tempdfunknown <- tempdfunknown[tempdfunknown$seurat_clusters==cluster, ]
    print(paste0("STATUS number of unknowns in cluster ", cluster, " is ", dim(tempdfunknown)[1]))
    for (row in rownames(tempdfunknown)) {
        immune.combined_sd@meta.data[row, "garnett_cluster_extend_lw"] <- celltypefinal
    }
}
Idents(immune.combined_sd) <- immune.combined_sd$garnett_cluster_extend_lw
DimPlot(immune.combined_sd, reduction = "umap")
ggsave("3.SC_Celltype_classification/2.garnett_clusters_lw.png", width = 6.5, height = 4, dpi = 500)
saveRDS(immune.combined_sd, "3.Immune.combinedobject_norm_fr_ct.rds")
# immune.combined_sd <- readRDS("3.Immune.combinedobject_norm_fr_ct.rds")


###############################
#Perform DE analysis
###############################
## -----2ndAnalysis-Get lists of IFN stimulated genes
generate_folder("3.SC_integrated_analysis/2nd_analysis/")

## -- Donors combined DE
immune.combined_sd@meta.data$sample <- str_remove_all(immune.combined_sd@meta.data$orig.ident, "_D\\d+")
Idents(immune.combined_sd) <- immune.combined_sd@meta.data$sample

DEHIVIFNb_HIV <- get_DE_between_conditions("CD4-T-cell_HIV_IFNb_8hr", "CD4-T-cell_HIV_8hr", "HIVIFNb-HIV", immune.combined_sd, "./3.SC_integrated_analysis/2nd_analysis/", fontsize = 7)
DEIFNb_mock <- get_DE_between_conditions("CD4-T-cell_M_IFNb_8hr", "CD4-T-cell_M_8hr", "IFNb-mock", immune.combined_sd, "./3.SC_integrated_analysis/2nd_analysis/", fontsize = 5)


## -----3rdAnalysis-HIV Analysis
generate_folder("./3.SC_integrated_analysis/3rd_Analysis/")
Idents(immune.combined_sd) <- immune.combined_sd@meta.data$orig.ident
immune.combined_sd$sample_HIV <- paste(immune.combined_sd@meta.data$orig.ident, immune.combined_sd@meta.data$HIV, sep="_")

cells <- immune.combined_sd@meta.data[immune.combined_sd@meta.data == "CD4-T-cell_HIV_8hr_D1" | 
                                         immune.combined_sd@meta.data == "CD4-T-cell_HIV_8hr_D2",]
immune_subset1 <- immune.combined_sd[,rownames(cells)]
Idents(immune_subset1) <- immune_subset1$HIV

## --DE HIV only pooled
DEgenesHIV <- get_DE_between_conditions("HIV", "noHIV", "HIVonly", immune_subset1, "./3.SC_integrated_analysis/3rd_Analysis/", fontsize=5)
pl <- SC_DE_barplot(DEgenesHIV, "HIV", "noHIV", fontsize = 10, horizontal = TRUE)
ggsave("./3.SC_integrated_analysis/3rd_Analysis/DEgenes_sig_HIV_noHIV_HIVonlybarPlotH.png", width = 13, height = 4, units = "in", dpi = 300)

## --DE HIV INFb pooled
cells <- immune.combined_sd@meta.data[
    immune.combined_sd@meta.data == "CD4-T-cell_HIV_IFNb_8hr_D1" |
        immune.combined_sd@meta.data == "CD4-T-cell_HIV_IFNb_8hr_D2",
]
immune_subset2 <- immune.combined_sd[, rownames(cells)]

Idents(immune_subset2) <- immune_subset2$HIV
DEgenesINFb <- get_DE_between_conditions("HIV", "noHIV", "HIV_IFNb", immune_subset2, "./3.SC_integrated_analysis/3rd_Analysis/", fontsize=5)

##################################################################
#Splits HIV into high and low HIV expression and does DE analysis
##################################################################
results_path <- "3.DE_HIV_lowvshighHIV/"
generate_folder(results_path)

cells1 <- immune.combined_sd@meta.data[immune.combined_sd@meta.data == "CD4-T-cell_HIV_8hr_D1" |
    immune.combined_sd@meta.data == "CD4-T-cell_HIV_8hr_D2" | immune.combined_sd@meta.data == "CD4-T-cell_M_8hr_D1" |
    immune.combined_sd@meta.data == "CD4-T-cell_M_8hr_D2", ]
cells_HIV <- immune.combined_sd@meta.data[immune.combined_sd@meta.data == "CD4-T-cell_HIV_8hr_D1" |
    immune.combined_sd@meta.data == "CD4-T-cell_HIV_8hr_D2", ]
cells_mock <- immune.combined_sd@meta.data[immune.combined_sd@meta.data == "CD4-T-cell_M_8hr_D1" |
    immune.combined_sd@meta.data == "CD4-T-cell_M_8hr_D2", ]
cells2 <- immune.combined_sd@meta.data[immune.combined_sd@meta.data == "CD4-T-cell_HIV_IFNb_8hr_D1" |
    immune.combined_sd@meta.data == "CD4-T-cell_HIV_IFNb_8hr_D2" | immune.combined_sd@meta.data == "CD4-T-cell_M_IFNb_8hr_D1" |
    immune.combined_sd@meta.data == "CD4-T-cell_M_IFNb_8hr_D2", ]

immune_subset1 <- immune.combined_sd[, rownames(cells1)]
immune_subset2 <- immune.combined_sd[, rownames(cells2)]

# -- Set threshold for low and high HIV expression used 50% 
m <- GetAssayData(object = immune_subset1, assay = "RNA", slot = "data")
hiv <- m["hiv-positive", ]
nohiv <- hiv[hiv == 0]
hivexprs <- hiv[hiv > 0]
x <- quantile(hivexprs)
classificationhiv <- c()
for (c in hivexprs) {
    if (c > x[[3]]) {
        classificationhiv <- c(classificationhiv, "high_hiv")
    } else {
        classificationhiv <- c(classificationhiv, "low_hiv")
    }
}
names(classificationhiv) <- names(hivexprs)
nohivclass <- rep("no_hiv", length(names(nohiv)))
names(nohivclass) <- names(nohiv)
c <- as.data.frame(classificationhiv)
l <- as.data.frame(nohivclass)

colnames(c)[1] <- "classification"
colnames(l)[1] <- "classification"

final1 <- rbind(c, l)

finaltmp <- final1[colnames(immune_subset1), ]
names(finaltmp) <- colnames(immune_subset1)
finaltmp <- as.data.frame(finaltmp)
if (all.equal(colnames(immune_subset1), rownames(finaltmp))!=TRUE) {
    message("WARNING: CELLS ARE NOT IN THE SAME ORDER RESULTS WILL BE WRONG")
    # final1 <- final1[match()]
}
immune_subset1$HIV_subclass <- as.character(finaltmp$finaltmp)

Idents(immune_subset1) <- immune_subset1@meta.data$HIV_subclass
immune_subset1@meta.data$sample <- str_remove_all(immune_subset1@meta.data$orig.ident, "_D\\d+")
immune_subset1@meta.data$HIVsubclass_sample <- paste(immune_subset1@meta.data$HIV_subclass, immune_subset1@meta.data$sample, sep="_")
Idents(immune_subset1) <- immune_subset1$HIVsubclass_sample

DEHighHIV <- get_DE_between_conditions("high_hiv_CD4-T-cell_HIV_8hr", "no_hiv_CD4-T-cell_M_8hr", "high_hiv", immune_subset1, results_path, fontsize = 8)
DEHlowHIV <- get_DE_between_conditions("low_hiv_CD4-T-cell_HIV_8hr", "no_hiv_CD4-T-cell_M_8hr", "low_hiv", immune_subset1, results_path, fontsize = 12)

DEHighHIV_orig <- DEHighHIV
DELowHIV_orig <- DEHlowHIV
DELowHIV_orig <- DELowHIV_orig[order(DELowHIV_orig$avg_logFC), ]
DEHighHIV_orig <- DEHighHIV_orig[order(-DEHighHIV_orig$avg_logFC), ]


m2 <- GetAssayData(object = immune_subset2, assay = "RNA", slot = "data")
hiv <- m2["hiv-positive", ]
nohiv <- hiv[hiv == 0]
hivexprs <- hiv[hiv > 0]
classificationhiv <- c()
for (c in hivexprs) {
    if (c > x[[3]]) {
        classificationhiv <- c(classificationhiv, "high_hiv")
    } else {
        classificationhiv <- c(classificationhiv, "low_hiv")
    }
}
names(classificationhiv) <- names(hivexprs)
nohivclass <- rep("no_hiv", length(names(nohiv)))
names(nohivclass) <- names(nohiv)
c <- as.data.frame(classificationhiv)
l <- as.data.frame(nohivclass)

colnames(c)[1] <- "classification"
colnames(l)[1] <- "classification"

final1 <- rbind(c, l)

finaltmp <- final1[colnames(immune_subset2), ]
names(finaltmp) <- colnames(immune_subset2)
finaltmp <- as.data.frame(finaltmp)
if (all.equal(colnames(immune_subset2), rownames(finaltmp)) != TRUE) {
    message("WARNING: CELLS ARE NOT IN THE SAME ORDER RESULTS WILL BE WRONG")
    # final1 <- final1[match()]
}
immune_subset2$HIV_subclass <- as.character(finaltmp$finaltmp)
immune_subset2@meta.data$sample <- str_remove_all(immune_subset2@meta.data$orig.ident, "_D\\d+")
immune_subset2@meta.data$HIVsubclass_sample <- paste(immune_subset2@meta.data$HIV_subclass, immune_subset2@meta.data$sample, sep = "_")
Idents(immune_subset2) <- immune_subset2$HIVsubclass_sample
DefaultAssay(immune_subset2) <- "RNA"


DEHighHIV <- get_DE_between_conditions("high_hiv_CD4-T-cell_HIV_IFNb_8hr", "no_hiv_CD4-T-cell_M_IFNb_8hr", "high_hiv", immune_subset2, results_path, fontsize = 8)
DElowHIV <- get_DE_between_conditions("low_hiv_CD4-T-cell_HIV_IFNb_8hr", "no_hiv_CD4-T-cell_M_IFNb_8hr", "low_hiv", immune_subset2, results_path, fontsize = 12)


##########################################
# Generate final figs for manuscript
##########################################
fig_results <- "finalfigs"
generate_folder(fig_results)
HIV_relabeled <- str_replace_all(immune.combined_sd@meta.data$HIV, "noHIV", "HIV RNA -")
HIV_relabeled <- str_replace_all(HIV_relabeled, "^HIV$", "HIV RNA +")

immune.combined_sd@meta.data$HIV_relabeled <- HIV_relabeled
cells <- immune.combined_sd@meta.data[
    immune.combined_sd@meta.data$HIV == "HIV",
]

# ---Fig 6a
DimPlot(immune.combined_sd,
    reduction = "umap",
    # group.by = "HIV_relabeled",
    pt.size = 0.1,
    cells.highlight = rownames(cells),
    sizes.highlight = 0.1,
    # cols = c("#ddd2d2", "#2d0446"),
    label = FALSE, raster = FALSE
) + scale_color_manual(labels = c("HIV RNA -", "HIV RNA +"), values = c("#ddd2d2", "black")) +
    theme(
        legend.text = element_text(size = 7)
    ) +
    labs(x = "UMAP 1", y = "UMAP 2", title = "") +
    theme(axis.title.x = element_text(size = 13)) +
    theme(axis.title.y = element_text(size = 13)) +
    theme(axis.text = element_text(size = 10)) +
    theme(legend.key.size = unit(0.2, "cm")) +
    theme(legend.position = "none")
# theme(legend.direction = "horizontal") +
# guides(colour = guide_legend(override.aes = list(size = 1)))
ggsave(file.path(fig_results, "Fig6a.png"), width = 4, height = 4, dpi = 500)
# ggsave(file.path(fig_results, "Fig6a.svg"), width = 4, height = 4, dpi = 500)
# ggsave(file.path(fig_results, "Fig6a.pdf"), width = 4, height = 4, dpi = 500)


# ---Fig 6B
library(RColorBrewer)
DefaultAssay(immune.combined_sd) <- "RNA"
FeaturePlot(immune.combined_sd,
    order = TRUE, reduction = "umap",
    shape.by = "HIV_relabeled", slot = "scale.data", features = c("hiv-positive"),
    pt.size = 0.6, min.cutoff = 0, max.cutoff = 11
) +
    # theme(legend.text = element_text(size = 7)) +
    theme(legend.text = element_text(size = 10)) +
    theme(axis.title.x = element_text(size = 13)) +
    theme(axis.title.y = element_text(size = 13)) +
    theme(axis.text = element_text(size = 10)) +
    labs(x = "UMAP 1", y = "UMAP 2", title = "") +
    guides(shape = guide_legend(override.aes = list(size = 3))) +
    scale_colour_gradientn(colours = c("#e6d8d8", "orange", "red")) # rev(brewer.pal(n = 11, name = "RdBu")))
ggsave(file.path(fig_results, "Fig6b_HIV_scaled.png"), width = 6, height = 4.5, dpi = 500)
ggsave(file.path(fig_results, "Fig6b_HIV_scaled.svg"), width = 6, height = 4.5, dpi = 500)
ggsave(file.path(fig_results, "Fig6b_HIV_scaled.pdf"), width = 6, height = 4.5, dpi = 500)


# --Fig 6E
enrichtable <- read.table("enrichment_results_wg_result1612297088.txt", stringsAsFactors = FALSE, sep = "\t", header = TRUE, row.names = 1)
enrichtable$description[3] <- "establishment of protein\nlocalization to membrane"
enrichtable$description[5] <- "protein localization\nto endoplasmic reticulum"
enrichtable$description[7] <- "ribonucleoprotein\ncomplex biogenesis"

enrichtable1 <- enrichtable[order(enrichtable$enrichmentRatio), ]
x <- as.character(enrichtable1$description)
enrichtable1$description <- factor(enrichtable1$description, levels = x)
ggplot(data = enrichtable1, aes(x = description, y = enrichmentRatio)) +
    geom_bar(stat = "identity", fill = "#868282") +
    ylim(c(0, 45)) +
    xlab("GO Biological Process Term") +
    ylab("Enrichment Ratio") +
    theme_Publication() +
    theme(axis.text = element_text(size = 10)) +
    theme(axis.title = element_text(size = 13)) +
    coord_flip()
ggsave(file.path(fig_results, "FigE.png"), width = 4.5, height = 4, units = "in", dpi = 500)
# ggsave(file.path(fig_results, "FigE.pdf"), width = 4.5, height = 4, units = "in", dpi = 500)

## -- Fig 6D
p <- DotPlot(immune_subset1, features = union(rownames(DEHighHIV_orig), rownames(DELowHIV_orig)), assay = "RNA", cols = c("blue", "#eb3306"), scale = TRUE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8))
ggsave(file.path(fig_results, "FigC_dotplot_scaled.png"), width = 16, height = 3.5, units = "in", dpi = 500)
write.csv(p$data, file.path(fig_results, "FigC_dotplot_scaled.csv"))


## --SupFig C
DefaultAssay(immune.combined_sd) <- "RNA"
## ------pull out cells for each sample
Idents(immune.combined_sd) <- immune.combined_sd$orig.ident
M_8hr_D1 <- WhichCells(object = immune.combined_sd, idents = "CD4-T-cell_M_8hr_D1")
M_8hr_D2 <- WhichCells(object = immune.combined_sd, idents = "CD4-T-cell_M_8hr_D2")
M_INFb_8hr_D1 <- WhichCells(object = immune.combined_sd, idents = "CD4-T-cell_M_IFNb_8hr_D1")
M_INFb_8hr_D2 <- WhichCells(object = immune.combined_sd, idents = "CD4-T-cell_M_IFNb_8hr_D2")
HIV_8hr_D1 <- WhichCells(object = immune.combined_sd, idents = "CD4-T-cell_HIV_8hr_D1")
HIV_8hr_D2 <- WhichCells(object = immune.combined_sd, idents = "CD4-T-cell_HIV_8hr_D2")
HIV_INFb_8hr_D1 <- WhichCells(object = immune.combined_sd, idents = "CD4-T-cell_HIV_IFNb_8hr_D1")
HIV_INFb_8hr_D2 <- WhichCells(object = immune.combined_sd, idents = "CD4-T-cell_HIV_IFNb_8hr_D2")

table1 <- immune.combined_sd$garnett_cluster_extend_lw[names(immune.combined_sd$garnett_cluster_extend_lw) %in% M_8hr_D1]
table2 <- immune.combined_sd$garnett_cluster_extend_lw[names(immune.combined_sd$garnett_cluster_extend_lw) %in% M_8hr_D2]
table3 <- immune.combined_sd$garnett_cluster_extend_lw[names(immune.combined_sd$garnett_cluster_extend_lw) %in% M_INFb_8hr_D1]
table4 <- immune.combined_sd$garnett_cluster_extend_lw[names(immune.combined_sd$garnett_cluster_extend_lw) %in% M_INFb_8hr_D2]
table5 <- immune.combined_sd$garnett_cluster_extend_lw[names(immune.combined_sd$garnett_cluster_extend_lw) %in% HIV_8hr_D1]
table6 <- immune.combined_sd$garnett_cluster_extend_lw[names(immune.combined_sd$garnett_cluster_extend_lw) %in% HIV_8hr_D2]
table7 <- immune.combined_sd$garnett_cluster_extend_lw[names(immune.combined_sd$garnett_cluster_extend_lw) %in% HIV_INFb_8hr_D1]
table8 <- immune.combined_sd$garnett_cluster_extend_lw[names(immune.combined_sd$garnett_cluster_extend_lw) %in% HIV_INFb_8hr_D2]

# ## ------ get percentages for each cell type (no zika info)
# table1p <- get_percentages(table1)
# table2p <- get_percentages(table2)
# table3p <- get_percentages(table3)
# table4p <- get_percentages(table4)
# table5p <- get_percentages(table5)
# table6p <- get_percentages(table6)
# table7p <- get_percentages(table7)
# table8p <- get_percentages(table8)

# ## ------integrate cell type information
# countstable <- merge(table1p, table2p, by.x = "tabletemp_orig", by.y = "tabletemp_orig", all = TRUE)
# colnames(countstable) <- c("cell", "CD4-T-cell_M_8hr_D1", "CD4-T-cell_M_8hr_D2")
# countstable <- merge(countstable, table3p, by.x = "cell", by.y = "tabletemp_orig", all = TRUE)
# colnames(countstable)[4] <- "CD4-T-cell_M_INFb_8hr_D1"
# countstable <- merge(countstable, table4p, by.x = "cell", by.y = "tabletemp_orig", all = TRUE)
# colnames(countstable)[5] <- "CD4-T-cell_M_INFb_8hr_D2"
# countstable <- merge(countstable, table5p, by.x = "cell", by.y = "tabletemp_orig", all = TRUE)
# colnames(countstable)[6] <- "CD4-T-cell_HIV_8hr_D1"
# countstable <- merge(countstable, table6p, by.x = "cell", by.y = "tabletemp_orig", all = TRUE)
# colnames(countstable)[7] <- "CD4-T-cell_HIV_8hr_D2"
# countstable <- merge(countstable, table7p, by.x = "cell", by.y = "tabletemp_orig", all = TRUE)
# colnames(countstable)[8] <- "CD4-T-cell_HIV_INFb_8hr_D1"
# countstable <- merge(countstable, table8p, by.x = "cell", by.y = "tabletemp_orig", all = TRUE)
# colnames(countstable)[9] <- "CD4-T-cell_HIV_INFb_8hr_D2"
# countstable[is.na(countstable)] <- 0
# counts_melt <- melt(countstable)

# ggplot(data = counts_melt, aes(x = cell, y = value, fill = variable)) +
#     geom_bar(stat = "identity", position = position_dodge(width = 1)) +
#     theme_Publication() +
#     theme(legend.text = element_text(size = 8)) +
#     scale_color_manual(values = c("firebrick1", "firebrick2", "firebrick3", "firebrick4", "steelblue1", "steelblue2", "steelblue3", "steelblue4")) +
#     scale_fill_manual(values = c("firebrick1", "firebrick2", "firebrick3", "firebrick4", "steelblue1", "steelblue2", "steelblue3", "steelblue4")) +
#     xlab("") +
#     ylim(c(0, 1)) +
#     ylab("% of Cells") +
#     theme(
#         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
#         axis.text.y = element_text(size = 12)
#     )
# ggsave(file.path(fig_results, "SupFigA.png"), height = 4, width = 6, dpi = 500)


## ------get zika percentages for each cell type
HIV_D1 <- get_virus_percentages(table5, "CD4-T-cell_HIV_8hr_D1")
HIV_D2 <- get_virus_percentages(table6, "CD4-T-cell_HIV_8hr_D2")
HIV_INFb_D1 <- get_virus_percentages(table7, "CD4-T-cell_HIV_INFb_8hr_D1")
HIV_INFb_D2 <- get_virus_percentages(table8, "CD4-T-cell_HIV_INFb_8hr_D2")
total_hiv <- rbind(HIV_D1, HIV_D2, HIV_INFb_D1, HIV_INFb_D2)

HIV_D1 <- get_virus_percentages_total(table5, "CD4-T-cell_HIV_8hr_D1")
HIV_D2 <- get_virus_percentages_total(table6, "CD4-T-cell_HIV_8hr_D2")
HIV_INFb_D1 <- get_virus_percentages_total(table7, "CD4-T-cell_HIV_INFb_8hr_D1")
HIV_INFb_D2 <- get_virus_percentages_total(table8, "CD4-T-cell_HIV_INFb_8hr_D2")
total_hiv_100 <- rbind(HIV_D1, HIV_D2, HIV_INFb_D1, HIV_INFb_D2)
samples <- str_replace_all(total_hiv_100$sample, "CD4-T-cell_HIV_8hr_D1", "HIV Infec. D1")
samples <- str_replace_all(samples, "CD4-T-cell_HIV_8hr_D2", "HIV Infec. D2")
samples <- str_replace_all(samples, "CD4-T-cell_HIV_INFb_8hr_D1", "IFN treat. HIV Infec. D1")
samples <- str_replace_all(samples, "CD4-T-cell_HIV_INFb_8hr_D2", "IFN treat. HIV Infec. D2")
total_hiv_100$sample <- samples

ggplot(data = total_hiv_100, aes(x = cell, y = Freq, fill = factor(z, levels = c("noHIV", "HIV")))) +
    geom_bar(stat = "identity") +
    theme_Publication() +
    facet_wrap(~sample, ncol = 4) +
    theme(strip.text = element_text(size = 10)) +
    theme(strip.background = element_rect(fill = c("white"))) +
    theme(legend.text = element_text(size = 8)) +
    theme(legend.key.size = unit(0.2, "cm")) +
    theme(legend.position = "top") +
    theme(legend.direction = "horizontal") +
    scale_color_manual(values = c("#ddd2d2", "red")) +
    scale_fill_manual(values = c("#ddd2d2", "red")) +
    xlab("") +
    labs(fill = "") +
    ylim(c(0, 1)) +
    ylab("% of Cells") +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11)
    )
ggsave(file.path(fig_results, "SupFig5c.png"), height = 4.5, width = 10, dpi = 500)
ggsave(file.path(fig_results, "SupFig5c.svg"), height = 4.5, width = 10, dpi = 500)
ggsave(file.path(fig_results, "SupFig5c.pdf"), height = 4.5, width = 10, dpi = 500)


## --Sup Fig 6
DefaultAssay(immune.combined_sd) <- "integrated"
Idents(immune.combined_sd) <- immune.combined_sd@meta.data$orig.ident
newidents <- str_replace_all(immune.combined_sd@meta.data$orig.ident, "CD4-T-cell_M_8hr_D1", "Mock D1")
newidents <- str_replace_all(newidents, "CD4-T-cell_M_8hr_D2", "Mock D2")
newidents <- str_replace_all(newidents, "CD4-T-cell_HIV_8hr_D1", "HIV Infec. D1")
newidents <- str_replace_all(newidents, "CD4-T-cell_HIV_8hr_D2", "HIV Infec. D2")
newidents <- str_replace_all(newidents, "CD4-T-cell_M_IFNb_8hr_D1", "IFN treat. D1")
newidents <- str_replace_all(newidents, "CD4-T-cell_M_IFNb_8hr_D2", "IFN treat. D2")
newidents <- str_replace_all(newidents, "CD4-T-cell_HIV_IFNb_8hr_D1", "IFN treat. HIV Infec. D1")
newidents <- str_replace_all(newidents, "CD4-T-cell_HIV_IFNb_8hr_D2", "IFN treat. HIV Infec. D2")

immune.combined_sd@meta.data$new.idents <- newidents
Idents(immune.combined_sd) <- immune.combined_sd@meta.data$new.idents
immune.combined_sd@active.ident <- factor(immune.combined_sd@active.ident, levels = c(
    "HIV Infec. D1", "HIV Infec. D2", "Mock D1", "Mock D2",
    "IFN treat. HIV Infec. D1", "IFN treat. HIV Infec. D2", "IFN treat. D1", "IFN treat. D2"
))
genes <- c(intersect(rownames(DEIFNb_mock), rownames(DEHIVIFNb_HIV)), setdiff(rownames(DEIFNb_mock), rownames(DEHIVIFNb_HIV)), 
    setdiff(rownames(DEHIVIFNb_HIV), rownames(DEIFNb_mock)))
DoHeatmap(immune.combined_sd,
    features = genes,
    slot = "scale.data", label = TRUE, draw.lines = TRUE, size = 4
) + theme(axis.text.y = element_text(size = 4.7)) + NoLegend()
ggsave(file.path(fig_results, "SupFig6.pdf"), height = 10, width = 8, units = "in", dpi = 300)
