# Function: Script to generate Oncogrid plot

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(tidyverse)
  library(circlize)
})

# set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analyses_dir <- file.path(root_dir, "analyses", "oncoplots")
input_dir <- file.path(analyses_dir, "results")
output_dir <- file.path(analyses_dir, "results", "cascade_plots")
dir.create(output_dir, recursive = T, showWarnings = F)

# matrix
mat = read.table(file.path(input_dir, "oncoprint.txt"),  header = TRUE, stringsAsFactors=FALSE, sep = "\t",check.names = FALSE)
mat[is.na(mat)] = ""
rownames(mat) = mat[, 1]
mat = mat[, -1]
mat = t(as.matrix(mat))

# remove RNA related annotations
mat[mat %in% c("OVE", "UNE", "FUS", "GAI", "LOS")] <- ""
mat = gsub(";UNE|;OVE|;FUS", "", mat)
mat = gsub("GAI|LOS|;LOS|GAI;LOS|LOS;GAI", "", mat)
mat[mat != ""] <- "Mutation"

keep = apply(mat, 1, FUN = function(x) length(unique(x)) == 2)
mat = mat[keep,]

# get unique values
apply(mat, 1, unique) %>%  unlist %>% unique()

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w, h, gp = gpar(fill = "#ffffff",col= "#595959"))
  },
  Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.3, "mm"), gp = gpar(fill = "black", col = NA))
  }
)
col = c("Mutation" = "black")

# read annotation and TMB info
annot_info <- read.delim(file.path(input_dir, "annotation.txt"), header = TRUE, check.names = TRUE)
annot_info <- annot_info %>%
  filter(Sample %in% colnames(mat),
         Sequencing_Experiment != "RNA-Seq") %>%
  remove_rownames() %>%
  column_to_rownames('Sample') %>%
  as.data.frame()
samples_to_use <- intersect(rownames(annot_info), colnames(mat))
mat <- mat[,samples_to_use]
annot_info <- annot_info[samples_to_use,]

# only plot top 20 genes 
snv_genes_to_keep = apply(mat, 1, FUN = function(x) (length(grep("Mutation", x))/ncol(mat))*100)
snv_genes_to_keep <- names(sort(snv_genes_to_keep, decreasing = TRUE)[1:20])
mat <- mat[which(rownames(mat) %in% snv_genes_to_keep),]

# annotation 1
col_fun_tmb = colorRamp2(c(0, max(annot_info$TMB, na.rm = T)), c("white", "magenta3"))
annot_info$Age <- factor(annot_info$Age, levels = c("[0,15]", "(15,26]", "(26,40]"))
annot_info$Molecular_Subtype <- as.character(annot_info$Molecular_Subtype)
annot_info$Cancer_Group <- as.character(annot_info$Cancer_Group)

ha = HeatmapAnnotation(df = annot_info %>% dplyr::select(-c(Sequencing_Experiment)), col = list(
  TMB = col_fun_tmb,
  Diagnosis = c("High-grade glioma/astrocytoma (WHO grade III/IV)" = "lightseagreen",
                "Diffuse Midline Glioma (WHO grade III/IV)" = "darkgreen",
                "Astrocytoma;Oligoastrocytoma (WHO grade III)" = "mediumorchid2",
                "Astrocytoma (WHO grade III/IV)" = "#5fff57", 
                "Glioblastoma (WHO grade IV)" = "#f268d6",
                "Pleomorphic xanthoastrocytoma (WHO grade II/III)" = "#005082"),
  Molecular_Subtype = c("DMG, H3 K28" = "#053061",
                        "DHG, H3 G35, TP53" = "#A6761D",
                        "HGG, H3 wildtype" = "#4393c3",
                        "HGG, H3 wildtype, TP53" = "darkgreen",
                        "DMG, H3 K28, TP53" = "#BC80BD",
                        "HGG, IDH, TP53" = "#FFFF99",
                        "IHG, NTRK-altered, TP53"  = "#E7298A",
                        "IHG, NTRK-altered" = "#f4a582",
                        "IHG, ROS1-altered" = "#d6604d",
                        "IHG, ALK-altered" = "#E31A1C",
                        "PXA" = "#67001f",
                        "HGG, IDH" = "#B3DE69",
                        "NA" = "#f1f1f1"),
  Diagnosis_Type = c("Initial CNS Tumor" = "#cee397",
                     "Progressive" = "#827397",
                     "Recurrence" = "#363062",
                     "Second Malignancy" = "#005082"),
  Tumor_Location = c("Cortical" = "#D4806C",
                     "Other/Multiple locations/NOS" = "#7C8F97",
                     "Midline" = "#344C68",
                     "Cerebellar" = "#94004C"),
  CNS_region = c("Posterior fossa" = "#D4806C",
                 "Other" = "#7C8F97",
                 "Midline" = "#344C68",
                 "Hemispheric" = "#94004C",
                 "Mixed" = "darkgreen"),
  Cancer_Group = c("DMG" = "#053061",
                   "DHG" = "#A6761D",
                   "HGG" = "#4393c3",
                   "IHG" = "#E7298A",
                   "NA" = "#f1f1f1"),
  Sex = c("Male" = "#0707CF",
          "Female" = "#CC0303"),
  Age = c("[0,15]" = "#C7E9C0",
          "(15,26]" = "#74C476",
          "(26,40]" = "#238B45")),
  annotation_name_gp = gpar(fontsize = 10),
  gp = gpar(col = "#595959"), 
  simple_anno_size = unit(4, "mm"), 
  annotation_name_side = "left")

# oncoprint
ht = oncoPrint(mat, get_type = function(x)strsplit(x, ";")[[1]],
               alter_fun = alter_fun, 
               col = col, 
               show_column_names = TRUE, 
               column_names_gp = gpar(fontsize = 10),
               row_names_gp = gpar(fontsize = 10),
               column_names_side = "top",
               top_annotation = ha,
               right_annotation = NULL,
               row_names_side = "left",
               pct_side = "right",
               remove_empty_rows = TRUE,
               heatmap_legend_param = list(title = "Alteration", nrow = 9, title_position = "topleft", direction = "horizontal",
                                           at = c("Mutation"),
                                           labels = c("SNV")
               ))

pdf(file = file.path(output_dir, "cascade_plot.pdf"), width = 15, height = 8) 
draw(ht,merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# order by Age 
ht = oncoPrint(mat, get_type = function(x)strsplit(x, ";")[[1]],
               alter_fun = alter_fun, 
               column_split = annot_info$Age,
               col = col, 
               show_column_names = TRUE, 
               column_names_gp = gpar(fontsize = 10),
               row_names_gp = gpar(fontsize = 10),
               column_names_side = "top",
               top_annotation = ha,
               right_annotation = NULL,
               row_names_side = "left",
               pct_side = "right",
               remove_empty_rows = TRUE,
               heatmap_legend_param = list(title = "Alteration", nrow = 9, title_position = "topleft", direction = "horizontal",
                                           at = c("Mutation"),
                                           labels = c("SNV")
               ))
pdf(file = file.path(output_dir, "cascade_orderby_age.pdf"), width = 15, height = 8) 
draw(ht,merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# order by Sex
ht = oncoPrint(mat, get_type = function(x)strsplit(x, ";")[[1]],
               alter_fun = alter_fun, 
               column_split = annot_info$Sex,
               col = col, 
               show_column_names = TRUE, 
               column_names_gp = gpar(fontsize = 10),
               row_names_gp = gpar(fontsize = 10),
               column_names_side = "top",
               top_annotation = ha,
               right_annotation = NULL,
               row_names_side = "left",
               pct_side = "right",
               remove_empty_rows = TRUE,
               heatmap_legend_param = list(title = "Alteration", nrow = 9, title_position = "topleft", direction = "horizontal",
                                           at = c("Mutation"),
                                           labels = c("SNV")
               ))
pdf(file = file.path(output_dir, "cascade_orderby_sex.pdf"), width = 15, height = 8) 
draw(ht,merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
