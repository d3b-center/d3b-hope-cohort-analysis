# Author: Krutika Gaonkar
#
# Read in consensus snv calls to gather alterations in TP53 and NF1
# to evaluate classifier
# @params snvConsensus multi-caller consensus snv calls
# @params cnvConsensus multi-caller consensus cnv calls
# @params histologyFile histology file: histologies.tsv
# @params outputFolder output folder for alteration file
# @params gencode cds bed file from gencode

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("GenomicRanges"))

#### Source functions ----------------------------------------------------------
# We can use functions from the `snv-callers` module of the OpenPedCan project
# TODO: if a common util folder is established, use that instead
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
scratch_dir <- file.path(root_dir, "scratch")

source(file.path(root_dir, 
                 "utils",
                 "tmb_functions.R"))

#### Parse command line options ------------------------------------------------

option_list <- list(
  make_option(c("-s", "--snvConsensus"),type="character",
              help="Consensus snv calls (.tsv) "),
  make_option(c("-t", "--snvTumorOnly"),type="character",
              help="Tumor only snv calls (.tsv) "),
  make_option(c("-n","--cnvTumorOnly"),type="character",
              help="consensus cnv calls (.rds) "),
  make_option(c("-c","--cnv"),type="character",
              help="consensus cnv calls (rds) "),
  make_option(c("-e","--expr"),type="character",
              help="rna expression tpm or fpkm (.rds) "),
  make_option(c("-h","--histologyFile"),type="character",
              help="histology file for all samples (.tsv)"),
  make_option(c("-o","--outputFolder"),type="character",
              help="output folder for results "),
  make_option(c("-g","--gencode"),type="character",
              help="cds gencode bed file"),
  make_option(c("-r","--cohort"),type="character",
              help="list of cohorts to subset the files to (.tsv)")
)

opt <- parse_args(OptionParser(option_list=option_list,add_help_option = FALSE))
snvConsensusFile <- opt$snvConsensus
snvTumorOnlyFile <- opt$snvTumorOnly
expFile <- opt$expr
histologyFile <- opt$histologyFile
outputFolder <- opt$outputFolder
gencodeBed <- opt$gencode
cnvTumorOnlyFile <- opt$cnvTumorOnly
cnvFile <- opt$cnv

if(!dir.exists(outputFolder)){
  dir.create(outputFolder)
}

#### Generate files with TP53, NF1 mutations -----------------------------------

# read in consensus SNV files
keep_columns <- c("Chromosome",
                  "Start_Position",
                  "End_Position",
                  "Strand",
                  "Variant_Classification",
                  "Tumor_Sample_Barcode",
                  "Hugo_Symbol")

tumoronly_snv <- data.table::fread(snvTumorOnlyFile, select = keep_columns, tmpdir = scratch_dir) %>%
  dplyr::rename("Kids_First_Biospecimen_ID" = "Tumor_Sample_Barcode") 

consensus_snv <- data.table::fread(snvConsensusFile, select = keep_columns, tmpdir = scratch_dir) %>%
  dplyr::rename("Kids_First_Biospecimen_ID" = "Tumor_Sample_Barcode") %>%
  bind_rows(tumoronly_snv)

# read in CNV tumor only file
cnv_tumor_only <- readr::read_rds(cnvTumorOnlyFile) %>%
  as.data.frame() %>%
  dplyr::filter(!grepl('chrX|chrY', chr)) %>%
  dplyr::select(chr, start, end, 
                gene_symbol,
                Kids_First_Biospecimen_ID,
                status) 

# read in CNV file
cnv <- readr::read_rds(cnvFile) %>%
  as.data.frame() %>%
  dplyr::filter(!grepl('chrX|chrY', chr)) %>%
  dplyr::select(chr, start, end, 
                gene_symbol,
                Kids_First_Biospecimen_ID,
                status)

# combine two cnv files together
cnvConsensus <- cnv %>% 
  bind_rows(cnv_tumor_only) %>% 
  distinct() %>% 
  dplyr::rename(Chromosome = chr) %>% 
  dplyr::rename(Start_Position = start) %>% 
  dplyr::rename(End_Position = end) %>% 
  mutate(Chromosome = paste0("chr", Chromosome))

# gencode cds region BED file
gencode_cds <- read_tsv(gencodeBed, col_names = FALSE)

# histology file
histology <- read_tsv(histologyFile, guess_max = 100000)

# filter the MAF data.frame to only include entries that fall within the
# CDS bed file regions
coding_consensus_snv <- snv_ranges_filter(maf_df = consensus_snv,
                                          keep_ranges = gencode_cds)

# subset to TP53, removing silent mutations and mutations in introns
tp53_coding <- coding_consensus_snv %>%
  filter(Hugo_Symbol == "TP53") %>%
  filter(!(Variant_Classification %in% c("Silent", "Intron"))) %>% 
  dplyr::rename("Tumor_Sample_Barcode" = "Kids_First_Biospecimen_ID")
  

# subset to TP53 cnv loss and format to tp53_coding file format
tp53_loss<-cnvConsensus %>% 
  filter(gene_symbol == "TP53",
         status == "loss") %>%
  dplyr::rename("Tumor_Sample_Barcode" = "Kids_First_Biospecimen_ID",
                "Variant_Classification" = "status",
                "Hugo_Symbol" =  "gene_symbol")

# subset to NF1, removing silent mutations, mutations in introns, and missense
# mutations -- we exclude missense mutations because they are not annotated
# with OncoKB
# https://github.com/AlexsLemonade/OpenPBTA-analysis/pull/381#issuecomment-570748578
nf1_coding <- coding_consensus_snv %>%
  filter(Hugo_Symbol == "NF1") %>%
  filter(!(Variant_Classification %in% c("Silent",
                                         "Intron",
                                         "Missense_Mutation"))) %>% 
  dplyr::rename("Tumor_Sample_Barcode" = "Kids_First_Biospecimen_ID")

# subset to NF1 loss and format to nf1_coding file format
nf1_loss<-cnvConsensus %>% 
  filter(gene_symbol == "NF1",
         status == "loss") %>%
  dplyr::rename("Tumor_Sample_Barcode" = "Kids_First_Biospecimen_ID",
                "Variant_Classification" = "status",
                "Hugo_Symbol" =  "gene_symbol")

# include only the relevant columns from the MAF file and merge cnv loss dataframes as well
tp53_nf1_coding <- tp53_coding %>%
  bind_rows(tp53_loss,nf1_coding,nf1_loss)

# biospecimen IDs for tumor DNA-seq
bs_ids <- histology %>%
  filter(sample_type == "Tumor",
         experimental_strategy %in% c("WGS", "WXS")) %>%
  pull(Kids_First_Biospecimen_ID)


# all BS ids that are not in the data frame that contain the TP53 and NF1
# coding mutations should be labeled as not having either
bs_ids_without_mut <- setdiff(bs_ids,
                              unique(tp53_nf1_coding$Tumor_Sample_Barcode))

# add the TP53 and NF1 wildtype samples into the data.frame
tp53_nf1_coding <- bind_rows(tp53_nf1_coding,
                             data.frame(
                               Tumor_Sample_Barcode = bs_ids_without_mut,
                               Hugo_Symbol = "No_TP53_NF1_alt")
                             )

# save TP53 and NF1 SNV alterations
write_tsv(tp53_nf1_coding,
          file.path(outputFolder, "TP53_NF1_snv_alteration.tsv"))

# read in expression RDS file
rna <- readRDS(expFile)

# subset hist for those in rna matrix
hist_rna <- histology %>%
  filter(Kids_First_Biospecimen_ID %in% names(rna)) %>%
  ## add RNA-library as stranded for CPTAC samples
  mutate(RNA_library = case_when(grepl("^C3", Kids_First_Biospecimen_ID) ~ "stranded", 
                                 RNA_library == "poly-A stranded" ~ "poly-A-stranded",
                                 TRUE ~ as.character(RNA_library)))

# prepare rna-seq files
library_list <- hist_rna %>%
  pull(RNA_library) %>%
  unique()

# subset RNA file by library type
for (each in library_list){
  
  bs_ids <- hist_rna %>%
    filter(RNA_library == each) %>%
    pull(Kids_First_Biospecimen_ID)
  
  exp_subset <- rna[,bs_ids]
  
  write_rds(exp_subset, file.path(scratch_dir, 
                                  paste0("gene-expression-rsem-tpm-collapsed-", each, ".rds"))
  )
}
