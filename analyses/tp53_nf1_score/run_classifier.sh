#!/bin/bash

# Adapted from OpenPedCan tp53_nf1_score module with some adjustment

# K S Gaonkar

# This analysis applies the TP53 inactivation classifier  https://linkinghub.elsevier.com/retrieve/pii/S2211124718304376
# NF1 inactivation classifier https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3519-7
# Predicts TP53 and NF1 inactivation score per polya and stranded RNAseq samples

# The script takes one environment variable, `OPENPBTA_BASE_SUBTYPING`, if value is 1 then
# uses histologies-base.tsv for subtyping if value is 0 runs all modules with histologies.tsv(Default)

set -e
set -o pipefail

# The script can tak three environment variables:
# `OPENPBTA_BASE_SUBTYPING`: 
#     if value is 1, then uses `pbta-histologies-base.tsv` for subtyping. 
#     if value is 0 (DEFAULT), runs module with `pbta-histologies.tsv`
# `OPENPEDCAN_POLYA_STRAND`: 
#     if value is 1 (DEFAULT), runs the POLYA_STRAND steps
#     if value is 0, skips the POLYA_STRAND steps

POLYA_STRAND=${OPENPEDCAN_POLYA_STRAND:-1}
RUN_FOR_SUBTYPING=${OPENPBTA_BASE_SUBTYPING:-0}


# Temporary solution to for the python3 rpy2 package to work
# rpy2 not loading R is a shared library (libR.so) when imported
# currently no viable solution for versions for rpy2 in the Opendpedcan docker debian release
export R_HOME='/usr/local/lib/R'
export LD_LIBRARY_PATH=${R_HOME}/lib:${LD_LIBRARY_PATH}


input_dir="input"
data_dir="../../data"
scratch_dir="../../scratch"
# cds gencode bed file  
cds_file="${input_dir}/gencode.v39.primary_assembly.annotation.bed"
snvconsensus_file="${data_dir}/Hope-snv-consensus-plus-hotspots.maf.tsv.gz"
snvtumoronly_file="${data_dir}/Hope-tumor-only-snv-mutect2.maf.tsv.gz"
cnvconsensus_file="${data_dir}/Hope-cnv-controlfreec-tumor-only.rds"
cnv_file="${data_dir}/Hope-cnv-controlfreec.rds"
collapsed_rna_file="${data_dir}/Hope-and-CPTAC-GBM-gene-expression-rsem-tpm-collapsed.rds"
histology_file="${data_dir}/Hope-GBM-histologies-base.tsv" 


# Convert GTF to BED file
# Here we are only extracting lines with as a CDS i.e. are coded in protein
# gunzip -c ${data_dir}/gencode.v39.primary_assembly.annotation.gtf.gz \
#  | awk '$3 ~ /CDS/' \
#  | convert2bed --do-not-sort --input=gtf - \
#  > $cds_file

# Prep the SNV consensus data for evaluation downstream
Rscript --vanilla 00-tp53-nf1-alterations.R \
  --snvConsensus ${snvconsensus_file} \
  --snvTumorOnly ${snvtumoronly_file} \
  --cnvTumorOnly ${cnvconsensus_file} \
  --cnv ${cnv_file} \
  --histologyFile ${histology_file} \
  --outputFolder results \
  --gencode ${cds_file} \
  --expr ${collapsed_rna_file}
  
  
echo "### finish 00-tp53-nf1-alterations.R ###"

# Define RNA library files, which result from the script above
collapsed_stranded="${scratch_dir}/gene-expression-rsem-tpm-collapsed-stranded.rds"
collapsed_polya_stranded="${scratch_dir}/gene-expression-rsem-tpm-collapsed-poly-A-stranded.rds"
collapsed_exome_capture="${scratch_dir}/gene-expression-rsem-tpm-collapsed-exome_capture.rds"

# Run classifier and ROC plotting for RNA data - currently, we have 3 types of RNA libraries. 
# No polyA in Hope cohort
# We should add to this if we get more types.
python3 01-apply-classifier.py -f ${collapsed_stranded}
echo "### finish 01-apply-classifier.py 01 ###"

python3 01-apply-classifier.py -f ${collapsed_exome_capture}
echo "### finish 01-apply-classifier.py 03 ###"

python3 01-apply-classifier.py -f ${collapsed_polya_stranded}
echo "### finish 01-apply-classifier.py 04 ###"

# check correlation expression and scores
Rscript -e "rmarkdown::render('02-qc-rna_expression_score.Rmd',params=list(base_run = $RUN_FOR_SUBTYPING, ci_run = $POLYA_STRAND))"
echo "02-qc"

# subset cnv where tp53 is lost
Rscript -e "rmarkdown::render('03-tp53-cnv-loss-domain.Rmd',params=list(base_run = $RUN_FOR_SUBTYPING))"
echo "03 tp53"

# subset SV where tp53 is lost
Rscript -e "rmarkdown::render('04-tp53-sv-loss.Rmd',params=list(base_run = $RUN_FOR_SUBTYPING))"
echo "04 tp53"

# gather TP53 altered status
Rscript -e "rmarkdown::render('05-tp53-altered-annotation.Rmd',params=list(base_run = $RUN_FOR_SUBTYPING))"
echo "05 tp53"

# evaluate classifer scores for stranded data
# Skip poly-A, poly-A stranded, exome capture steps in CI
if [ "$POLYA_STRAND" -gt "0" ]; then
  python3 06-evaluate-classifier.py -s results/tp53_altered_status.tsv -f results/gene-expression-rsem-tpm-collapsed-poly-A-stranded_classifier_scores.tsv -c ${histology_file} -o polya_stranded
  echo "6"
  python3 06-evaluate-classifier.py -s results/tp53_altered_status.tsv -f results/gene-expression-rsem-tpm-collapsed-exome_capture_classifier_scores.tsv -c ${histology_file} -o exome_capture
  echo "7"
  python3 06-evaluate-classifier.py -s results/tp53_altered_status.tsv -f results/gene-expression-rsem-tpm-collapsed-stranded_classifier_scores.tsv -c ${histology_file} -o stranded
  echo "8"
fi
