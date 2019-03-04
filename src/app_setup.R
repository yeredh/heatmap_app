# 4/23/2017
#
# save objects for probe and gene counts heatmaps


rm(list=ls())
options(stringsAsFactors = F)
ody_dir = "/Volumes/irizarryfs01/ypitajuarez/"
# ody_dir = "/net/irizarryfs01/srv/export/irizarryfs01_backed_up/share_root/ypitajuarez/"


# ==== Orthologs ====
mah_genes = readRDS("Irizarry/Project II/data/lin_ortholog_features.RDS")

library(preprocessCore)
K = 1
# ==== Lin2 Gene Counts ====
## SAMPLE ANNOTATION
lin2_pd = readRDS(paste0(ody_dir,"Encode/data/enc_src.RDS"))
# remove problematic samples
lin2_pd = lin2_pd[!(lin2_pd$Biosample.term.name %in% c("pancreas","spleen")),]
lin2_pd = lin2_pd[lin2_pd$Technical.replicate == 2,]
lin2_pd = lin2_pd[!duplicated(lin2_pd$Experiment.accession),]

## HUMAN GENE COUNTS
lin2_hs_gene_count = readRDS(paste0(ody_dir,"Encode/ENCSR/Lin2/hs_gene_countF02"))
rownames(lin2_hs_gene_count$counts) = gsub("([A-Z0-9]+)\\..*","\\1",rownames(lin2_hs_gene_count$counts))

## MOUSE GENE COUNTS
lin2_mm_gene_count = readRDS(paste0(ody_dir,"Encode/ENCSR/Lin2/mm_gene_countF02"))

## PAIR BY ORTHOLOGS
hs_ind02 = match(mah_genes$hs,rownames(lin2_hs_gene_count$counts))
mm_ind02 = match(mah_genes$mm,rownames(lin2_mm_gene_count$counts))
lin2_gene_count = cbind(
  lin2_hs_gene_count$counts[hs_ind02,],
  lin2_mm_gene_count$counts[mm_ind02,]
)
rownames(lin2_gene_count) = paste(
  rownames(lin2_hs_gene_count$counts)[hs_ind02],
  rownames(lin2_mm_gene_count$counts)[mm_ind02],
  sep=":"
)
identical(lin2_pd$Experiment.accession,colnames(lin2_gene_count))

## NORMALIZATION
# log-transformation
lin2_gene_log2trans = log2(lin2_gene_count+K)
# Lin et al.: transform, then quantile normalize
lin2_gene_QNlog2trans = normalize.quantiles(log2(lin2_gene_count+K))
dimnames(lin2_gene_QNlog2trans) = dimnames(lin2_gene_count)

# alignment setting
al=14

# ==== Lin2: Probe Counts ====
## HUMAN PROBE COUNTS
lin2_hs_probe_count = readRDS(paste0(ody_dir,"Encode/output/al_tests/union_exons/counts/lin2_hs_fcLst_union_exons_AE3_",al,".RDS"))

## MOUSE PROBE COUNTS
lin2_mm_probe_count = readRDS(paste0(ody_dir,"Encode/output/al_tests/union_exons/counts/lin2_mm_fcLst_union_exons_AE3_",al,".RDS"))

## PAIR BY ORTHOLOGS
if(identical(rownames(lin2_hs_probe_count$counts),rownames(lin2_mm_probe_count$counts))){
  lin2_probe_counts = cbind(
    lin2_hs_probe_count$counts,
    lin2_mm_probe_count$counts
  )
}else{
  orth_ind = match(rownames(lin2_hs_probe_count$counts),rownames(lin2_mm_probe_count$counts))
  lin2_probe_counts = cbind(
    lin2_hs_probe_count$counts,
    lin2_mm_probe_count$counts[orth_ind,]
  )
}
if(!identical(colnames(lin2_probe_counts),lin2_pd$Experiment.accession)){
  sind = match(lin2_pd$Experiment.accession,colnames(lin2_probe_counts))
  lin2_probe_counts = lin2_probe_counts[,sind]
}

## NORMALIZATION
# log-transformation
lin2_probe_log2trans = log2(lin2_probe_counts+K)
# Lin et al.: transform, then quantile normalize
lin2_probe_QNlog2trans = normalize.quantiles(log2(lin2_probe_counts+K))
dimnames(lin2_probe_QNlog2trans) = dimnames(lin2_probe_counts)


# ==== Save Data ====
# format sample annotation
lin2_pd = lin2_pd[,c("Experiment.accession","Biosample.organism","Biosample.term.name")]
colnames(lin2_pd) = c("experiment","species","tissue")

# orthologs
saveRDS(mah_genes,"Shiny/paper_figures/data/mah_genes.RDS")
# sample annotation
saveRDS(lin2_pd,"Shiny/paper_figures/data/lin2_pd.RDS")
# count matrices
saveRDS(lin2_probe_QNlog2trans,"Shiny/paper_figures/data/lin2_probe_QNlog2trans.RDS")
saveRDS(lin2_gene_QNlog2trans,"Shiny/paper_figures/data/lin2_gene_QNlog2trans.RDS")
