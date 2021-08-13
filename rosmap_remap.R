#-=-=-=-=-=-=-=-=-=-=-=-=-=--==--=-=-=-=-=-=-=-=---=-=--=-=-=-=--=-=-=-=--=-=-=-
# Path and libraries
#-=-=-=-=-=-=-=-=-=-=-=-=-=--==--=-=-=-=-=-=-=-=---=-=--=-=-=-=--=-=-=-=--=-=-=-
devtools::load_all("~/Documents/packages//MethReg/")
library(dplyr)
library(MethReg)
library(GenomicRanges)
library(SummarizedExperiment)
library(microbenchmark)

# path.dropbox <- dir("~",pattern = "Dropbox",full.names = TRUE)
path.dropbox <- "~/TBL Dropbox/Tiago Silva/"
path.project <- file.path(path.dropbox,"/PanCancer/MethReg-useCase/")
path.data <- file.path(path.project,"data")
path.usecase <- file.path(path.project,"/UseCase2_rosmap_remap//")
path.plots <- file.path(path.usecase,"/plots/")
path.tables <- file.path(path.usecase,"/results/")
path.tables.promoter <- file.path(path.usecase,"/results/promoter")
path.ewas.analysis <-  file.path(path.dropbox,"/coMethDMR_metaAnalysis/")
setwd(path.project)

#-=-=-=-=-=-=-=-=-=-=-=-=-=--==--=-=-=-=-=-=-=-=---=-=--=-=-=-=--=-=-=-=--=-=-=-
# Data from other projects
#-=-=-=-=-=-=-=-=-=-=-=-=-=--==--=-=-=-=-=-=-=-=---=-=--=-=-=-=--=-=-=-=--=-=-=-
london.blood <- readr::read_csv("~/Dropbox (BBSR)/Tiago Silva/coMethDMR_metaAnalysis/code_validation/Meta_analysis_code/London_blood_brain_correlation_results/London_blood_brain_correlation_cpgs.csv")
colnames(london.blood) <- paste0("London_blood_",colnames(london.blood))
colnames(london.blood)[1] <- "probeID"

#-=-=-=-=-=-=-=-=-=-=-=-=-=--==--=-=-=-=-=-=-=-=---=-=--=-=-=-=--=-=-=-=--=-=-=-
# Aux functions
#-=-=-=-=-=-=-=-=-=-=-=-=-=--==--=-=-=-=-=-=-=-=---=-=--=-=-=-=--=-=-=-=--=-=-=-

add_annot_cpgs <- function(results){
  results$UCSC_RefGene_Group <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Other[results$probeID,]$UCSC_RefGene_Group
  results$UCSC_RefGene_Name <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Other[results$probeID,]$UCSC_RefGene_Name
  results$Relation_to_Island <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Islands.UCSC[results$probeID,]$Relation_to_Island
  results
}


update_met_IQR <- function(results, dnam){
  iqr <- MethReg:::calculate_IQR(dnam)
  results$res.met.IQR <- results$met.IQR
  results$met.IQR <- iqr$IQR[match(results$regionID, iqr$ID)]
  results <- results %>% dplyr::relocate(dplyr::contains("IQR"), .after = last_col())
  return(results)
}


add_percent_zero_q1_q4 <- function(results, dnam, exp){
  
  aux <- plyr::adply(
    unique(results[,c("probeID","target")]),
    .margins = 1,
    .fun = function(row) {
      rna.target <- exp[rownames(exp) == row$target, , drop = FALSE]
      met <- dnam[rownames(dnam) == as.character(row$probeID), ]
      data <- data.frame(
        rna.target = rna.target %>% as.numeric,
        met = met %>% as.numeric
      )
      quant.met <-  quantile(data$met,na.rm = TRUE)
      low.cutoff <- quant.met[2]
      upper.cutoff <- quant.met[4]
      
      data.high.low <- data %>% filter(.data$met <= low.cutoff | .data$met >= upper.cutoff)
      data.high.low$metGrp <- ifelse(data.high.low$met <= low.cutoff, 0, 1)
      pct.zeros.in.quant.samples <- sum(
        data.high.low$rna.target == 0,
        na.rm = TRUE) / nrow(data.high.low)
      
      data.frame("% of 0 target genes (Q1 and Q4)" = paste0(round(pct.zeros.in.quant.samples * 100,digits = 2)," %"))
    }
  )
  results$`% of 0 residual target genes (Q1 and Q4)` <- results$`% of 0 target genes (Q1 and Q4)`
  results$`% of 0 target genes (Q1 and Q4)` <- aux$X..of.0.target.genes..Q1.and.Q4.[
    match(paste0(results$probeID,results$target),paste0(aux$probeID,aux$target))
    ]
  return(results)
}

#-------------------------------------------------------------------------------
# Read EWAS results
#-------------------------------------------------------------------------------
meta.analysis.cpgs <- readxl::read_xlsx(
  file.path(
    path.ewas.analysis,
    "/DRAFT_REVISION_NatComm_9-15-2020/Supplementary Tables - ALL_8-7-2020-no-highlights.xlsx"
  ), 
  skip = 3
)

#-------------------------------------------------------------------------------
# Read Rosmap data
#-------------------------------------------------------------------------------
load(file = file.path(path.usecase,"data/matched_normalized_and_residuals_data.rda"))

iqr <- MethReg:::calculate_IQR(matched.dnam)
cpgs.iqr.higher.0.03 <- iqr$ID[iqr$IQR > 0.03]

meta.analysis.cpgs <- meta.analysis.cpgs[meta.analysis.cpgs$cpg %in% cpgs.iqr.higher.0.03,]

resid.met.cpg.awas <- residuals.matched.met[rownames(residuals.matched.met) %in% meta.analysis.cpgs$cpg,]
all(colnames(residuals.matched.exp) == colnames(resid.met.cpg.awas))

hm450.hg19 <- MethReg:::get_met_probes_info(genome = "hg19", arrayType = "450k")
rownames(resid.met.cpg.awas) <- MethReg:::make_names_from_granges(hm450.hg19[rownames(resid.met.cpg.awas) ,])

rnaseq.tf.es.gsva.ensg <- rnaseq.tf.es.gsva
rownames(rnaseq.tf.es.gsva.ensg) <- MethReg:::map_symbol_to_ensg(rownames(rnaseq.tf.es.gsva))
rnaseq.tf.es.gsva.ensg <- rnaseq.tf.es.gsva.ensg[!is.na(rownames(rnaseq.tf.es.gsva.ensg)),]


matched.dnam.with.regions <- matched.dnam
hm450.hg19 <- MethReg:::get_met_probes_info(genome = "hg19", arrayType = "450k")
rownames(matched.dnam.with.regions) <- MethReg:::make_names_from_granges(hm450.hg19[rownames(matched.dnam.with.regions) ,])

#-------------------------------------------------------------------------------
# Get triplets using remap
#-------------------------------------------------------------------------------
library(ReMapEnrich)
remapCatalog2018hg19 <- downloadRemapCatalog(path.data, assembly = "hg19")
remapCatalog <- bedToGranges(remapCatalog2018hg19)

#-------------------------------------------------------------------------------
# Promoter analysis, triplets using remap
#-------------------------------------------------------------------------------
triplet.promoter.ewas <- create_triplet_distance_based(
  region = hm450.hg19[meta.analysis.cpgs$cpg,],
  genome = "hg19",
  target.method =  "genes.promoter.overlap",
  TF.peaks.gr = remapCatalog,
  motif.search.window.size = 500
) 
dim(triplet.promoter.ewas) # 69928 triplets
triplet.promoter.ewas <- triplet.promoter.ewas %>% dplyr::filter(.data$TF %in% rownames(rnaseq.tf.es.gsva.ensg))
dim(triplet.promoter.ewas) # 57715 triplets
triplet.promoter.ewas$probeID <- names(hm450.hg19)[match(triplet.promoter.ewas$regionID,make_names_from_granges(hm450.hg19))]


triplet.promoter.ewas$TF %>% unique %>% length # 286
triplet.promoter.ewas$regionID %>% unique %>% length # 1404
triplet.promoter.ewas$target %>% unique %>% length # 1164


file.promoter <- file.path(path.tables.promoter, "ROSMAP_and_remap_promoter_analysis_using_TF_es_gsva_all_triplet.csv")

if (!file.exists(file.promoter)) {
  cores <- 10
  results.promoter.analysis <- 
    triplet.promoter.ewas %>% 
    cor_tf_target_gene(
      exp = residuals.matched.exp,
      tf.activity.es = rnaseq.tf.es.gsva.ensg,
      cores = cores
    ) %>% cor_dnam_target_gene(
      dnam = resid.met.cpg.awas,
      exp = residuals.matched.exp,
      cores = cores,
      filter.results = FALSE, 
      min.cor.estimate = 0.2,
      min.cor.pval = 0.05
    ) %>%  interaction_model(
      dnam = resid.met.cpg.awas,
      exp = residuals.matched.exp,
      tf.activity.es = rnaseq.tf.es.gsva.ensg,
      cores = cores,
      filter.correlated.tf.exp.dnam = FALSE,
      filter.triplet.by.sig.term = FALSE,
      sig.threshold = 0.05,
      fdr = TRUE,
      stage.wise.analysis = TRUE
    ) %>%  stratified_model(
      dnam = resid.met.cpg.awas,
      exp = residuals.matched.exp,
      tf.activity.es = rnaseq.tf.es.gsva.ensg,
      cores = cores
    )
  
  results.promoter.analysis <- results.promoter.analysis %>% add_annot_cpgs() %>%  
    add_percent_zero_q1_q4(dnam = matched.dnam, exp = matched.exp.log2) %>%
    update_met_IQR(dnam = matched.dnam.with.regions)
  
  results.promoter.analysis$RLM_TF_fdr <- p.adjust(results.promoter.analysis$RLM_TF_pvalue,method = "fdr")
  results.promoter.analysis$RLM_DNAmGroup_fdr <- p.adjust(results.promoter.analysis$RLM_DNAmGroup_pvalue,method = "fdr")
  results.promoter.analysis$`RLM_DNAmGroup:TF_fdr` <- p.adjust(results.promoter.analysis$`RLM_DNAmGroup:TF_pvalue`,method = "fdr")
  
  readr:::write_csv(
    x = results.promoter.analysis,
    file = file.promoter
  )
  results.promoter.analysis.sig.fdr.int <- results.promoter.analysis %>%
    dplyr::filter(`RLM_DNAmGroup:TF_fdr` < 0.05)
  
  results.promoter.analysis.sig.fdr.int.with.blood <- merge(results.promoter.analysis.sig.fdr.int, london.blood)
  readr:::write_csv(
    x = results.promoter.analysis.sig.fdr.int.with.blood,
    file = gsub("all_triplet","sig_fdr_int_triplet_with_london_blood",file.promoter)
  )
  readr:::write_csv(
    x = results.promoter.analysis.sig.fdr.int,
    file = gsub("all_triplet","sig_fdr_int_triplet",file.promoter)
  )
  
  results.promoter.analysis.sig.int <- results.promoter.analysis %>%
    dplyr::filter(`RLM_DNAmGroup:TF_triplet_stage_wise_adj_pvalue` < 0.05)
  
  readr:::write_csv(
    x = results.promoter.analysis.sig.int,
    file = gsub("all_triplet","sig_stage_wise_fdr_int_triplet",file.promoter)
  )
  
  results.promoter.analysis.sig <- results.promoter.analysis %>%
    filter_at(vars(contains("triplet_stage")), any_vars(. < 0.05))
  
  readr:::write_csv(
    x = results.promoter.analysis.sig,
    file = gsub("all_triplet", "sig_any_stage_wise_triplet", file.promoter)
  )
} else {
  results.promoter.analysis <- readr::read_csv(
    file = file.promoter,
    col_types = c(`RLM_DNAmGroup:TF_triplet_stage_wise_adj_pvalue` = "n")
  )
  
  results.promoter.analysis$is.enhancer.EhnancerAtlas.Brain.neuro <- results.promoter.analysis$probeID %in% enhancer.probes
}


results.promoter.analysis.sig.int <- results.promoter.analysis.sig.int[order(results.promoter.analysis.sig.int$`quant_triplet_stage_wise_adj_pval_metGrp:es.tf`),]



#-------------------------------------------------------------------------------
# Distal analysis, triplets using remap
#-------------------------------------------------------------------------------
triplet.distal.ewas <- create_triplet_distance_based(
  region = hm450.hg19[meta.analysis.cpgs$cpg,],
  genome = "hg19",
  target.method =  "nearby.genes",
  TF.peaks.gr = remapCatalog,
  motif.search.window.size = 500,
  target.rm.promoter.regions.from.distal.linking = TRUE
) 
dim(triplet.distal.ewas) # 1,422,184 triplets
triplet.distal.ewas <- triplet.distal.ewas %>% dplyr::filter(.data$TF %in% rownames(rnaseq.tf.es.gsva.ensg))
dim(triplet.distal.ewas) # 946,848  triplets
triplet.distal.ewas$probeID <- names(hm450.hg19)[match(triplet.distal.ewas$regionID,make_names_from_granges(hm450.hg19))]


triplet.distal.ewas$TF %>% unique %>% length # 288
triplet.distal.ewas$regionID %>% unique %>% length # 1756
triplet.distal.ewas$target %>% unique %>% length # 9,934


file.distal <- file.path(path.tables, "ROSMAP_and_remap_distal_analysis_using_TF_es_gsva_all_triplet.csv")

if (!file.exists(file.distal)) {
  cores <- 4
  results.distal.analysis <- 
    triplet.distal.ewas %>% 
    cor_tf_target_gene(
      exp = residuals.matched.exp,
      tf.activity.es = rnaseq.tf.es.gsva.ensg,
      cores = cores
    ) %>% cor_dnam_target_gene(
      dnam = resid.met.cpg.awas,
      exp = residuals.matched.exp,
      cores = cores,
      filter.results = FALSE, 
      min.cor.estimate = 0.2,
      min.cor.pval = 0.05
    ) %>%  interaction_model(
      dnam = resid.met.cpg.awas,
      exp = residuals.matched.exp,
      tf.activity.es = rnaseq.tf.es.gsva.ensg,
      cores = cores,
      filter.correlated.tf.exp.dnam = FALSE,
      filter.triplet.by.sig.term = FALSE,
      sig.threshold = 0.05,
      fdr = TRUE,
      stage.wise.analysis = TRUE
    ) %>%  stratified_model(
      dnam = resid.met.cpg.awas,
      exp = residuals.matched.exp,
      tf.activity.es = rnaseq.tf.es.gsva.ensg,
      cores = cores
    )
  
  #results.distal.analysis <- results.distal.analysis[results.distal.analysis$probeID %in% cpgs.iqr.higher.0.03,]
  results.distal.analysis <- results.distal.analysis %>% 
    add_annot_cpgs() %>%  
    add_percent_zero_q1_q4(dnam = matched.dnam, exp = matched.exp.log2)  %>%
    update_met_IQR(dnam = matched.dnam.with.regions)
  
  
  results.distal.analysis$RLM_TF_fdr <- p.adjust(results.distal.analysis$RLM_TF_pvalue,method = "fdr")
  results.distal.analysis$RLM_DNAmGroup_fdr <- p.adjust(results.distal.analysis$RLM_DNAmGroup_pvalue,method = "fdr")
  results.distal.analysis$`RLM_DNAmGroup:TF_fdr` <- p.adjust(results.distal.analysis$`RLM_DNAmGroup:TF_pvalue`,method = "fdr")
  
  readr:::write_csv(
    x = results.distal.analysis,
    file = file.distal
  )
  
  results.distal.analysis.sig.fdr.int <- results.distal.analysis %>%
    dplyr::filter(`RLM_DNAmGroup:TF_fdr` < 0.05)

  results.distal.analysis.sig.fdr.int.with.blood <- merge(results.distal.analysis.sig.fdr.int, london.blood)
  readr:::write_csv(
    x = results.distal.analysis.sig.fdr.int.with.blood,
    file = gsub("all_triplet","sig_fdr_int_triplet_with_london_blood",file.distal)
  )
  
  
  readr:::write_csv(
    x = results.distal.analysis.sig.fdr.int,
    file = gsub("all_triplet","sig_fdr_int_triplet",file.distal)
  )
  
  results.distal.analysis.sig.int <- results.distal.analysis %>%
    dplyr::filter(`RLM_DNAmGroup:TF_triplet_stage_wise_adj_pvalue` < 0.05)
  
  readr:::write_csv(
    x = results.distal.analysis.sig.int,
    file = gsub("all_triplet","sig_stage_wise_fdr_int_triplet",file.distal)
  )
  
  results.distal.analysis.sig <- results.distal.analysis %>%
    filter_at(vars(contains("triplet_stage")), any_vars(. < 0.05))
  
  readr:::write_csv(
    x = results.distal.analysis.sig,
    file = gsub("all_triplet", "sig_any_stage_wise_triplet", file.distal)
  )
}  else {
  results.distal.analysis <- readr::read_csv(
    file = file.distal,
    col_types = c(`RLM_DNAmGroup:TF_triplet_stage_wise_adj_pvalue` = "n")
  )
  
  results.distal.analysis$is.enhancer.EhnancerAtlas.Brain.neuro <- results.distal.analysis$probeID %in% enhancer.probes
}


results.distal.analysis.sig.int <- results.distal.analysis.sig.int[order(results.distal.analysis.sig.int$`RLM_DNAmGroup:TF_triplet_stage_wise_adj_pvalue`),]

plots <- plot_interaction_model(
  triplet.results = results.distal.analysis.sig.int[1:20,], 
  dnam = resid.met.cpg.awas, 
  exp = residuals.matched.exp,
  tf.activity.es = rnaseq.tf.es.gsva.ensg
)

# Merge plots into one file 
plots.one.page <- gridExtra::marrangeGrob(plots, nrow = 1, ncol = 1)

ggplot2::ggsave(
  filename = file.path(path.plots, paste0("Distal_top_20_remap_tf.es.gsva_rna_residuals_dnam_residuals",".pdf")),
  plot = plots.one.page,
  width = 11,
  height = 13
)  



#-------------------------------------------------------------------------------
# Regulon analysis, triplets using remap
#-------------------------------------------------------------------------------
library(tidyr)
regulons.brain <- readr::read_tsv(file.path(path.plots, "../data//brain.TFs.gmt.txt"),col_names = F)
regulons.brain$X2 <- NULL
colnames(regulons.brain)[1] <- "tf"
regulons.brain <- regulons.brain %>%  pivot_longer(colnames(regulons.brain)[-1], names_to = "column", values_to = "target")
regulons.brain$column <- NULL


triplet.regulon.ewas <- create_triplet_regulon_based(
  region = hm450.hg19[meta.analysis.cpgs$cpg,],
  tf.target = regulons.brain,
  genome = "hg19",
  TF.peaks.gr = remapCatalog,
  motif.search.window.size = 500
) 
gc()
dim(triplet.regulon.ewas) # 76,203 triplets
triplet.regulon.ewas <- triplet.regulon.ewas %>% dplyr::filter(.data$TF %in% rownames(rnaseq.tf.es.gsva.ensg))
dim(triplet.regulon.ewas) # 75,839   triplets
triplet.regulon.ewas$probeID <- names(hm450.hg19)[match(triplet.regulon.ewas$regionID,make_names_from_granges(hm450.hg19))]


triplet.regulon.ewas$TF %>% unique %>% length # 285
triplet.regulon.ewas$regionID %>% unique %>% length # 3465
triplet.regulon.ewas$target %>% unique %>% length # 6,750

file.regulon <- file.path(path.tables, "ROSMAP_and_remap_regulon_analysis_using_TF_es_gsva_all_triplet.csv")

if (!file.exists(file.regulon)) {
  cores <- 1
  results.regulon.analysis <- 
    triplet.regulon.ewas %>% 
    cor_tf_target_gene(
      exp = residuals.matched.exp,
      tf.activity.es = rnaseq.tf.es.gsva.ensg,
      cores = cores
    ) %>% cor_dnam_target_gene(
      dnam = resid.met.cpg.awas,
      exp = residuals.matched.exp,
      cores = cores,
      filter.results = FALSE, 
      min.cor.estimate = 0.2,
      min.cor.pval = 0.05
    ) %>%  interaction_model(
      dnam = resid.met.cpg.awas,
      exp = residuals.matched.exp,
      tf.activity.es = rnaseq.tf.es.gsva.ensg,
      cores = cores,
      filter.correlated.tf.exp.dnam = FALSE,
      filter.triplet.by.sig.term = FALSE,
      sig.threshold = 0.05,
      fdr = TRUE,
      stage.wise.analysis = TRUE
    ) %>%  stratified_model(
      dnam = resid.met.cpg.awas,
      exp = residuals.matched.exp,
      tf.activity.es = rnaseq.tf.es.gsva.ensg,
      cores = cores
    )
  
  # results.regulon.analysis <- results.regulon.analysis[results.regulon.analysis$probeID %in% cpgs.iqr.higher.0.03,]
  results.regulon.analysis <- results.regulon.analysis %>% 
    add_annot_cpgs() %>%  
    add_percent_zero_q1_q4(dnam = matched.dnam, exp = matched.exp.log2)  %>%
    update_met_IQR(dnam = matched.dnam.with.regions)
  
  
  results.regulon.analysis$RLM_TF_fdr <- p.adjust(results.regulon.analysis$RLM_TF_pvalue,method = "fdr")
  results.regulon.analysis$RLM_DNAmGroup_fdr <- p.adjust(results.regulon.analysis$RLM_DNAmGroup_pvalue,method = "fdr")
  results.regulon.analysis$`RLM_DNAmGroup:TF_fdr` <- p.adjust(results.regulon.analysis$`RLM_DNAmGroup:TF_pvalue`,method = "fdr")
  
  readr:::write_csv(
    x = results.regulon.analysis,
    file = file.regulon
  )
  results.regulon.analysis.sig.fdr.int <- results.regulon.analysis %>%
    dplyr::filter(`RLM_DNAmGroup:TF_fdr` < 0.05)
  results.regulon.analysis.sig.fdr.int <- results.regulon.analysis.sig.fdr.int[order(results.regulon.analysis.sig.fdr.int$`RLM_DNAmGroup:TF_fdr`),]
  
  readr:::write_csv(
    x = results.regulon.analysis.sig.fdr.int,
    file = gsub("all_triplet","sig_fdr_int_triplet",file.regulon)
  )
  
  results.regulon.analysis.sig.fdr.int.with.blood <- merge(results.regulon.analysis.sig.fdr.int, london.blood)
  readr:::write_csv(
    x = results.regulon.analysis.sig.fdr.int.with.blood,
    file = gsub("all_triplet","sig_fdr_int_triplet_with_london_blood",file.regulon)
  )
  
  results.regulon.analysis.sig.int <- results.regulon.analysis %>%
    dplyr::filter(`RLM_DNAmGroup:TF_triplet_stage_wise_adj_pvalue` < 0.05)
  
  results.regulon.analysis.sig.int <- results.regulon.analysis.sig.int[order(results.regulon.analysis.sig.int$`RLM_DNAmGroup:TF_triplet_stage_wise_adj_pvalue`),]
  readr:::write_csv(
    x = results.regulon.analysis.sig.int,
    file = gsub("all_triplet","sig_stage_wise_fdr_int_triplet",file.regulon)
  )
  
  results.regulon.analysis.sig <- results.regulon.analysis %>%
    filter_at(vars(contains("triplet_stage")), any_vars(. < 0.05))
  
  readr:::write_csv(
    x = results.regulon.analysis.sig,
    file = gsub("all_triplet", "sig_any_stage_wise_triplet", file.regulon)
  )
} else {
  results.regulon.analysis <- readr::read_csv(
    file = file.regulon,
    col_types = c(`quant_triplet_stage_wise_adj_pval_metGrp:es.tf` = "n")
  )
  results.regulon.analysis$is.enhancer.EhnancerAtlas.Brain.neuro <- results.regulon.analysis$probeID %in% enhancer.probes
}



plots <- plot_interaction_model(
  triplet.results = results.regulon.analysis.sig.int[1:20,], 
  dnam = resid.met.cpg.awas, 
  exp = residuals.matched.exp,
  tf.activity.es = rnaseq.tf.es.gsva.ensg
)

# Merge plots into one file 
plots.one.page <- gridExtra::marrangeGrob(plots, nrow = 1, ncol = 1)

ggplot2::ggsave(
  filename = file.path(path.plots, paste0("Regulon_top_20_remap_tf.es.gsva_rna_residuals_dnam_residuals",".pdf")),
  plot = plots.one.page,
  width = 11,
  height = 13
)  


plot_in_one_page <- function(results, out){
  plots <- plot_interaction_model(
    triplet.results = results, 
    dnam = matched.dnam.with.regions, 
    exp = residuals.matched.exp,
    genome = "hg19",
    label.dnam = "beta-value",
    label.exp = "residuals",
    tf.activity.es = rnaseq.tf.es.gsva.ensg
  )
  
  # Merge plots into one file 
  plots.one.page <- gridExtra::marrangeGrob(plots, nrow = 1, ncol = 1)
  
  ggplot2::ggsave(
    filename = file.path(path.plots, out),
    plot = plots.one.page,
    width = 11,
    height = 13
  )  
}


plot_in_one_page_resid <- function(results, out){
  plots <- plot_interaction_model(
    triplet.results = results, 
    dnam = resid.met.cpg.awas, 
    exp = residuals.matched.exp,
    genome = "hg19",
    label.dnam = "residuals",
    label.exp = "residuals",
    tf.activity.es = rnaseq.tf.es.gsva.ensg
  )
  
  # Merge plots into one file 
  plots.one.page <- gridExtra::marrangeGrob(plots, nrow = 1, ncol = 1)
  
  ggplot2::ggsave(
    filename = file.path(path.plots, out),
    plot = plots.one.page,
    width = 11,
    height = 8
  )  
}

to.plot <- readxl::read_xlsx("~/TBL Dropbox/Tiago Silva//Draft_TF/DRAFT_TABLES_FIGURES_10-7-2020/Main_Table 3 ROSMAP promoter-distance-based - top10-V2.xlsx", skip = 5)
to.plot <- to.plot[!is.na(to.plot$...2),]
to.plot <- to.plot[,c(1,3,4)]
colnames(to.plot) <- c("probeID","TF_symbol", "target_gene")
to.plot$tf <- map_symbol_to_ensg(to.plot$TF_symbol)
to.plot$target <- map_symbol_to_ensg(to.plot$target_gene)
to.plot$regionID <- MethReg:::make_names_from_granges(hm450.hg19[to.plot$probeID,])
all.results <- plyr::rbind.fill(
  readr::read_csv(file.promoter),
  readr::read_csv(file.distal),
  readr::read_csv(file.regulon)
)


to.plot.aux <- all.results[match(paste0(to.plot$probeID,to.plot$target,to.plot$tf),paste0(all.results$probeID,all.results$target,all.results$TF)),]
plot_in_one_page(results = to.plot.aux, out = "Main_Table 3 ROSMAP promoter-distance-based_beta_values.pdf") 
plot_in_one_page_resid(results = to.plot.aux, out = "Main_Table 3 ROSMAP promoter-distance-based_residuals.pdf") 

plot_in_one_page_resid(results.promoter.analysis.sig.fdr.int, "Promoter_fdr_Fixed_IQR_top_3_remap_tf.es.gsva_rna_residuals_dnam_residuals.pdf") 
plot_in_one_page_resid(results.regulon.analysis.sig.fdr.int[1:20,], "Regulon_fdr_Fixed_IQR_top_20_remap_tf.es.gsva_rna_residuals_dnam_residuals.pdf") 
plot_in_one_page_resid(results.distal.analysis.sig.fdr.int[1:20,], "Distal_fdr_Fixed_IQR_top_20_remap_tf.es.gsva_rna_residuals_dnam_residuals.pdf") 


hm450.hg19 <- MethReg:::get_met_probes_info(genome = "hg19", arrayType = "450k")
load(file.path(path.data,"EnhancerAtlas_v2/Species_enhancer.RData"))
enhancers <- HS %>% dplyr::select(c(
  "Enh_peaks",
  "Astrocyte",
  "hNCC",
  "KELLY",
  "BE2C",
  "ESC_NPC",
  "NGP",
  "SK.N.MC",
  "Cerebellum",
  "Gliobla",
  "SH.SY5Y",
  "SK.N.SH_RA",
  "T98G",
  "Fetal_brain",
  "U87","H54",
  "SK.N.SH",
  "SK.N.MC",
  "NH.A",
  "Macrophage", 
  "ESC_neuron"
))

enhancers <- enhancers[rowSums(enhancers[,-1]) > 0,]
enhancers.gr <- enhancers %>% 
  tidyr::separate(col = "Enh_peaks", into = c("chr","start", "end")) %>% 
  makeGRangesFromDataFrame()

enhancer.probes <- names(subsetByOverlaps(hm450.hg19 + 250, enhancers.gr))

