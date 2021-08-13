library(ELMER)
library(dplyr)
library(SummarizedExperiment)


setwd("~/Dropbox (BBSR)/Draft_TF/simulation_3-18-2021_confounding/")
# path.dropbox <- dir("~",pattern = "Dropbox",full.names = TRUE)
# setwd(file.path(path.dropbox,"Draft_TF/simulation_3-16-2021"))

library (MASS)
library (ggplot2)
library (dplyr)
library (tidyr)
library (pROC)
source ("./FUNCTION_dataGen_v3.R")
source ("./FUNCTION_fit_models_v2.R")

dna.met.chr21.filtered  <- readRDS("./data/dna.met.chr21.filtered.RDS")
gene.exp.chr21.filtered <- readRDS("./data/gene.exp.chr21.filtered.RDS")

dim (dna.met.chr21.filtered)  # 72 38
dim (gene.exp.chr21.filtered) # 119 38

###########################

n.reps <- 1000

parms.beta.tf <- c(0:10)
#parms.beta.tf <- 0

res.all <- data.frame (matrix(nrow = 0, ncol = 9))

data <-  plyr::alply(.data = parms.beta.tf,.margins = 1,.fun = function(p){
  
  row.met <- sample (
    1:nrow(dna.met.chr21.filtered),
    size = n.reps,
    replace = TRUE
  )
  
  row.tf <- sample (
    1:nrow(gene.exp.chr21.filtered),
    size = n.reps,
    replace = TRUE
  )
  
  met          <- dna.met.chr21.filtered [row.met[1] ,] %>% as.numeric 
  met_ordered  <- sort (met, decreasing = TRUE)
  tf_ordered   <- gene.exp.chr21.filtered[row.tf[1]  ,] %>% as.numeric %>% sort 
  
  met.tf <- data.frame (
    met = met_ordered, 
    rna.tf = tf_ordered 
  )
  
  met.tf <- subset (met.tf, met.tf$rna.tf > 0 )
  
  # add variables low.met, high.met, and met.grp only, 
  # not using parms mu, theta, beta.tf in this function
  sim.data.trueNeg <- dataGen_confounding(
    data = met.tf, 
    mu = 10.59, theta = 17.78,
    beta.tf = parms.beta.tf [p]
  )
  
  sim.data.trueNeg
})




devtools::load_all("~/Documents/packages/MethReg/")


library(ggplot2)
list.plots <- lapply(2:length(data), function(beta){
  df <- data[[beta]]
  df$TF.group <- NA
  df[df$rna.tf < quantile(df$rna.tf)[2],]$TF.group <- "1 - Low"
  df[df$rna.tf > quantile(df$rna.tf)[4],]$TF.group <- "2 - High"
  df$met.grp <- ifelse(df$met.grp == "low","DNAm low","DNAm high")
  df$met.grp <- factor(df$met.grp, levels = c("DNAm low","DNAm high"))
  
  get_scatter_plot_results(
    df[!is.na(df$met.grp),],
    x = "rna.tf",
    y = "rna.target",
    facet.by = "met.grp",
    color = "black",
    ylab = "Target gene expression",
    xlab = "TF expression"
  ) + ggtitle(paste0("Beta = ", beta - 1))
})

ggarrange(plotlist = list.plots)