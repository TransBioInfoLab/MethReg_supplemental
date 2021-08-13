# evaluate test size of corrected p-values (rlm.binary model) computed from empirical null distribution

setwd("C:/Users/lxw391/Dropbox (BBSR)/Draft_TF/simulation_3-16-2021")

library (dplyr)
##############################################################################################
# (1) simulate 1000 true postive genes for each scenario parms.beta.tf <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
# - fit rlm.binary model
# - compute z-scores

# (2) simulate 10,000 null genes for parm.beta.tf = 0
# - fit rlm.binary model
# - compute z-scores

# (3) combine data from (1) and (2)
# - locfdr to compute location and scale parameter & corrected p-values

# (4) estimate test size
########################################################################################################

# res.pvals.tall <- readRDS("C:/Users/lxw391/Dropbox (BBSR)/Draft_TF/simulation/res.pvals.tall.RDS")
# 
# res.pvals.tall$beta.tf <- as.numeric (as.character(res.pvals.tall$beta.tf))
# 
# all_pvals <- subset (res.pvals.tall, model == "rlm.binary")

source ("./FUNCTION_dataGen_v2.R")
source ("./FUNCTION_fit_models_v2.R")

library(MethReg)

dna.met.chr21.filtered  <- readRDS("./data/dna.met.chr21.filtered.RDS")
gene.exp.chr21.filtered <- readRDS("./data/gene.exp.chr21.filtered.RDS")

dim (dna.met.chr21.filtered)  # 72 38
dim (gene.exp.chr21.filtered) # 119 38

dna.met <- dna.met.chr21.filtered
gene.exp <- gene.exp.chr21.filtered

#### using full dataset

# load("C:/Users/lxw391/Dropbox (BBSR)/PanCancer/useCase/data/COAD_READ_matched_data_dnam_rna_cnv.rda")
# 
# library(SummarizedExperiment)
# 
# dnam <- filter_dnam_by_quant_diff(
#   dnam = dnam,
#   diff.mean.th = 0.0,
#   cores = 5
# )
# 
# dna.met <- assays(dnam)[[1]] [, 1:38]
# 
# # Remove genes that has 0 for more than half of the samples
# rnaseq.filtered <- rnaseq %>%
#   filter_genes_zero_expression(
#     max.samples.percentage = 0.3
#   )
# 
# gene.exp <- rnaseq.filtered[, 1:38]

#--------------------------------- (1) true positives ----------------------------------- 

n.reps <- 1000

set.seed (7262020)
row.met <- sample (
  1:nrow(dna.met),
  size = n.reps,
  replace = TRUE
)

set.seed (7262021)
row.tf <- sample (
  1:nrow(gene.exp),
  size = n.reps,
  replace = TRUE
)

parms.beta.tf <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

res.all <- data.frame (matrix(nrow = 0, ncol = 9))

set.seed (7262021)


for (p in 1:length(parms.beta.tf)){
  
  for (r in 1:n.reps){
    
    met.tf <- data.frame(
      met    = dna.met [row.met[r] ,] %>% as.numeric,
      rna.tf = gene.exp[row.tf[r]  ,] %>% as.numeric
    )
    
    # sim.data <- dataGen_normal_distribution (
    #   data = met.tf,
    #   mu = 20, sigma = 1,
    #   beta.tf = parms.beta.tf [p]
    # )
    
    sim.data <- dataGen (
      data = met.tf,
      mu = 10.59, theta = 17.78,
      beta.tf = parms.beta.tf [p]
    )
    
    res.rlm.binary <- fit.rlm.binary.effectSize (sim.data)

    res.one <- data.frame (
      beta.tf = parms.beta.tf [p],
      rep = r,
      row.met = row.met[r],
      row.rna = row.tf [r],
      res.rlm.binary
    )
    
    ########## plots, diagnositics
    
    ## plots
    
    # hist(sim.data$met)
    # hist (sim.data$rna.tf)
    # hist (sim.data$rna.target)
    # table (sim.data$low.met)
    #
    # plot.low.others (sim.data)

    res.all <- rbind (res.all, res.one)
    
  }
}

#-------------------------- (2) true negatives - null genes -------------------------------------
n.reps <- 10000

set.seed (7262020)
row.met <- sample (
  1:nrow(dna.met),
  size = n.reps,
  replace = TRUE
)

set.seed (7262021)
row.tf <- sample (
  1:nrow(gene.exp),
  size = n.reps,
  replace = TRUE
)

parms.beta.tf <- 0

res.null <- data.frame (matrix(nrow = 0, ncol = 9))

for (p in 1:length(parms.beta.tf)){
  
  for (r in 1:n.reps){
    
    set.seed (7262021 + r)
    
    met.tf <- data.frame(
      met    = dna.met [row.met[r] ,] %>% as.numeric,
      rna.tf = gene.exp[row.tf[r]  ,] %>% as.numeric
    )
    
    # sim.data <- dataGen_normal_distribution (
    #   data = met.tf,
    #   mu = 20, sigma = 1,
    #   beta.tf = parms.beta.tf [p]
    # )
    
    sim.data <- dataGen (
      data = met.tf,
      mu = 10.59, theta = 17.78,
      beta.tf = parms.beta.tf [p]
    )
    
    res.rlm.binary <- fit.rlm.binary.effectSize (sim.data)
    
    
    res.one <- data.frame (
      beta.tf = parms.beta.tf [p],
      rep = r,
      row.met = row.met[r],
      row.rna = row.tf [r],
      res.rlm.binary
    )
    
    ########## plots, diagnositics
    
    ## plots
    
    # hist(sim.data$met)
    # hist (sim.data$rna.tf)
    # hist (sim.data$rna.target)
    # table (sim.data$low.met)
    #
    # plot.low.others (sim.data)
    
    
    
    res.null <- rbind (res.null, res.one)
    
  }
}

# hist(res.null$pval)
# sum(res.null$pval <0.05)/length(res.null$pval)

#saveRDS (res.null, "res.null.RDS")

#--------------------------------- (4) Efron's locfdr approach  ------------------------------

rlmBinary_corrected_pval <- function (zvalue){
  
  require (locfdr)
  
  m <- median (zvalue, na.rm = TRUE)
  
  zvalue.centered <- zvalue - m
  
  w <- data.frame( locfdr(zvalue.centered)$fp0 )
  
  delta <- w["mlest", "delta"]  #location parameter
  
  sigma <- w["mlest", "sigma"]  #scale parameter
  
  print (paste0("delta: ", delta))
  print (paste0("sigma: ", sigma))
  
  zvalue.corrected <- (zvalue - m - delta)/sigma
  
  pvalue <- 1- pnorm(zvalue.corrected)
  
}

all.test.size <- data.frame (matrix(nrow = 0, ncol = 2))

all.pvals <- data.frame (matrix(nrow=0, ncol=12))


for (i in 1:9){
  
  truePos  <- subset (res.all, beta.tf == i)
  
  res.both <- rbind (truePos, res.null)
  
  
  ###### --------------- computing corrrected p-values ------------------------
  
  # take out missing values
  res.both <- res.both[!is.na(res.both$tValue) ,]
  
  # maximum zvalue when it's not Inf
  noInf <- subset (res.both, !is.infinite(qnorm(pt(res.both$tValue, res.both$df))))
  
  max <- max(qnorm(pt(noInf$tValue, noInf$df)))
  
  # set to max if Inf
  res.both$zvalue <- ifelse (
    is.infinite(qnorm(pt(res.both$tValue, res.both$df))),
    max,
    qnorm(pt(res.both$tValue, res.both$df))
  )
  
  res.both$pval.corrected <- rlmBinary_corrected_pval (res.both$zvalue)
  
  
  ### compute test size -----------------------------
  res.test.size <- subset (res.both, beta.tf == 0)
  
  test.size <- sum(res.test.size$pval.corrected < 0.05)/length(res.test.size$pval.corrected)
  
  # test size for one scenario
  
  one.test.size <- data.frame (
    beta.tf = i, 
    test.size = test.size
  )
  
  ### output dataset
  
  all.test.size <- rbind (all.test.size, one.test.size)
  
  all.pvals <- rbind (all.pvals, res.both)
  
}

all.test.size$beta.tf <- as.character(all.test.size$beta.tf)

avg <- data.frame ( 
  beta.tf = "avg", 
  test.size = mean (all.test.size$test.size)
)

all <- rbind (all.test.size, avg)

all.pvals$model = "rlm.binary.en"


saveRDS (all, "rlmBinary.en_all_test_size.RDS")

saveRDS (all.pvals, "rlmBinary.en_all_pvals.RDS")

