
setwd("C:/Users/lxw391/Dropbox (BBSR)/Draft_TF/simulation_3-18-2021_confounding")
# path.dropbox <- dir("~",pattern = "Dropbox",full.names = TRUE)
# setwd(file.path(path.dropbox,"Draft_TF/simulation_3-16-2021"))

library (MASS)
library (ggplot2)
library (dplyr)
library (tidyr)
library (pROC)

######################## 1. simulate data  #############################
# (1) file 1_get_data.R
# - filter dnam & exp data
# - rna data: estimate mean and variance,
#             compute mu & theta for negative binomial distribution
#
# (2) simulate_data.R, function simData
# - generate target gene values
#       in low DNAm samples:  target ~ negBin ( beta_tf *TF, theta)
#       in other samples:     target ~ negBin ( mu , theta)
#
# (3) fit models
# - rlm.cont, rlm.binary, lm.cont, lm.binary, main effect models
#
# (4) save resulting p-values
#
########################################################################

source ("./FUNCTION_dataGen_v3.R")
source ("./FUNCTION_fit_models_v2.R")

dna.met.chr21.filtered  <- readRDS("./data/dna.met.chr21.filtered.RDS")
gene.exp.chr21.filtered <- readRDS("./data/gene.exp.chr21.filtered.RDS")

dim (dna.met.chr21.filtered)  # 72 38
dim (gene.exp.chr21.filtered) # 119 38

###########################

n.reps <- 1000

parms.beta.tf <- c(3, 6, 9)

res.all <- data.frame (matrix(nrow = 0, ncol = 9))

set.seed (7262021)

# doParallel::registerDoParallel(cores = 4)
# plyr::a_ply(.data = 1:length(parms.beta.tf),.margins = 1,.fun = function(p){

for (p in 1:length(parms.beta.tf)){
  
         ### generate row numbers in dnam and rnaseq data matrices
        set.seed (7262020 + p)
        row.met <- sample (
          1:nrow(dna.met.chr21.filtered),
          size = n.reps,
          replace = TRUE
        )
        
        set.seed (7262021 + p)
        row.tf <- sample (
          1:nrow(gene.exp.chr21.filtered),
          size = n.reps,
          replace = TRUE
        )
        
        # set.seed (3252021 + p)
        # row.target <- sample (
        #   1:nrow(gene.exp.chr21.filtered),
        #   size = n.reps,
        #   replace = TRUE
        # )
 
    for (r in 1:n.reps){

       ### 1. (a) true positives - generate data  
      
        # met.tf <- data.frame(
        #     met    = dna.met.chr21.filtered [row.met[r] ,] %>% as.numeric,
        #     rna.tf = gene.exp.chr21.filtered[row.tf[r]  ,] %>% as.numeric
        # )
        # 
        # sim.data <- dataGen (
        #     data = met.tf,
        #     mu = 10.59, theta = 17.78,
        #     beta.tf = parms.beta.tf [p]
        # )
        # 
        # #### 1. (b) true positives - fit models 
        # 
        # res.rlm.cont   <- fit.rlm.cont (sim.data)
        # res.rlm.binary <- fit.rlm.binary (sim.data)
        # 
        # res.main       <- fit.main.effect.models.met (sim.data)
        # 
        # res.lm.cont    <- fit.lm.cont (sim.data)
        # res.lm.binary  <- fit.lm.binary(sim.data)
        # 
        # res <- rbind (
        #   res.rlm.binary,
        #   res.rlm.cont,
        #   res.main,
        #   res.lm.cont,
        #   res.lm.binary
        # )
        # 
        # res.one <- data.frame (
        #   beta.tf = parms.beta.tf [p],
        #   rep = r,
        #   row.met = row.met[r],
        #   row.rna = row.tf [r],
        #   res
        # )
        
        
        #### 2. generate data - true negatives (confounding case)
        # - order rna.tf according to -met
       
        met          <- dna.met.chr21.filtered [row.met[r] ,] %>% as.numeric 
        met_ordered  <- sort (met, decreasing = TRUE)
        tf_ordered   <- gene.exp.chr21.filtered[row.tf[r]  ,] %>% as.numeric %>% sort 
         
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
        
        ## test 
        # low  <- subset(data, data$low.met == 1)
        # high <- subset(data, data$low.met == 0)
        # 
        # cor.test (low$rna.tf, low$rna.target, method = "spearman")
        # cor.test (high$rna.tf, high$rna.target, method = "spearman")
        # 
        # plot (low$rna.tf, low$rna.target)
        # plot (high$rna.tf, high$rna.target)
        
        #### 2. (b) true negatives - fit models 
        
        res.rlm.cont   <- fit.rlm.cont (sim.data.trueNeg)
        res.rlm.binary <- fit.rlm.binary (sim.data.trueNeg)
        
        res.main       <- fit.main.effect.models.met (sim.data.trueNeg)
        
        res.lm.cont    <- fit.lm.cont (sim.data.trueNeg)
        res.lm.binary  <- fit.lm.binary(sim.data.trueNeg)
        
        res <- rbind (
          res.rlm.binary,
          res.rlm.cont,
          res.main,
          res.lm.cont,
          res.lm.binary
        )
        
        res.one <- data.frame (
          beta.tf = parms.beta.tf [p],
          rep = r,
          row.met = row.met[r],
          row.rna = row.tf [r],
          res
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

#},.progress = "time",.parallel = TRUE)

#### computing test sizes & power for each model

# add max(pval.binary, pval.cont)

res.pvals <- subset (
    res.all,
    select = c(beta.tf, rep, row.met, row.rna, model, pval)
)

res.pvals.wide <- spread (
    data = res.pvals,
    key = model,
    value = pval
)

# res.pvals.wide$max.rlm.pval <- ifelse (
#     res.pvals.wide$rlm.binary > res.pvals.wide$rlm.cont,
#     res.pvals.wide$rlm.binary,
#     res.pvals.wide$rlm.cont
# )


saveRDS (res.pvals.wide, "./results/res.pvals.wide.RDS")

write.csv (res.pvals.wide, "./results/res.pvals.wide.csv", row.names = FALSE)

# wide to long format
res.pvals.tall <- res.pvals.wide %>% pivot_longer(
    cols = lm.binary:wilcox.main.met.tf,
    names_to = "model",
    values_to = "pval"
)

res.pvals.tall$beta.tf <- as.factor (res.pvals.tall$beta.tf)

res.pvals.tall$sig.pval <- ifelse (res.pvals.tall$pval < 0.05, 1, 0)

saveRDS (res.pvals.tall, "./results/res.pvals.tall.RDS")
