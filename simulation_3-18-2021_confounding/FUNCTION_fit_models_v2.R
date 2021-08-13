
####################################################
# input:
#    - data has variables for met, rna.tf, met.grp
# output:
#    - model, pval
#####################################################

#data <- sim.data


### fit lm.continuous model

fit.lm.cont <- function (data){

    res.lm.cont <- lm (
        rna.target ~ met + rna.tf + rna.tf * met,
        data = data
    ) %>% summary %>% coef %>% data.frame

    data.frame (
        model = "lm.cont",
        pval = res.lm.cont ["met:rna.tf", "Pr...t.." ]
    )
}

### fit lm.binary model

fit.lm.binary <- function (data){
  
    data <- data[!is.na(data$met.grp) ,]

    res.lm.binary <- lm (
        rna.target ~ met.grp + rna.tf + rna.tf * met.grp,
        data = data
    ) %>% summary %>% coef %>% data.frame

    data.frame (
        model = "lm.binary",
        pval  = res.lm.binary ["met.grplow:rna.tf", "Pr...t.." ]
    )
}

### robust linear models, continuous met

fit.rlm.cont <- function (data){

    res.rlm <- rlm (
        rna.target ~ met + rna.tf + rna.tf * met,
        data = data,
        psi = psi.bisquare,
        maxit = 100) %>% summary %>% coef %>% data.frame

    res.rlm$pval <- 2 * (1 - pt( abs(res.rlm$t.value), df = nrow(data) - 4 ) )

    data.frame (
        model = "rlm.cont",
        pval = res.rlm["met:rna.tf" , "pval"]
    )
}

### robust linear model, binary met - output only pvalues

fit.rlm.binary <- function (data){

    require (MASS)

    data <- data[complete.cases(data) ,]
    
    #data$met.grp <- factor(data$met.grp)

    res.rlm.binary <- rlm (
        rna.target ~ met.grp + rna.tf + rna.tf * met.grp,
        data = data,
        psi = psi.bisquare,
        maxit = 1000) %>% summary %>% coef %>% data.frame

    res.rlm.binary$pval <- 2 * (1 - pt( abs(res.rlm.binary$t.value), df = nrow(data) - 4 ) )

    data.frame (
        model = "rlm.binary",
        pval = res.rlm.binary["met.grplow:rna.tf" , "pval"]
    )
}


### robust linear model, binary met - output pvals + test stat

fit.rlm.binary.effectSize <- function (data){

    require (MASS)

    data <- data[!is.na(data$met.grp) ,]

    res.rlm.binary <- rlm (
        rna.target ~ met.grp + rna.tf + rna.tf * met.grp,
        data = data,
        psi = psi.bisquare,
        maxit = 1000) %>% summary %>% coef %>% data.frame

    res.rlm.binary$pval <- 2 * (1 - pt( abs(res.rlm.binary$t.value), df = nrow(data) - 4 ) )

    data.frame (
        model    = "rlm.binary",
        estimate = res.rlm.binary["met.grplow:rna.tf" , "Value"],
        stdErr   = res.rlm.binary["met.grplow:rna.tf" , "Std..Error"],
        tValue   = res.rlm.binary["met.grplow:rna.tf" , "t.value"],
        df       = nrow(data) - 4,
        pval     = res.rlm.binary["met.grplow:rna.tf" , "pval"]
    )
}

### main effect model:
#   (1) target.rna ~ DNAm;  cor (target.rna, dnam)
#   (2) target.rna ~ tf;    cor (target.rna, tf)

fit.main.effect.models <- function (data){

    # predictor = met
    # lm.main.met <- lm (
    #     rna.target ~ met,
    #     data = data
    # ) %>% summary %>% coef %>% data.frame
    #
    # pval.lm.main.met <- data.frame (model = "lm.main.met",
    #                                 pval = lm.main.met ["met", "Pr...t.."])

    # predictor = tf
    lm.main.tf <- lm (
        rna.target ~ rna.tf,
        data = data
    ) %>% summary %>% coef %>% data.frame

    pval.lm.main.tf <- data.frame (
        model = "lm.main.tf",
        pval = lm.main.tf ["rna.tf", "Pr...t.."]
    )

    cor.main.tf <- cor.test (data$rna.target, data$rna.tf, method = "spearman", exact = FALSE)
    pval.cor.main.tf <- data.frame (
        model = "spearman.corr.tf",
        pval  = cor.main.tf$p.value
    )

    rbind (pval.lm.main.tf, pval.cor.main.tf)
}


### main effect model:
#   (1) target.rna ~ DNAm (high/low)
#   (2) cor (target.rna, dnam)

fit.main.effect.models.met <- function (data){
    require(rstatix)
    
    wilcox.main.met  <- wilcox_test (rna.target ~ met.grp, data = data, exact = TRUE) %>% data.frame
    
    pval.wilcox.main.met <- data.frame (
        model = "wilcox.main.met", 
        pval  = wilcox.main.met[, "p"]
    ) 
    
    cor.main.met <- cor.test (data$rna.target, data$met, method = "spearman", exact = FALSE)
    pval.cor.main.met <- data.frame (
        model = "spearman.corr.met",
        pval  = cor.main.met$p.value
    )
    
    ## correlating met with tf
    wilcox.main.met.tf  <- wilcox_test (rna.tf ~ met.grp, data = data, exact = TRUE) %>% data.frame
    
    pval.wilcox.main.met.tf <- data.frame (
        model = "wilcox.main.met.tf", 
        pval  = wilcox.main.met.tf[, "p"]
    ) 
    
    cor.main.met.tf <- cor.test (data$rna.tf, data$met, method = "spearman", exact = FALSE)
    pval.cor.main.met.tf <- data.frame (
        model = "spearman.corr.met.tf",
        pval  = cor.main.met.tf$p.value
    )
 
    rbind (pval.wilcox.main.met,    pval.cor.main.met, 
           pval.wilcox.main.met.tf, pval.cor.main.met.tf)
}


# summary (lm (rna.target ~ met,
#              data = sim.data))
#
# cor.test (sim.data$rna.target, sim.data$met, method = "spearman")
# cor.test (sim.data$rna.target, sim.data$rna.tf, method = "spearman")
#
# plot (sim.data$met, sim.data$rna.target)


