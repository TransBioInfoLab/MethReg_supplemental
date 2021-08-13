
setwd("C:/Users/lxw391/Dropbox (BBSR)/Draft_TF/simulation")

library(MethReg)

# get datasets ready

### (1) filter data
#  - dnam: require dnam.q4 - dnam.q1 > 0.2
#  - rna: require pct.zeros in rna < 0.25
#  - rna: require fold change for q4/q1 > 1.5

# DNAm
data("dna.met.chr21")
dna.met.chr21.filtered <- filter_regions_by_mean_quantile_difference(dna.met.chr21)
saveRDS (dna.met.chr21.filtered, file = "./data/dna.met.chr21.filtered.RDS" )

# RNA
filter_genes_zero_expression <- function(exp, max.samples.percentage = 0.25){
    genes.keep <- (rowSums(exp == 0) / ncol(exp) <= max.samples.percentage) %>% which %>% names
    message("Removing ", nrow(exp) - length(genes.keep), " out of ", nrow(exp), " genes")
    exp[genes.keep,]
}

data("gene.exp.chr21")
genes.few.zeros <- filter_genes_zero_expression (gene.exp.chr21)
gene.exp.chr21.filtered <- filter_genes_by_quantile_mean_fold_change(genes.few.zeros)
saveRDS (gene.exp.chr21.filtered, file = "./data/gene.exp.chr21.filtered.RDS")

#### (2) rna data:
#   - estimate mean and variance,
#   - compute mu & theta for negative binomial distribution

med.mean <- median(apply (gene.exp.chr21.filtered, 1, mean))  # median of mean gene expression

med.std <- median (apply (gene.exp.chr21.filtered, 1, sd))  # median of standard deviations

med.var <- med.std^2

## for neg bin distribution

mu <- med.mean  #mu = 10.59295

# variance = mu*mu/theta + mu
theta <- (mu*mu)/(med.var - mu)  # 17.7809

# check
mu + mu*mu/theta  # same as med.var
