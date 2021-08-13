
# 2. Estimate parameters from neg binomial distribution
# - fit negative binomial distribution to rna data
# - estimate mu and theta -> histogram, median

library(MASS)

dat <- readRDS("./data/gene.exp.chr21.filtered.RDS")

# these didn't work because data did not follow negative
# binomial distribution for this small sample

# for (i in 1:10){
#
#     onerow <- data.frame (rna = dat[6,])
#     onerow$rna <- trunc(onerow$rna)
#     hist(onerow$rna)
#
#     f <- glm.nb( rna ~ 1 , data = onerow, link = identity)
# }
#
# summary(fm)
#
# fm <- glm.nb(Days ~ 1, data = quine, link = identity)
# summary(fm)


med.mean <- median(apply (dat, 1, mean))  # median of mean gene expression

med.std <- median (apply (dat, 1, sd))  # median of standard deviations

med.var <- med.std^2

## for neg bin distribution

mu <- med.mean  #mu = 10.59295

theta <- (mu*mu)/(med.var - mu)  # 17.7809

# check
mu + mu*mu/theta  # same as med.var
