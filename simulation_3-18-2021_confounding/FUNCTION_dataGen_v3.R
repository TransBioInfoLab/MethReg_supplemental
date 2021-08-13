
# Simulate data
# input:
# data = dataframe, include variables met, rna.tf
# theta = dispersion parameter from neg binomial distribution
# beta.df = beta coefficient, strength of association between target gene and TF

# mu <- 10.59
# theta <- 17.78
# data <- met.tf
# beta.tf <- 0.5

dataGen <- function (data, mu, theta, beta.tf){

  quant.met <-  quantile(data$met,na.rm = TRUE)
  q2.met <- quant.met [2]
  q4.met <- quant.met [4]

  # variable for met.grp
  data$low.met  <- ifelse (data$met < q2.met, 1, 0)
  data$high.met <- ifelse (data$met > q4.met, 1, 0)

  data$met.grp <- ifelse (data$low.met ==1 , "low",
                          ifelse(data$high.met == 1, "high", NA) )

  # simulate targe gene expressions

  data$target.1 <- rnbinom(n = nrow(data),
                           size = theta,
                           mu = mu + beta.tf * data$rna.tf )

  data$target.0 <- rnbinom (n = nrow(data),
                            size = theta,
                            mu = mu) #randomly generated
  
  # normalize data so that target.1 and target.0 have the same mean and sd
  # mean1 <- mean (data$target.1)
  # sd1   <- sd (data$target.1)
  # 
  # mean0 <- mean(data$target.0)
  # sd0   <- sd (data$target.0)
    
  # data$target.1.norm <- (data$target.1 - mean1)
  # 
  # data$target.0.norm <- (data$target.0 - mean0)
  # 
  # data$rna.target <- ifelse (data$low.met == 1, data$target.1.norm, data$target.0.norm)

  data$rna.target <- ifelse (data$low.met == 1, data$target.1, data$target.0)
  
  data
}


### simulate scenario with confounding effect
#   target gene ~ met
#   TF ~ met 

dataGen_confounding <- function (data, mu, theta, beta.met, beta.tf){
  
  quant.met <-  quantile(data$met,na.rm = TRUE)
  q2.met <- quant.met [2]
  q4.met <- quant.met [4]
  
  # variable for met.grp
  data$low.met  <- ifelse (data$met < q2.met, 1, 0)
  data$high.met <- ifelse (data$met > q4.met, 1, 0)
  
  data$met.grp <- ifelse (data$low.met ==1 , "low",
                          ifelse(data$high.met == 1, "high", NA) )
  
  
  # simulate targe gene expressions
  data$rna.target <- rnbinom(n = nrow(data),
                             size = theta,
                             mu = mu + beta.tf * data$rna.tf )
  
  data
}

### simulate data from normal distribution
# the normal distribution has mean = mu, std = sigma

dataGen_normal_distribution <- function (data, mu, sigma, beta.tf){
  
  quant.met <-  quantile(data$met,na.rm = TRUE)
  q2.met <- quant.met [2]
  q4.met <- quant.met [4]
  
  # variable for met.grp
  data$low.met  <- ifelse (data$met < q2.met, 1, 0)
  data$high.met <- ifelse (data$met > q4.met, 1, 0)
  
  data$met.grp <- ifelse (data$low.met ==1 , "low",
                          ifelse(data$high.met == 1, "high", NA) )
  
  # simulate targe gene expressions
  
  data$target.1 <- rnorm (n = nrow(data),
                          mean = mu + beta.tf * data$rna.tf, 
                          sd = sigma)
  
  data$target.0 <- rnorm (n = nrow(data),
                          mean = mu, 
                          sd = sigma) #randomly generated
  
  data$rna.target <- ifelse (data$low.met == 1, data$target.1, data$target.0)
  
  
  data
}



