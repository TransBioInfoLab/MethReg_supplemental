setwd("C:/Users/lxw391/Dropbox (BBSR)/Draft_TF/simulation_3-18-2021_confounding/results")

# setwd("C:/Users/lxw391/Dropbox (BBSR)/Draft_TF/simulation/results")
# path.dropbox <- dir("~",pattern = "Dropbox",full.names = TRUE)
# setwd(file.path(path.dropbox,"Draft_TF/simulation"))
library (MASS)
library (ggplot2)
library (dplyr)
library (tidyr)
library (pROC)

########################################################################
# (1) combine results for rlm.binary.en model & all other models
#
# (2) compute power, test sizes, roc, auc
#
########################################################################

res.pvals.tall <- readRDS("C:/Users/lxw391/Dropbox (BBSR)/Draft_TF/simulation_3-18-2021_confounding/results/res.pvals.tall.RDS")

res.pvals.tall <- subset (res.pvals.tall, !res.pvals.tall$model %in% c("spearman.corr.met.tf", "wilcox.main.met.tf"))


confounding <- aggregate (sig.pval ~ beta.tf + model, data = res.pvals.tall, FUN = "mean") %>% rename (pct.sig.pval = sig.pval)

saveRDS   (confounding, "confound_pct_sig.RDS")
write.csv (confounding, "confound_pct_sig.csv", row.names = FALSE)

pdf ("confounding_scenario.pdf", height = 6, width = 11)

library(ggplot2)
library(colortools)

colors <- c(pals("dream")[2:4], pals("mystery")[5],  
            pals("cheer")[4:5]
            )

p <- ggplot(confounding, aes(x = beta.tf, y = pct.sig.pval, fill =  model)) +
    geom_bar(stat="identity", position = position_dodge()) +
    scale_fill_manual("legend", values = colors[1:8]) +
    labs(y= "proportion of significant p-values", x = "effect size of TF association with target gene expression") +
    theme_bw()
p

dev.off()

############### 

# histgrams
res.pvals.wide.test.size <- subset (res.pvals.wide, beta.tf == 0)

models <- c(
    "lm.binary",
    "lm.cont",
    "wilcox.main.met",
    "rlm.binary",
    "rlm.cont",
    "spearman.corr.met"
    #"max.rlm.pval"
)

pdf ("histgram under null.pdf")

for (i in 1:length(models))
    hist(res.pvals.wide.test.size[, models[i]],
         main = models[i], xlab = " ")
dev.off ()

############################ roc & auc ----------------------------------------------------

models <- c(
    "lm.binary",
    "lm.cont",
    "wilcox.main.met",
    "rlm.binary",
    "rlm.binary.en",
    "rlm.cont",
    "spearman.corr.met"
)


rlmBinary.en_all_pvals$pval <- rlmBinary.en_all_pvals$pval.corrected
rlmBinary.en_all_pvals$sig.pval <- ifelse (rlmBinary.en_all_pvals$pval.corrected < 0.05, 1, 0)

one <- subset (rlmBinary.en_all_pvals, 
                 select = c(beta.tf, rep, row.met, row.rna, model,pval, sig.pval))

two <- res.pvals.tall

res.pvals.tall <- rbind (one, two)



res.pvals.tall$actual <- as.factor(
    ifelse (res.pvals.tall$beta.tf == 0,
            "true.negative",
            "true.positive")
)

auc.all <- data.frame (matrix(nrow = 0, ncol = 2))


lroc <- plyr::alply( 1:length(models), .margins = 1, function(i){
    one <- subset (
        res.pvals.tall,
        model == models[i]
    )
    roc (response = one$actual, one$pval, auc = TRUE)

})
library(pROC)
names(lroc) <- models
plot <- ggroc(lroc) +
    ggpubr::theme_pubr() +
    ggtitle("ROC curve") +
    theme(legend.position = "right") +
    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1),
                 color = "grey",
                 linetype = "dashed"
    ) + labs(color = "Approaches")
ggsave(filename = "ROC_curve_all_models_small.png", plot = plot, height = 3, width = 5)
ggsave(filename = "ROC_curve_all_models_small.pdf", plot = plot, height = 3, width = 5)

colors <- c("red","green", "blue", "purple", "red", "pink","orange", "black")
for (i in 1:length(models)){

    one <- subset (
        res.pvals.tall,
        model == models[i]
    )

    res <- roc (response = one$actual, one$pval, auc = TRUE)

    if(i == 1){
        plot (res, col = colors[i])

    } else {
        plot (res, main = models[i], col = colors[i], add = TRUE)

    }

    one_auc <- data.frame (
        model = models[i],
        auc = res$auc
    )

    auc.all <- rbind (auc.all, one_auc)

}

legend(
    x = 0,
    y = 0.8,
    legend = models,
    title = "Models",
    lty = c(1,1),
    lwd = c(2,2),
    col = colors
)

write.csv (auc.all, "auc.all.csv", row.names = FALSE)



## specificity - does not look good
# spec <- subset(res.pvals.tall, beta.tf == 0 )
# spec$specificity <- spec$accuracy
#
# p <- ggplot(spec, aes(x = model, y = specificity)) +
#     geom_boxplot()
#
# p


# FDRs = wrongly rejected / # of rejected
#      = true negative & sig / sig

# pdf ("fdr.pdf")
# 
# rejected <- subset (res.pvals.tall, sig.pval == 1)
# 
# rejected$false.rejected <- ifelse (rejected$beta.tf == 0, 1, 0)
# 
# fdr <- aggregate (false.rejected ~  model2, data = rejected, FUN = "mean")
# 
# fdr <- fdr[fdr$model2 != "max.rlm",]
# 
# p <- ggplot(fdr, aes(x = model2, y = false.rejected)) +
#     geom_bar(stat="identity", fill = "lightblue") +
#     geom_text(aes(label= round(false.rejected, digits = 3)), vjust = 1.6, size=3.5) +
#     theme_bw() +
#     labs(y= "false discovery rate", x = " ") +
#     theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# p
# 
# dev.off()
