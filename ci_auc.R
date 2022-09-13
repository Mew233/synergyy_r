library(pROC)
library(ggplot2)
mpreds_trans = read.csv('data/predicts_transynergy_liu_DrugComb.txt')

roc.list <- roc(actuals ~ predicts_transynergy_liu, data = mpreds_trans)


ciobj <-ci.se(roc.list, specificities = seq(0, 1, 0.1), conf.level=0.95, boot.n=1000,)

# dat.ci.list <- lapply(ciobj, function(ciobj) 
dat.ci.list <-  data.frame(x = as.numeric(rownames(ciobj)),
             lower = ciobj[, 1],
             upper = ciobj[, 3])

p <- ggroc(roc.list) + theme_minimal() + geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=0.7, color = "grey") + coord_equal()

for(i in 1:3) {
  p <- p + geom_ribbon(
    # data = dat.ci.list[[i]],
    data = dat.ci.list,
    aes(x = x, ymin = lower, ymax = upper),
    fill = i + 1,
    alpha = 0.2,
    inherit.aes = F) 
} 

p
