#https://www.nature.com/articles/ncomms7921
#https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0269570

library(msigdbr)
library(igraph)
library(clusterProfiler)
library(enrichplot)
library(ggridges)
library(dplyr)
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggupset)
library("RColorBrewer")
library(reshape2)
library(readgmt)

hm.palette <- colorRampPalette(rev(brewer.pal(9, 'YlOrRd')), space='Lab')

#使用pid,c2,stautlab for GSEA
sigDB_t2g = as_tibble(fread("data/SignatureDB_030920.csv")) %>%
  dplyr::select(Signature,symbol)
colnames(sigDB_t2g) = c("gene_set","gene")
C2_t2g <- read_gmt("data/c2.cp.pid.v7.5.1.symbols.gmt",tidy = TRUE) #c2.cp.pid.v7.5.1.symbols.gmt, c2.cp.kegg.v7.3.symbols.gmt
Reactome_t2g <- read_gmt("data/c2.cp.reactome.v7.5.1.symbols.gmt",tidy = TRUE)
# C2_t2g = rbind(C2_t2g,Reactome_t2g,sigDB_t2g)
C2_t2g = rbind(C2_t2g,sigDB_t2g)

#输入shapley, 正确预测並且为协同
#shap = as_tibble(fread("data/multitask_shap.csv"))
shap = as_tibble(fread("data/transynergy_shap.csv"))
syn = shap %>% filter(actuals == 1)
# %>% slice_min(proba, n = 50)
#include drugs in synergistic pairs
#shap %>% filter(Drug2 %in% c(pull(syn, Drug1),pull(syn, Drug2)))
nonsyn = shap %>% filter(actuals == 0) %>% slice_min(proba, n = 100)
syn = rbind(syn,nonsyn)


#调换drug1,drug2的顺序
syn= syn %>% mutate(key = paste0(pmin(Drug1, Drug2), pmax(Drug1, Drug2), sep = ""))
keys = syn %>% mutate(key = paste0(pmin(Drug1, Drug2), pmax(Drug1, Drug2), sep = "")) %>% distinct(key)
# keys = c()
# for (i in  1:dim(syn)[1]) {
#   keys = append(keys,syn[i,]$Drug1)
#   key = syn[i,]$Drug2
#   print(key %in% keys)
#   if (key %in% keys) {
#     c2swap <- c("Drug1", "Drug2")
#     syn[i,][c2swap] <- syn[i,][rev(c2swap)]
#   }else{
#     next
#   }
# }
#畫圖: 统计有多少drug 和 drug combination
# ggplot(syn, aes_(x = ~key)) + geom_bar(fill='azure4') +  
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

tp = select(syn, V1, Drug2)
colnames(tp) = c("V1","Drug1")
syn_wide = as.data.frame(rbind(select(syn, V1,Drug1),tp))
syn_wide %>% count(Drug1,sort = TRUE)
#fig1, 畫每個drug出現多少次
syn_wide <- within(syn_wide, 
                   Drug1 <- factor(Drug1, 
                                      levels=names(sort(table(Drug1), 
                                                        decreasing=TRUE))))
ggplot(syn_wide, aes_(x=~Drug1)) + geom_bar(fill='azure4') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_continuous(limits = c(0, 15)) +
  labs(title  = "Single drug apperances")

## fig2, Pazopanib hydrochloride出現次數最多

syn1 = syn_wide %>% filter(Drug1 %in% c('Pazopanib hydrochloride'))
syn2 = syn_wide %>% filter(V1 %in% syn1$V1)

syn_wide = as.data.frame(syn2[-4,])
syn_wide %>% count(Drug1,sort = TRUE)

res = res <- tibble::tibble(Description = split(syn_wide[,2], syn_wide[,1]))

ggplot(res, aes_(x = ~Description)) + geom_bar(fill='azure4')+
  #theme_dose(font.size = 12) +
  xlab(NULL) + ylab(NULL) +
  ggupset::scale_x_upset(order_by = "degree") +
  labs(title  = "Drug2 in the combo of Pazopanib hydrochloride")
ggsave("saved.pdf")

#把每个gsea结果都保存
list_gses <- vector('list', dim(syn)[1])
for(i in  1:dim(syn)[1]) {
  shap1 = syn[i,]
  pairs_exp = shap1 %>% dplyr::select(c(2:(2+2750*3-1)))
  #normalization by z-score
  # pairs_exp=t(scale(t(pairs_exp)))
  pairs_exp = shap1 %>% dplyr::select(c(5501:(5501+2750-1))) 
  
  original_gene_list = t(pairs_exp)
  gene_list<-na.omit(original_gene_list)

  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  names(gene_list) = rownames(original_gene_list)
  
  gse <- GSEA(gene_list, 
               TERM2GENE = C2_t2g,#sigDB_t2g, C2_t2g
               nPerm = 100,
               minGSSize = 30,
               maxGSSize = 500,
               pvalueCutoff = 0.05,
               verbose = TRUE,
               seed=TRUE,
               pAdjustMethod = "none")
  print(i)
  list_gses[[i]] <- gse

}

##utilities
list2df <- function(inputList) {
  # ldf <- lapply(1:length(inputList), function(i) {
  ldf <- lapply(seq_len(length(inputList)), function(i) {
    data.frame(
      categoryID=rep(names(inputList[i]),length(inputList[[i]])),
      Gene=inputList[[i]])
  })
  
  do.call('rbind', ldf)
}

#extract gene sets
my_df_list = list()
for(i in 1:length(list_gses)){
  print(i)
  geneSets <- geneInCategory(list_gses[[i]])
  # 有些sample没有map到pathway
  if(length(geneSets) > 0){
    d_ <- list2df(geneSets)
    #enrichment score
    y_ = as.data.frame(list_gses[[i]])[,c(1,4,5)]
    d_2 = merge(d_, y_, by = "CustomerId",by.x ='categoryID', by.y='ID')
    #normalized shapley score
    y_2 = as.data.frame(gse@geneList)
    y_2$ID = rownames(y_2)
    d = merge(d_2, y_2, by = "CustomerId",by.x ='Gene', by.y='ID')
    
    d['V1'] = rep(paste0(syn[i,]['Drug1'],"_",syn[i,]['Drug2'],"_",syn[i,]['DepMap_ID']),dim(d)[1])
    d['actuals'] = rep(paste0(syn[i,]['actuals']),dim(d)[1])
    d['proba'] = rep(paste0(syn[i,]['proba']),dim(d)[1])
    d['Drug1'] = rep(syn[i,]['Drug1'],dim(d)[1])
    d['Drug2'] = rep(syn[i,]['Drug2'],dim(d)[1])
    d['DepMap_ID'] = rep(syn[i,]['DepMap_ID'],dim(d)[1])
    my_df_list[[i]] = d
  }
}

long_hp = reduce(my_df_list, bind_rows)
# from long to wide
# 只保留pathway enriched
#畫單獨的cell line,"ACH-000267" "ACH-000817" "ACH-000551"
long_hp_pathway = long_hp[!duplicated(long_hp[,c('categoryID','V1','actuals')]),]
# long_hp_pathway1 = long_hp_pathway %>% filter(DepMap_ID %in% c("ACH-000551"))
long_hp_pathway1 = long_hp_pathway
#normalized enrichement score
wide_hp = dcast(long_hp_pathway1, V1 ~ categoryID,value.var="NES")
wide_hp[is.na(wide_hp)] <- 0
wide_hp <- data.frame(wide_hp[,-1], row.names=wide_hp[,1])
#fig3, 畫enriched pathway 的distribution
ggplot(long_hp_pathway1, aes(x = NES, group =factor(V1), color = factor(actuals))) + geom_density(alpha = 0.05)+
  labs(title="enriched score distribution of samples", x ="NES", y="", color ="synergistic")


annotation = long_hp_pathway1[!duplicated(long_hp_pathway1[,c('V1','actuals')]),] %>%
    select(V1,actuals,DepMap_ID,Drug1,Drug2)
v2 = c()
for (i in  1:dim(annotation)[1]) {
  if (annotation[i,]$Drug2 %in% c("Pazopanib hydrochloride")) {#Pazopanib hydrochloride
    c2swap <- c("Drug1", "Drug2")
    annotation[i,][c2swap] <- annotation[i,][rev(c2swap)]
    v2 <- c(v2,paste0(annotation[i,]['Drug2']))
  }
  #print(i)
  else if (annotation[i,]$Drug1 %in% c("Pazopanib hydrochloride")){
    v2 <- c(v2,paste0(annotation[i,]['Drug2']))
  }
  else {
    v2 <- c(v2,"")
  }
}

annotation['Drug2_of_Pazopanib'] = v2
n = pull(annotation,V1)
annotation = annotation %>% select(actuals,DepMap_ID,Drug2_of_Pazopanib)
rownames(annotation) = n

# fig4, 熱圖
out <- pheatmap::pheatmap(wide_hp,fontsize_row=2, fontsize_col=4, annotation_row = annotation,annotation_legend=F,
                   filename="test.pdf") 
colnames(wide_hp[,out$tree_col[["order"]]])

#fig5, 畫韋恩圖overlap genes in pathways, by interactive venn
cat(paste(c(geneSets[1]),sep="\n"))


#fig6, 畫shapley value distributions of core enriched genes for GSEA enriched categories
nameofgeneset = "scGC_cluster1_LZ"
long_hp1 = long_hp %>% filter(categoryID %in% c(nameofgeneset))
genes_in_pathway = long_hp1$Gene
syn_lg = syn %>% dplyr::select(c(5501:(5501+2750-1))) %>% select(genes_in_pathway)

syn_lgks = syn %>% dplyr::select(c(5501:(5501+2750-1))) %>% select(genes_in_pathway)
syn_lgks$actuals = syn$actuals
#計算ks test
syn_lgks.df = syn_lgks %>% 
  group_by(actuals) %>%
  summarise(across(everything(), mean))
ks.test = ks.test(as.numeric(syn_lgks.df[,-1][1,]), as.numeric(syn_lgks.df[,-1][2,]))
print(ks.test)

dist_df_list = list()
for(i in  1:dim(syn_lg)[1]) {
  original_gene_list = t(syn_lg[i,])
  dist_d = as.data.frame(original_gene_list)
  dist_d['rowval'] = rep(paste0(syn[i,]['actuals']),dim(dist_d)[1])
  dist_d['sample'] = rep(paste0(i),dim(dist_d)[1])
  
  colnames(dist_d) = c("V1","rowval","sample")
  dist_df_list[[i]] = dist_d
  #disp = ggplot(as.data.frame(original_gene_list), aes(x=V1,..scaled..)) + geom_density(alpha=.2, fill="#FF6666") 
  }
gs2val.df = reduce(dist_df_list, bind_rows)

ggplot(gs2val.df, aes(x = V1, group =factor(sample), color = factor(rowval))) + stat_ecdf(size = 0.2) +
  labs(title  = nameofgeneset, x ="shapley", y="", color ="synergistic") + 
  annotate("text", x = 0.003, y = 0.25, label = paste0("KS.test,D-stat:",round(ks.test$statistic,4), "\n","p-value:",round(ks.test$p.value,4)))
#+ geom_density(alpha = 0.3)// stat_ecdf(size = 0.2)
ggsave("mtcars.pdf")

# ggplot(gs2val.df, aes(x=V1, y=factor(sample),color = factor(rowval))) + ggridges::geom_density_ridges()
# ggsave("mtcars.pdf")





#LR: =====================================================================
lg = dcast(long_hp_pathway, V1+actuals ~ categoryID,value.var="NES")
# lg = lg %>%
#   purrr::discard(~sum(is.na(.x))/length(.x)* 100 >=50)
lg[is.na(lg)] <- 0
lg <- data.frame(lg[,-1], row.names=lg[,1])
lg$actuals <- as.factor(lg$actuals)

library(glmnet)
#lasso
# mylogit <- glmnet(y=lg$actuals, x=as.matrix(lg[,-1]), alpha=1, family = "binomial")
cvfit = cv.glmnet(y=lg$actuals, x=as.matrix(lg[,-1]), family = "binomial", alpha=1, type.measure = "class")
coef(cvfit, s="lambda.min")

# 只保留gene enriched
# long_hp_gene = long_hp[!duplicated(long_hp[,c('Gene','V1')]),]
# wide_hp_ge = dcast(long_hp_gene, V1 ~ Gene,value.var="gse@geneList")
# wide_hp_ge[is.na(wide_hp_ge)] <- 0
# wide_hp_ge <- data.frame(wide_hp_ge[,-1], row.names=wide_hp_ge[,1])
# pheatmap::pheatmap(wide_hp_ge,fontsize_col=4, filename="test.pdf")
# axis(1,cex.axis=2)

# logistic regression to distinguish syn vs non-syn
# long_hp1 = long_hp %>% filter(categoryID %in% c("REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION","REACTOME_SIGNALING_BY_GPCR"))
long_hp1 = long_hp
long_hp_gene = long_hp1[!duplicated(long_hp1[,c('Gene','V1','actuals')]),]
# wide_hp_ge = dcast(long_hp_gene, actuals+V1 ~ Gene,value.var="gse@geneList")
# wide_hp_ge = wide_hp_ge %>% `rownames<-`(.[,2]) %>% select(-V1)
# wide_hp_ge$actuals <- as.numeric(wide_hp_ge$actuals)
#
# #removes columns with an x% of NAs(50%)
# wide_hp_ge = wide_hp_ge %>%
#   purrr::discard(~sum(is.na(.x))/length(.x)* 100 >=50)
# # 这里应该选择原本的而不是na
# wide_hp_ge[is.na(wide_hp_ge)] <- 0

# 选取genes做logistic regression
genes_in_pathway = long_hp_gene$Gene
syn_lg = syn %>% dplyr::select(c(5501:(5501+2750-1))) %>% select(genes_in_pathway)
#syn_lg$actuals <- as.numeric(syn$actuals)
#
# mylogit <- glm(actuals ~ .,  data = syn_lg, family = "binomial",maxit=50)
# summary(mylogit)
# coef(summary(mylogit))[,'Pr(>|z|)']
mylogit <- cv.glmnet(y=as.numeric(syn$actuals), x=as.matrix(syn_lg), family = "binomial", alpha=1, type.measure = "class")

df = as.data.frame(as.matrix(coef(mylogit, s="lambda.min")))
colnames(df) = "score"
df %>% filter(score != 0)  %>% arrange(desc(score))

# res <- tibble::tibble(Description = split(long_hp[,1], long_hp[,3]))
# ggplot(res, aes_(x = ~Description)) + geom_bar()+
#   #theme_dose(font.size = 12) +
#   xlab(NULL) + ylab(NULL) +
#   ggupset::scale_x_upset(order_by = "freq")

