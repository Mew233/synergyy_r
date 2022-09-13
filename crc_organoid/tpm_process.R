library(dplyr)
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)
library(rentrez)
library(pheatmap)

## for organoid
crc_tibble = as_tibble(fread("crc_organoid/CRC_list_of_samples_RNA_counts.csv"))
Gene_names <- crc_tibble$gene
Genes2 <- select(org.Hs.eg.db, keys = Gene_names,
                 columns = c('ENTREZID'), keytype = 'SYMBOL')
joined <- left_join(crc_tibble, Genes2, 
                           by = c("gene" = "SYMBOL"))
joined <- joined %>% filter(!is.na(ENTREZID))


#https://www.jianshu.com/p/a2351dacab51
length<- as_tibble(fread("crc_organoid/All_hg19gene_len.txt"))#自己在Excel中把网盘里的txt文件基因和长度共两列提取出来
joined_slen <- joined  %>%left_join(., length, 
                                    by = c("gene" = "Gene")) %>% filter(!is.na(Length))

mycounts<-joined_slen#最后两列Length是基因长度和entrezid
head(mycounts)

#TPM计算
kb <- mycounts$Length / 1000
kb
countdata <- mycounts[,2:69]

countdata <- mutate_all(countdata, function(x) as.numeric(as.character(x)))
rpk <- countdata / kb
rpk
tpm <- t(t(rpk)/colSums(rpk) * 1000000)

log2_tpm <- as.tibble(tpm)%>% mutate(across(where(is.numeric), ~log2(. + 1)))
rownames(log2_tpm) <- mycounts$ENTREZID
write.table(log2_tpm,file="crc_organoid/crc_exp.csv",sep=",",quote=F)

#==========================================================================
    #   optional-在xean下载counts,用python parse.
    #   import pandas as pd
    #   df = pd.read_csv('TCGA-DLBC.htseq_counts.tsv',sep='\t',index_col=0)
    #   df = df.iloc[:-5]
    #   df.transform(lambda x: (2**x)-1).to_csv("tgca_DLBC_RNA_counts.csv")

## for tgca
tgca_tibble = as_tibble(fread("crc_organoid/tgca_DLBC_RNA_counts.csv"))
Gene_names <- tgca_tibble$Ensembl_ID
ensLookup <- gsub("\\.[0-9]*$", "", Gene_names)
tgca_tibble$Ensembl_ID <- ensLookup
Genes2 <- select(org.Hs.eg.db, keys = ensLookup,
                 columns = c('ENTREZID','SYMBOL'), keytype = 'ENSEMBL')
joined <- left_join(tgca_tibble, Genes2, 
                    by = c("Ensembl_ID" = "ENSEMBL"))
joined <- joined %>% dplyr::filter(!is.na(ENTREZID))

length<- as_tibble(fread("crc_organoid/All_hg19gene_len.txt"))#自己在Excel中把网盘里的txt文件基因和长度共两列提取出来
joined_slen <- joined  %>%left_join(., length, 
                                    by = c("SYMBOL" = "Gene")) %>% dplyr::filter(!is.na(Length))

mycounts<- joined_slen%>% 
  # Base the removal on the "Age" column
  distinct(ENTREZID, .keep_all = TRUE) #最后三列entreid (去除duplicates), symbol 和length


#TPM计算
kb <- mycounts$Length / 1000
countdata <- mycounts[,2:49] #ncol(mycounts)-3

countdata <- mutate_all(countdata, function(x) as.numeric(as.character(x)))
rpk <- countdata / kb
tpm <- t(t(rpk)/colSums(rpk) * 1000000)

log2_tpm <- as.tibble(tpm)%>% mutate(across(where(is.numeric), ~log2(. + 1)))
rownames(log2_tpm) <- mycounts$ENTREZID
write.table(log2_tpm,file="crc_organoid/tcga_DLBC_exp.csv",sep=",",quote=F)


#look at expression profile on selected genes==========================================================================
# selected genes 2674
# input_df = input_df %>% filter(as.numeric(rownames(.)) %in% selected_genes$V1)
# row_mean <- as_tibble(rowMeans(input_df))
# row_mean$entreID <- rownames(input_df)
# 
# row_mean = row_mean  %>% arrange(desc(value))
# Genes2 <- select(org.Hs.eg.db, keys = row_mean$entreID,
#                  columns = c('SYMBOL'), keytype = 'ENTREZID')













