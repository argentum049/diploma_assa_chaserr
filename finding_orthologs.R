source("lib.R")
source("mart_files.R")



sum_07 <- data_frame(read.csv("ASO_G0272888_AD_07.DE_Summary",sep = "\t",))

raw_genes <- strsplit(sum_07$geneSymbol,";")
raw_genes <- lapply(raw_genes, "[[", 1)
genes <- substring(raw_genes,6)


df07 <- data.frame(str_split_fixed(sum_07$prmtrID,"_",4))
colnames(df07) <- c("chrom","prmStart","prmEnd","strand")
df07$score <- "0"
df07$name <- paste(sum_07$prmtrID, genes, sep = "@")
df07 <- df07[,c("chrom","prmStart","prmEnd", "name", "score","strand")]
                   
write.table(df07,"prm07.bed", sep = "\t", row.names = FALSE,col.names = FALSE,quote = FALSE)
                   

result <- getLDS(attributes = c("hgnc_symbol"),
       filters = "hgnc_symbol", values = genes, mart = human,
       attributesL = c("mgi_symbol"), martL = mouse)

temp <- result$HGNC.symbol
unique_orthologs <- temp[ave(temp,temp,FUN=length)==1]

mouse_orthologs_to_check <- result[result$HGNC.symbol %in% unique_orthologs, "MGI.symbol"]


result2 <- getLDS(attributes = c("mgi_symbol"),
                 filters = "mgi_symbol", values = mouse_orthologs_to_check, mart = mouse,
                 attributesL = c("hgnc_symbol"), martL = human)

temp2 <- result2$MGI.symbol
unique_mouse_orthologs <- temp2[ave(temp2,temp2,FUN=length)==1]

ortologs_07 <- result2[result2$MGI.symbol %in% unique_mouse_orthologs, ]
ortologs_07 <- arrange(data_frame(ortologs_07), HGNC.symbol)
ortologs_07$pair <- paste(ortologs_07$HGNC.symbol, ortologs_07$MGI.symbol, sep = "_")


write.csv(ortologs_07,"gene_ortologs_07.csv")




sum_10 <- data_frame(read.csv("ASO_G0272888_AD_10.DE_Summary",sep = "\t"))

raw_genes <- strsplit(sum_10$geneSymbol,";")
raw_genes <- lapply(raw_genes, "[[", 1)
genes <- substring(raw_genes,6)


df10 <- data.frame(str_split_fixed(sum_10$prmtrID,"_",4))
colnames(df10) <- c("chrom","prmStart","prmEnd","strand")
df10$score <- "0"
df10$name <- paste(sum_10$prmtrID, genes, sep = "@")
df10 <- df10[,c("chrom","prmStart","prmEnd", "name", "score","strand")]


write.table(df10,"prm10.bed", sep = "\t", row.names = FALSE,col.names = FALSE,quote = FALSE)



result <- getLDS(attributes = c("hgnc_symbol"),
                 filters = "hgnc_symbol", values = genes, mart = human,
                 attributesL = c("mgi_symbol"), martL = mouse)

temp <- result$HGNC.symbol
unique_orthologs <- temp[ave(temp,temp,FUN=length)==1]

mouse_orthologs_to_check <- result[result$HGNC.symbol %in% unique_orthologs, "MGI.symbol"]


result2 <- getLDS(attributes = c("mgi_symbol"),
                  filters = "mgi_symbol", values = mouse_orthologs_to_check, mart = mouse,
                  attributesL = c("hgnc_symbol"), martL = human)

temp2 <- result2$MGI.symbol
unique_mouse_orthologs <- temp2[ave(temp2,temp2,FUN=length)==1]

ortologs_10 <- result2[result2$MGI.symbol %in% unique_mouse_orthologs, ]
ortologs_10 <- arrange(data_frame(ortologs_10), HGNC.symbol)
ortologs_10$pair <- paste(ortologs_10$HGNC.symbol, ortologs_10$MGI.symbol, sep = "_")

write.csv(ortologs_10,"gene_ortologs_10.csv")


rm(list=ls())

