source("lib.R")

liftover_prm_07 <- read.table("hglft_genome_07.bed",
                                         col.names = c("chrom_name","prm_start","prm_end","name","score","strand"))
liftover_prm_10 <- read.table("hglft_genome_10.bed",
                                         col.names = c("chrom_name","prm_start","prm_end","name","score","strand"))

#файлы ortologs - попарные ортологи человека и мыши
ortologs_07 <- read.csv("gene_ortologs_07.csv", row.names = 1)
ortologs_10 <- read.csv("gene_ortologs_10.csv", row.names = 1)


TSS_intersect_07 <- read.table("hglft_genome_07.TSS.bed")
TSS_intersect_07 <- subset(TSS_intersect_07,select = c("V1","V8","V6","V4","V10"))
colnames(TSS_intersect_07) <- c("chromosome","TSS","strand","original_human_name","mouse_transcript_name")

temp <- strsplit(TSS_intersect_07$original_human_name,"@")
TSS_intersect_07$human_genes <- lapply(temp, "[[", 2)
temp <- strsplit(TSS_intersect_07$mouse_transcript_name,":")
TSS_intersect_07$mouse_genes <- lapply(temp, "[[", 2)
TSS_intersect_07$gene_pair <- paste(TSS_intersect_07$human_genes, TSS_intersect_07$mouse_genes, sep = "_")

TSS_intersect_07 <- TSS_intersect_07[TSS_intersect_07$gene_pair %in% ortologs_07$pair, ]
prm_coord_07 <- liftover_prm_07[c("prm_start","prm_end","name")]
rownames(prm_coord_07) <- prm_coord_07$name

temp <- TSS_intersect_07 %>%
  mutate(
    prm_coord_07[original_human_name,]
  )
TSS_intersect_07 <- subset(temp,select=-name)
TSS_intersect_07 <- TSS_intersect_07 %>%
  mutate(
    TSS_in=if_else(TSS>prm_start & TSS<prm_end,"yes","no")
  )
TSS_intersect_07 <- TSS_intersect_07 %>%
  mutate(
    TSS_dist=if_else(strand=="+",TSS-prm_end,prm_start-TSS)
  )

TSS_intersect_07 <- TSS_intersect_07[!(TSS_intersect_07$TSS_in=="no" & TSS_intersect_07$TSS_dist<0),]

temp <- TSS_intersect_07 %>%
  group_by(original_human_name) %>%
  mutate(
    best_dist=min(abs(TSS_dist))
  )

TSS_final_07 <- temp %>%
  group_by(original_human_name) %>%
  filter(abs(TSS_dist)==best_dist) %>%
  distinct(original_human_name, .keep_all = TRUE) %>%
  group_by(mouse_transcript_name) %>%
  distinct(mouse_transcript_name, .keep_all = TRUE)

TSS_final_07$name <- with(TSS_final_07, paste(original_human_name,mouse_transcript_name,sep="_"))


TSS_final_07 <- subset(TSS_final_07, select = c("chromosome","TSS","TSS","name","strand"))
TSS_final_07 <- add_column(TSS_final_07,score=0, .before = "strand")

write.table(TSS_final_07,"hglft_genome_07.only_TSS.bed", sep = "\t", row.names = FALSE,col.names = FALSE,quote = FALSE)



TSS_intersect_10 <- read.table("hglft_genome_10.TSS.bed")
TSS_intersect_10 <- subset(TSS_intersect_10,select = c("V1","V8","V6","V4","V10"))
colnames(TSS_intersect_10) <- c("chromosome","TSS","strand","original_human_name","mouse_transcript_name")

temp <- strsplit(TSS_intersect_10$original_human_name,"@")
TSS_intersect_10$human_genes <- lapply(temp, "[[", 2)
temp <- strsplit(TSS_intersect_10$mouse_transcript_name,":")
TSS_intersect_10$mouse_genes <- lapply(temp, "[[", 2)
TSS_intersect_10$gene_pair <- paste(TSS_intersect_10$human_genes, TSS_intersect_10$mouse_genes, sep = "_")

TSS_intersect_10 <- TSS_intersect_10[TSS_intersect_10$gene_pair %in% ortologs_10$pair, ]
prm_coord_10 <- liftover_prm_10[c("prm_start","prm_end","name")]
rownames(prm_coord_10) <- prm_coord_10$name

temp <- TSS_intersect_10 %>%
  mutate(
    prm_coord_10[original_human_name,]
  )
TSS_intersect_10 <- subset(temp,select=-name)

TSS_intersect_10 <- TSS_intersect_10 %>%
  mutate(
    TSS_in=if_else(TSS>prm_start & TSS<prm_end,"yes","no")
  )

TSS_intersect_10 <- TSS_intersect_10 %>%
  mutate(
    TSS_dist=if_else(strand=="+",TSS-prm_end,prm_start-TSS)
  )

TSS_intersect_10 <- TSS_intersect_10[!(TSS_intersect_10$TSS_in=="no" & TSS_intersect_10$TSS_dist<0),]

temp <- TSS_intersect_10 %>%
  group_by(original_human_name) %>%
  mutate(
    best_dist=min(abs(TSS_dist))
  )

TSS_final_10 <- temp %>%
  group_by(original_human_name) %>%
  filter(abs(TSS_dist)==best_dist) %>%
  distinct(original_human_name, .keep_all = TRUE) %>%
  group_by(mouse_transcript_name) %>%
  distinct(mouse_transcript_name, .keep_all = TRUE)

TSS_final_10$name <- with(TSS_final_10, paste(original_human_name,mouse_transcript_name,sep="_"))
TSS_final_10 <- subset(TSS_final_10, select = c("chromosome","TSS","TSS","name","strand"))
TSS_final_10 <- add_column(TSS_final_10,score=0, .before = "strand")

write.table(TSS_final_10,"hglft_genome_10.only_TSS.bed", sep = "\t", row.names = FALSE,col.names = FALSE,quote = FALSE)

rm(list=ls())




