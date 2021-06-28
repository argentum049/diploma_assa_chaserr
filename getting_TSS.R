source("lib.R")
source("mart_files.R")

#файлы ortologs - попарные ортологи человека и мыши
ortologs_07 <- read.csv("gene_ortologs_07.csv", row.names = 1)
ortologs_10 <- read.csv("gene_ortologs_10.csv", row.names = 1)

#файлы chrom_coord - столбцы $chrom_region data_frame'ов biomart_file
chrom_coord_07 <- read.table("chrom_coord_07.txt") 
chrom_coord_10 <- read.table("chrom_coord_10.txt")


TSS_07 <- getBM(attributes = c("chromosome_name","transcription_start_site","transcript_start",
                            "transcript_end","ensembl_transcript_id_version","mgi_symbol", "strand") ,  
             filters = c("chromosomal_region","mgi_symbol"), 
             values = list(chrom_coord_07$V1, ortologs_07$MGI.symbol),
             mart=mouse)



TSS_07$chromosome_name <- paste0("chr",TSS_07$chromosome_name)

TSS_07$start <- TSS_07$transcription_start_site
TSS_07$name <- paste(TSS_07$ensembl_transcript_id_version,TSS_07$mgi_symbol,sep=":")
TSS_07$strand_plus_minus <- ifelse(TSS_07$strand==1,"+","-")
TSS_07$score <- 0

TSS_07 <- TSS_07[c("chromosome_name","transcription_start_site","transcription_start_site",
               "name","score","strand_plus_minus")]


write.table(TSS_07,"TSS_07.bed", sep = "\t", row.names = FALSE,col.names = FALSE,quote = FALSE)







TSS_10 <- getBM(attributes = c("chromosome_name","transcription_start_site", "transcript_start",
                            "transcript_end","ensembl_transcript_id_version","mgi_symbol","strand") ,  
             filters = c("chromosomal_region","mgi_symbol"), 
             values = list(chrom_coord_10$V1, ortologs_10$MGI.symbol),
             mart=mouse)


TSS_10$chromosome_name <- paste0("chr",TSS_10$chromosome_name)

TSS_10$start <- TSS_10$transcription_start_site
TSS_10$name <- paste(TSS_10$ensembl_transcript_id_version,TSS_10$mgi_symbol,sep=":")
TSS_10$strand_plus_minus <- ifelse(TSS_10$strand==1,"+","-")
TSS_10$score <- 0

TSS_10 <- TSS_10[c("chromosome_name","transcription_start_site","transcription_start_site",
                   "name","score","strand_plus_minus")]

write.table(TSS_10,"TSS_10.bed", sep = "\t", row.names = FALSE,col.names = FALSE,quote = FALSE)


rm(list=ls())