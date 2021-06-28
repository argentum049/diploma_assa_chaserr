source("lib.R")


liftover_prm_07 <- data_frame(read.table("hglft_genome_07.bed",
                                      col.names = c("chrom_name","prm_start","prm_end","name","score","strand")))


biomart_file_07 <- liftover_prm_07
biomart_file_07$chrom <- substring(biomart_file_07$chrom_name ,4)

biomart_file_07$new_strand <- ifelse(biomart_file_07$strand=="+",1,-1)
biomart_file_07 <- subset(biomart_file_07,select = c("chrom","prm_start","prm_end","new_strand"))

biomart_file_07$chrom_region <- str_c(biomart_file_07$chrom,biomart_file_07$prm_start,
                                      biomart_file_07$prm_end,biomart_file_07$new_strand,sep = ":")

write.table(biomart_file_07$chrom_region,"chrom_coord_07.txt", row.names = FALSE,col.names = FALSE,quote = FALSE)




liftover_prm_10 <- data_frame(read.table("hglft_genome_10.bed",
                                      col.names = c("chrom_name","prm_start","prm_end","name","score","strand")))

biomart_file_10 <- liftover_prm_10
biomart_file_10$chrom <- substring(biomart_file_10$chrom_name ,4)

biomart_file_10$new_strand <- ifelse(biomart_file_10$strand=="+",1,-1)
biomart_file_10 <- subset(biomart_file_10,select = c("chrom","prm_start","prm_end","new_strand"))

biomart_file_10$chrom_region <- str_c(biomart_file_10$chrom,biomart_file_10$prm_start,
                                      biomart_file_10$prm_end,biomart_file_10$new_strand,sep = ":")

write.table(biomart_file_10$chrom_region,"chrom_coord_10.txt", row.names = FALSE,col.names = FALSE,quote = FALSE)


rm(list=ls())






