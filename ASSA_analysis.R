source("lib.R")

assa_07 <- read.csv("output_07.txt", row.names = 1, sep = "\t")
assa_07 <- assa_07 %>%
  filter(Pvalue<0.01)

rnaup_07 <- read.table("rnaup_sites_07.txt",row.names = 1)
colnames(rnaup_07) <- c("name1",	"name2",	"siteStart1",	"siteEnd1",	"siteStart2",	"siteEnd2",	"siteLen",	"siteCmpl",	"siteScore",	"regStart1",	"regLen1",	"regStart2",	"regLen2",	"regGC",	"deltaG")
rnaup_07 <- rnaup_07 %>%
  filter(name2 %in% assa_07$name2) %>%
  filter(deltaG<0)
rnaup_07$Pvalue <- assa_07$Pvalue[match(rnaup_07$name2,assa_07$name2)]
rnaup_07$Padj <- assa_07$Padj[match(rnaup_07$name2,assa_07$name2)]
write.table(rnaup_07,"rnaup_sites_07_only_values.txt",quote=F,sep = "\t")


assa_10 <- read.csv("output_10.txt", row.names = 1, sep = "\t")
assa_10 <- assa_10 %>%
  filter(Pvalue<0.01)

rnaup_10 <- read.table("rnaup_sites_10.txt",row.names = 1)
colnames(rnaup_10) <- c("name1",	"name2",	"siteStart1",	"siteEnd1",	"siteStart2",	"siteEnd2",	"siteLen",	"siteCmpl",	"siteScore",	"regStart1",	"regLen1",	"regStart2",	"regLen2",	"regGC",	"deltaG")
rnaup_10 <- rnaup_10 %>%
  filter(name2 %in% assa_10$name2) %>%
  filter(deltaG<0)
rnaup_10$Pvalue <- assa_10$Pvalue[match(rnaup_10$name2,assa_10$name2)]
rnaup_10$Padj <- assa_10$Padj[match(rnaup_10$name2,assa_10$name2)]
write.table(rnaup_10,"rnaup_sites_10_only_values.txt",quote=F,sep = "\t")

rnaup_total_mouse <- rbind(rnaup_07,rnaup_10)
rnaup_total_mouse <- rnaup_total_mouse %>%
  mutate(
    key=paste(name2, siteStart1)
  )
rnaup_total_mouse <- rnaup_total_mouse %>%
  distinct(key, .keep_all = TRUE)
write.table(rnaup_total_mouse,"rnaup_total_mouse.txt",quote=F,sep = "\t")






human_assa_07 <- read.csv("output_hg_07.txt",row.names = 1, sep = "\t")
human_assa_07 <- human_assa_07 %>%
  filter(Pvalue<0.01)

rnaup_hg_07 <- read.table("hg_sites_07.txt", row.names = 1)
colnames(rnaup_hg_07) <- c("name1",	"name2",	"siteStart1",	"siteEnd1",	"siteStart2",	"siteEnd2",	"siteLen",	"siteCmpl",	"siteScore",	"regStart1",	"regLen1",	"regStart2",	"regLen2",	"regGC",	"deltaG")
rnaup_hg_07 <- rnaup_hg_07 %>%
  filter(name2 %in% human_assa_07$name2) %>%
  filter(deltaG<0)
rnaup_hg_07$Pvalue <- human_assa_07$Pvalue[match(rnaup_hg_07$name2,human_assa_07$name2)]
rnaup_hg_07$Padj <- human_assa_07$Padj[match(rnaup_hg_07$name2,human_assa_07$name2)]
write.table(rnaup_hg_07,"hg_sites_07_only_values.txt",quote=F,sep = "\t")


human_assa_10 <- read.csv("output_hg_10.txt",row.names = 1, sep = "\t")
human_assa_10 <- human_assa_10 %>%
  filter(Pvalue<0.01)

rnaup_hg_10 <- read.table("hg_sites_10.txt", row.names = 1)
colnames(rnaup_hg_10) <- c("name1",	"name2",	"siteStart1",	"siteEnd1",	"siteStart2",	"siteEnd2",	"siteLen",	"siteCmpl",	"siteScore",	"regStart1",	"regLen1",	"regStart2",	"regLen2",	"regGC",	"deltaG")
rnaup_hg_10 <- rnaup_hg_10 %>%
  filter(name2 %in% human_assa_10$name2) %>%
  filter(deltaG<0)
rnaup_hg_10$Pvalue <- human_assa_10$Pvalue[match(rnaup_hg_10$name2,human_assa_10$name2)]
rnaup_hg_10$Padj <- human_assa_10$Padj[match(rnaup_hg_10$name2,human_assa_10$name2)]
write.table(rnaup_hg_10,"hg_sites_10_only_values.txt",quote=F,sep = "\t")

rnaup_total_human <- rbind(rnaup_hg_07,rnaup_hg_10)
rnaup_total_human <- rnaup_total_human %>%
  mutate(
    key=paste(name2, siteStart1)
  )
rnaup_total_human <- rnaup_total_human %>%
  distinct(key, .keep_all = TRUE)
write.table(rnaup_total_human,"rnaup_total_human.txt",quote=F,sep = "\t")








#assa_07_wo_first <- read.csv("output_07_wo_first.txt", row.names = 1, sep = "\t")
#assa_pvalue_07_wo_first <- assa_07_wo_first %>%
#  filter(Pvalue<0.01)
#
#assa_10_wo_first <- read.csv("output_10_wo_first.txt", row.names = 1, sep = "\t")
#assa_pvalue_10_wo_first <- assa_10_wo_first %>%
#  filter(Pvalue<0.01)
#
#
#rnaup_07_wo_first <- read.table("rnaup_sites_07_wo_first.txt",row.names = 1)
#colnames(rnaup_07_wo_first) <- c("name1",	"name2",	"siteStart1",	"siteEnd1",	"siteStart2",	"siteEnd2",	"siteLen",	"siteCmpl",	"siteScore",	"regStart1",	"regLen1",	"regStart2",	"regLen2",	"regGC",	"deltaG")
#rnaup_07_wo_first <- rnaup_07_wo_first %>%
#  filter(name2 %in% assa_pvalue_07_wo_first$name2)
#
#rnaup_10_wo_first <- read.table("rnaup_sites_10_wo_first.txt",row.names = 1)
#colnames(rnaup_10_wo_first) <- c("name1",	"name2",	"siteStart1",	"siteEnd1",	"siteStart2",	"siteEnd2",	"siteLen",	"siteCmpl",	"siteScore",	"regStart1",	"regLen1",	"regStart2",	"regLen2",	"regGC",	"deltaG")
#rnaup_10_wo_first <- rnaup_10_wo_first %>%
#  filter(name2 %in% assa_pvalue_10_wo_first$name2)





#Chaserr_exons <- read.fasta("Mus_musculus_Chaserr_202_sequence.fa",as.string = T)
#Chaserr <- read.fasta("Chaserr.fna",as.string = T)



#последовательность изоформы Chaserr - ENSMUST00000184554
#Chaserr <- "CGCGCGCGATGTGGGAGCTCGCGGCCGGAGCGCCCGGGGAGGCCGGGCCCACGACGCCCGTGGCTCCGCTGCGGCAGCGGCGGTGCTGGTGTCGCGCGCGGCCGGGAGGCGGCTTCGCGCCGCGGGCGGGAGGCTGCGGCGGGCGACCCGTCCTGGACACGCGAGGAAGAGCGAGCCGATGGCGGCAGGGGCCGCGCTTCGACCCGGTAACTTAGAAGATGATAATTAATGTGGTTGCTGATAATTCTGAATAAATACAGCTTTTATCCCAGGTGTGCCATTTTGAAGACTGAGACCATAGAGTTCTAAGAATAAAGGAAAGAGCCCTTGGGAAATTATTATATATAGCAAAAATGTGAATCCTCAGATGGAATGAAAGGCCTGCACCATAGACATCGAAGCATTTAAATTTTTTTCTTCTAATTTTTTATGAAGCACCCCGCTTGAAGAGTTTGAAATGGACTTTACCACTGAGAAATCAAGATGGCAGCCCATTATGGGGAATTGAGGAAAATGGATTAATGCAAGAATGCTGTAATATTATACAACCAACACAGGATTCTTTTAATGTGGATTCCATGAAATGAATGATTCTTACCCAACACAAATGGACAGTGGAATTTACTTCCTAAAGACTTGTTACATGTCATGTACATTTTTGACATCTGGAGAAGACTCTACAATTCTACAAATGGTAGTTTGTATTCCTGGAATTTCTTGCAGTTTGATCTGAAGTGACCTTATGGAATGTTAACTTTAATAAAATCTCTAAAACTTAAAAA"

#экзоны, их последовательности и границы
#first_exon <- "CGCGCGCGATGTGGGAGCTCGCGGCCGGAGCGCCCGGGGAGGCCGGGCCCACGACGCCCGTGGCTCCGCTGCGGCAGCGGCGGTGCTGGTGTCGCGCGCGGCCGGGAGGCGGCTTCGCGCCGCGGGCGGGAGGCTGCGGCGGGCGACCCGTCCTGGACACGCGAGGAAGAGCGAGCCGATGGCGGCAGGGGCCGCGCTTCGACCCG"
#first_len <- nchar(first_exon)
#second_exon <- "GTAACTTAGAAGATGATAATTAATGTGGTTGCTGATAATTCTGAATAAATACAGCTTTTATCCCAG"
#second_len <- nchar(second_exon)
#third_exon <- "GTGTGCCATTTTGAAGACTGAGACCATAGAGTTCTAAGAATAAAGGAAAGAGCCCTTGGGAAATTATTATATATAGCAA"
#third_len <- nchar(third_exon)
#fourth_exon <- "AAATGTGAATCCTCAGATGGAATGAAAGGCCTGCACCATAGACATCGAAGCATTTA"
#fourth_len <- nchar(fourth_exon)
#fifth_exon <- "AATTTTTTTCTTCTAATTTTTTATGAAGCACCCCGCTTGAAGAGTTTGAAATGGACTTTACCACTGAGAAATCAAGATGGCAGCCCATTATGGGGAATTGAGGAAAATGGATTAATGCAAGAATGCTGTAATATTATACAACCAACACAGGATTCTTTTAATGTGGATTCCATGAAATGAATGATTCTTACCCAACACAAATGGACAGTGGAATTTACTTCCTAAAGACTTGTTACATGTCATGTACATTTTTGACATCTGGAGAAGACTCTACAATTCTACAAATGGTAGTTTGTATTCCTGGAATTTCTTGCAGTTTGATCTGAAGTGACCTTATGGAATGTTAACTTTAATAAAATCTCTAAAACTTAAAAA"
#fifth_len <- nchar(fifth_exon)





