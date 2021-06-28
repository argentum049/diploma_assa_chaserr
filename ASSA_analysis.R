source("lib2.R")
library(seqinr)

a <- read.csv("output_07_201.txt", row.names = 1, sep = "\t")
assa_a <- a %>%
  filter(Pvalue<0.01)


assa_07 <- read.csv("output_07.txt", row.names = 1, sep = "\t")
assa_pvalue_07 <- assa_07 %>%
  filter(Pvalue<0.01)

assa_10 <- read.csv("output_10.txt", row.names = 1, sep = "\t")
assa_pvalue_10 <- assa_10 %>%
  filter(Pvalue<0.01)


rnaup_07 <- read.table("rnaup_sites_07.txt",row.names = 1)
colnames(rnaup_07) <- c("name1",	"name2",	"siteStart1",	"siteEnd1",	"siteStart2",	"siteEnd2",	"siteLen",	"siteCmpl",	"siteScore",	"regStart1",	"regLen1",	"regStart2",	"regLen2",	"regGC",	"deltaG")
rnaup_07 <- rnaup_07 %>%
  filter(name2 %in% assa_pvalue_07$name2)

rnaup_10 <- read.table("rnaup_sites_10.txt",row.names = 1)
colnames(rnaup_10) <- c("name1",	"name2",	"siteStart1",	"siteEnd1",	"siteStart2",	"siteEnd2",	"siteLen",	"siteCmpl",	"siteScore",	"regStart1",	"regLen1",	"regStart2",	"regLen2",	"regGC",	"deltaG")
rnaup_10 <- rnaup_10 %>%
  filter(name2 %in% assa_pvalue_10$name2)

assa_07_wo_first <- read.csv("output_07_wo_first.txt", row.names = 1, sep = "\t")
assa_pvalue_07_wo_first <- assa_07_wo_first %>%
  filter(Pvalue<0.01)

assa_10_wo_first <- read.csv("output_10_wo_first.txt", row.names = 1, sep = "\t")
assa_pvalue_10_wo_first <- assa_10_wo_first %>%
  filter(Pvalue<0.01)


rnaup_07_wo_first <- read.table("rnaup_sites_07_wo_first.txt",row.names = 1)
colnames(rnaup_07_wo_first) <- c("name1",	"name2",	"siteStart1",	"siteEnd1",	"siteStart2",	"siteEnd2",	"siteLen",	"siteCmpl",	"siteScore",	"regStart1",	"regLen1",	"regStart2",	"regLen2",	"regGC",	"deltaG")
rnaup_07_wo_first <- rnaup_07_wo_first %>%
  filter(name2 %in% assa_pvalue_07_wo_first$name2)

rnaup_10_wo_first <- read.table("rnaup_sites_10_wo_first.txt",row.names = 1)
colnames(rnaup_10_wo_first) <- c("name1",	"name2",	"siteStart1",	"siteEnd1",	"siteStart2",	"siteEnd2",	"siteLen",	"siteCmpl",	"siteScore",	"regStart1",	"regLen1",	"regStart2",	"regLen2",	"regGC",	"deltaG")
rnaup_10_wo_first <- rnaup_10_wo_first %>%
  filter(name2 %in% assa_pvalue_10_wo_first$name2)











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





