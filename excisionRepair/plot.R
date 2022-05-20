
#length distribution
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)

#read all bed file
xlist <- list.files(pattern = "*.txt")
for(i in xlist) {
  x <- read.table((i), header = T)
  assign(i,x)
}
y <- tibble(xlist)

file_lent_dis.txt <- mutate(file_lent_dis.txt, percent=(Reads/sum(file_lent_dis.txt$Reads))*100, cellline= "Hela", time = "3h", repair = "Edu")
file_lent_dis.txt_p <- ggplot(data=file_lent_dis.txt, aes(Length, percent)) + geom_bar(stat='identity', position = 'identity', color='blue', width=0.5) + theme_tufte() + geom_rangeframe() + 
  theme(axis.text = element_text(size = 18), axis.title.x = (element_text(size = 20)), axis.title.y = (element_text(size = 20))) + ylab("Read Count (% of Total)") +
  ggtitle("3h") + theme(plot.title = element_text(size = 20)) +xlab("Sequence length (nt)") +
  theme(strip.text = element_text(size = 12)) +
  geom_vline(xintercept = c(26), linetype="dashed", color = "black", size=0.5) +
  scale_x_continuous(breaks=c(10, 20, 26, 30, 40, 50), labels= c("10","20", "26", "30", "40", "50"))
file_lent_dis.txt_p
ggsave('file_length.pdf', file_lent_dis.txt_p)

#nucleotide distribution
#!/usr/bin/env Rscript

library(ggplot2)

args<-commandArgs(TRUE)
input = args[1]
output = paste(input,'.pdf',sep='')
pdf(output)

df <- read.table(input, header = FALSE)
colnames(df) <- c("Length", "Position", "Base", "Percentage")
df$Base <- factor(df$Base, levels = c("G", "A", "C", "T"), ordered = TRUE)
ggplot(df, aes(x=Position, y=Percentage)) + 
  geom_bar(stat = "identity", aes(fill = Base)) + 
  scale_fill_manual(values = c("G" = "purple4", "C" = "dodgerblue4", "A" = "gree
n4", "T" = "orange")) + 
  facet_grid(Length ~ .) + 
  theme(plot.title = element_text(hjust = 0.5)) + ylab("Nucleotide frequency (of
 Total)") +
  theme(axis.text = element_text(size = 20), axis.title.x = (element_text(size =
                                                                            20)), axis.title.y = (element_text(size = 20))) +
  theme(legend.title=element_text(size=18), legend.text=element_text(size=18), s
        trip.text = element_text(size=20))
dev.off()

#genome-wide ts nts plot
file_TS.txti <- select(file_TS.txt, symbol, i, RPKM)
file_NTS.txti <- select(file_NTS.txt, symbol, i, RPKM)
file_TS.txtii <- summarise(group_by(file_TS.txti, i), mean(RPKM))
file_NTS.txtii <- summarise(group_by(file_NTS.txti, i), mean(RPKM))
file_TS.txtii <- mutate(file_TS.txtii, strand = "TS", cellline= "Hela", time = "48h")
file_NTS.txtii <- mutate(file_NTS.txtii, strand = "NTS", cellline= "Hela", time = "48h")
file <- rbind(file_TS.txtii, file_NTS.txtii)
colnames(file) <- c("i", "mean_RPKM", "strand", "cellline", "time")

file$strand <- factor(file$strand, levels = c("TS", "NTS"))
file$time <- factor(file$time, levels = c("6h", "9h", "12h", "24h", "48h"))
file_map <- ggplot(file, aes(i, mean_RPKM, color = strand)) + geom_line(size = 0.5) + ylab("Average Repair (RPKM)") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.key = element_rect(fill = "white", colour = "NA")) +
  scale_y_continuous(limits = c(0, 0.8)) + theme(axis.text = element_text(size = 8), axis.title.x = (element_text(size = 12)), axis.title.y = (element_text(size = 12)), axis.text.x = element_text(angle = 90, vjust=0.5)) +
  theme(plot.title = element_text(hjust = 0.5, size = 12)) + coord_fixed(ratio = 200) +
  scale_x_continuous(breaks=c(0, 25, 75, 125, 150), labels= c("-2kb", "TSS", "50%", "TES", "2kb")) + labs(x = NULL) +
  geom_vline(xintercept = c(25, 125), linetype="dashed", color = "black", size=0.5) + facet_grid(.~time) +
  theme(legend.position = c("bottom"), strip.text = element_text(size=12))

file_map
ggsave('file.pdf', file_map)

#chromatin states
library("ggplot2")
library("reshape2")
library(ggthemes)

#read all files
xlist <- list.files(pattern = "*.txt")
for(i in xlist) {
  x <- read.table((i))
  assign(i,x)
}
#data frame all file
library(dplyr)
y <- tibble(xlist)
#name all data
for(i in xlist) {
  i.tmp <- get(i)
  names(i.tmp) <- c("CHROM", "START", "END", "state", "a", "b", "c", "d", "e", "count") 
  assign(i, i.tmp)
}
rm(i.tmp, x, y)
rm(i, xlist)

file <- select(file_HMM.txt, CHROM, START, END, state, count)
file <- mutate(file, RPKM = (count*1000000000)/(1220244*(END-START)), column = c(1:641016), time = "3h", cell = "hela", rep = "R2")
file_a <- summarise(group_by(file, state), mean(RPKM))
file_a <- mutate(file_a, time="3h")
file$state <- factor(file$state, levels = c("1_Active_Promoter", "2_Weak_Promoter", "3_Poised_Promoter", "4_Strong_Enhancer", "5_Strong_Enhancer", "6_Weak_Enhancer", "7_Weak_Enhancer", "8_Insulator", "9_Txn_Transition", "10_Txn_Elongation", "11_Weak_Txn", "12_Repressed", "13_Heterochrom/lo", "14_Repetitive/CNV", "15_Repetitive/CNV"), ordered = TRUE)
file$time <- factor(file$time, levels=c("6h", "9h", "12h", "24h", "48h"))
file_p <- ggplot(data=file, aes(state, `mean(RPKM)`, fill = state)) + geom_bar(stat='identity', position = 'identity', color="black", width=0.75) +
  theme(axis.text = element_text(size = 10), axis.title.x = (element_text(size = 10)), axis.text.x =element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  ylab("RPKM") + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), legend.key = element_rect(fill = "white", colour = "NA")) +
  coord_fixed(ratio=4) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_fill_manual(values=c("red", "orange", "purple", "gold3","gold3", "yellow", "yellow","steelblue1",
                             "green4", "green4", "palegreen", "grey50", "grey", "grey", "grey")) +
  theme(axis.text.y = element_text(size = 12)) + facet_grid(time~.) +
  theme(plot.title = element_text(hjust = 0.5, size = 10)) +ylim(0,1.5) + theme(legend.title=element_text(size=8), legend.text=element_text(size=8), strip.text = element_text(size=12)) 

file_p
ggsave('file_p.pdf', file_p)

#DNase plot
library(dplyr)
library(wesanderson)
library(ggplot2)
#rm(list = ls())
#read all bed file
xlist <- list.files(pattern = "*.txt")
for(i in xlist) {
  x <- read.table((i))
  assign(i,x)
}
#data frame all file

y <- tibble(xlist)
#name all data
for(i in xlist) {
  i.tmp <- get(i)
  names(i.tmp) <- c("i", "percent") 
  assign(i, i.tmp)
}

rm(i.tmp, x, y)
rm(i, xlist)
file_DNase_Genea.txt <- mutate(file_DNase_Genea.txt, RPKM=(percent/1224133)*1000000000)
file_DNase_Gene.txti <- summarise(group_by(file_DNase_Genea.txt, i), mean(RPKM))
file_DNase_Gene.txti <- mutate(file_DNase_Gene.txti, cellline= "Hela", time = "3h", region = "Genes")
file_DNase_intergenica.txt <- mutate(file_DNase_intergenica.txt, RPKM=(percent/1224133)*1000000000)
file_DNase_intergenic.txti <- summarise(group_by(file_DNase_intergenica.txt, i), mean(RPKM))
file_DNase_intergenic.txti <- mutate(file_DNase_intergenic.txti, cellline= "Hela", time = "3h", region = "intergenic")
file$time <- factor(file$time, levels = c("6h","9h", "12h", "24h", "48h"))
file_map <- ggplot(file, aes(i, `mean(RPKM)`, color= region)) + geom_line(size = 1) + ylab("Average Repair (RPKM)") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.key = element_rect(fill = "white", colour = "NA")) +
  scale_y_continuous(limits = c(0, 2)) + theme(axis.text = element_text(size = 12), axis.title.x = (element_text(size = 12)), axis.title.y = (element_text(size = 12)), axis.text.x = element_text(angle = 90, vjust=0.5)) +
  theme(plot.title = element_text(hjust = 0.5, size = 15)) + coord_fixed(ratio = 20) +
  scale_x_continuous(breaks=c(0, 25, 50), labels= c("-1.5kb", "Peak Center", "1.5kb")) + labs(x = NULL) +
  scale_color_manual(values=c("#3182bd", "#9ecae1")) + facet_grid(.~time) +
  theme(legend.position="bottom", legend.text = element_text(size=12), legend.title = element_text(size=12), strip.text = element_text(size=12))  
file_map
ggsave('file.pdf', file_map)




