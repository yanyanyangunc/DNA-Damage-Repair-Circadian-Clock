#read all bed file
xlist <- list.files(pattern = "*.bed")
for(i in xlist) {
  x <- read.table((i))
  assign(i,x)
}
#data frame all file
library(dplyr)
y <- tbl_df(xlist)
#name all data
for(i in xlist) {
  i.tmp <- get(i)
  names(i.tmp) <- c("CHROM", "START", "END",  "SYMBOL", "ID", "strand", "reads" ) 
  assign(i, i.tmp)
}
rm(i.tmp, x, y)
rm(i, xlist)

#read all bed file
xlist <- list.files(pattern = "*.txt")
for(i in xlist) {
  x <- read.table((i))
  assign(i,x)
}
#data frame all file
library(dplyr)
y <- tbl_df(xlist)
#name all data
for(i in xlist) {
  i.tmp <- get(i)
  names(i.tmp) <- c("CHROM", "START", "END",  "SYMBOL", "ID", "strand", "reads" ) 
  assign(i, i.tmp)
}
rm(i.tmp, x, y)
rm(i, xlist)

WT5minnm1 <- mutate(MsmegWT5min_neg_smeg_XRseq_minus.bed, STRAND = "NTS", TIME = "5min", type = "WT", RPKM = (reads*1000000000)/(278237*(END-START)))
WT5minnp1 <- mutate(MsmegWT5min_neg_smeg_XRseq_plus.bed, STRAND = "TS", TIME = "5min", type = "WT", RPKM = (reads*1000000000)/(278237*(END-START)))
WT5minpm1 <- mutate(MsmegWT5min_pos_smeg_XRseq_minus.bed, STRAND = "TS", TIME = "5min", type = "WT", RPKM = (reads*1000000000)/(278237*(END-START)))
WT5minpp1 <- mutate(MsmegWT5min_pos_smeg_XRseq_plus.bed, STRAND = "NTS", TIME = "5min", type = "WT", RPKM = (reads*1000000000)/(278237*(END-START)))

UvrD10minnm1 <- mutate(MsmegUvrD10min_neg_smeg_XRseq_minus.bed, STRAND = "NTS", TIME = "10min", type = "UvrD", RPKM = (reads*1000000000)/(606989*(END-START)))
UvrD10minnp1 <- mutate(MsmegUvrD10min_neg_smeg_XRseq_plus.bed, STRAND = "TS", TIME = "10min", type = "UvrD", RPKM = (reads*1000000000)/(606989*(END-START)))
UvrD10minpm1 <- mutate(MsmegUvrD10min_pos_smeg_XRseq_minus.bed, STRAND = "TS", TIME = "10min", type = "UvrD", RPKM = (reads*1000000000)/(606989*(END-START)))
UvrD10minpp1 <- mutate(MsmegUvrD10min_pos_smeg_XRseq_plus.bed, STRAND = "NTS", TIME = "10min", type = "UvrD", RPKM = (reads*1000000000)/(606989*(END-START)))


UvrD10min_NTS <- mutate(MsmegUvrD10min_smeg_XRseq_NTS.txt, STRAND = "NTS", TIME = "10min", type = "UvrD", RPKM = (reads*1000000000)/(606989*(END-START)))
UvrD10min_TS <- mutate(MsmegUvrD10min_smeg_XRseq_TS.txt, STRAND = "TS", TIME = "10min", type = "UvrD", RPKM = (reads*1000000000)/(606989*(END-START)))

WT5min_NTS <- mutate(MsmegWT5min_smeg_XRseq_NTS.txt, STRAND = "NTS", TIME = "5min", type = "WT", RPKM = (reads*1000000000)/(278237*(END-START)))
WT5min_TS <- mutate(MsmegWT5min_smeg_XRseq_TS.txt, STRAND = "TS", TIME = "5min", type = "WT", RPKM = (reads*1000000000)/(278237*(END-START)))

WT5min_TS1 <- rbind(WT5minpm1, WT5minnp1)
WT5min_NTS1 <- rbind(WT5minnm1, WT5minpp1)

UvrD10min_TS1 <- rbind(UvrD10minpm1, UvrD10minnp1)
UvrD10min_NTS1 <- rbind(UvrD10minnm1, UvrD10minpp1)

WT5min <- mutate(WT5min_TS, TS=RPKM, NTS=WT5min_NTS$RPKM)
WT5min <- select(WT5min, CHROM, START, END,  SYMBOL, ID, strand, TS, NTS)
WT5min <- mutate(WT5min, ratio=TS/NTS)

UvrD10min <- mutate(UvrD10min_TS, TS=RPKM, NTS=UvrD10min_NTS$RPKM)
UvrD10min <- select(UvrD10min, CHROM, START, END,  SYMBOL, ID, strand, TS, NTS)
UvrD10min <- mutate(UvrD10min, ratio=TS/NTS)

WT5min_top25 <- inner_join(WT5min, smeg_top25, by = c("ID"="ID"))  
UvrD10min_top25 <- inner_join(UvrD10min, smeg_top25, by = c("ID"="ID"))
WT5min_top25a <- select(WT5min_top25, ratio)
WT5min_top25a <- WT5min_top25a[!is.infinite(rowSums(WT5min_top25a)),]
library(IDPmisc)
WT5min_top25a <- NaRV.omit(WT5min_top25a)
mean(WT5min_top25a)
mean(UvrD10min_top25$ratio)

WT5min_top25_p <- ggplot(WT5min_top25, aes(log2(ratio))) + geom_histogram(fill=I("deepskyblue3")) +
  geom_vline(xintercept = c(0), linetype="dashed", color = "black", size=0.5) +
  ylab("Number of genes") + xlab("log2 TS/NTS Repair") +
  scale_x_continuous(limits = c(-5, 5)) +
  theme(axis.text = element_text(size = 15), axis.title.x = (element_text(size = 15)), axis.title.y = (element_text(size = 15))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_fixed(ratio=0.025) + scale_y_continuous(limits = c(0, 350))

WT5min_top25_p
ggsave('WT5min_top25_TSNTSratio.pdf', WT5min_top25_p)

UvrD10min_top25_p <- ggplot(UvrD10min_top25, aes(log2(ratio))) + geom_histogram(fill=I("deepskyblue3")) +
  geom_vline(xintercept = c(0), linetype="dashed", color = "black", size=0.5) +
  ylab("Number of genes") + xlab("log2 TS/NTS Repair") +
  scale_x_continuous(limits = c(-5, 5)) +
  theme(axis.text = element_text(size = 15), axis.title.x = (element_text(size = 15)), axis.title.y = (element_text(size = 15))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_fixed(ratio=0.025) + scale_y_continuous(limits = c(0, 350))

UvrD10min_top25_p
ggsave('UvrD10min_top25_TSNTSratio.pdf', UvrD10min_top25_p)


#total
WT5min_p <- ggplot(WT5min, aes(log2(ratio))) + geom_histogram() +
  geom_vline(xintercept = c(0), linetype="dashed", color = "black", size=0.5) +
  ylab("Number of genes") + xlab("log2 TS/NTS Repair") +
  scale_x_continuous(limits = c(-5, 5)) +
  theme(axis.text = element_text(size = 15), axis.title.x = (element_text(size = 15)), axis.title.y = (element_text(size = 15))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_fixed(ratio=0.0075)

WT5min_p
ggsave('WT5min_TSNTSratio.pdf', WT5min_p)

UvrD10min_p <- ggplot(UvrD10min, aes(log2(ratio))) + geom_histogram() +
  geom_vline(xintercept = c(0), linetype="dashed", color = "black", size=0.5) +
  ylab("Numbers of Genes") + xlab("log2 TS/NTS Repair") +
  scale_x_continuous(limits = c(-5, 5)) +
  theme(axis.text = element_text(size = 15), axis.title.x = (element_text(size = 15)), axis.title.y = (element_text(size = 15))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_fixed(ratio=0.006)

UvrD10min_p
ggsave('UvrD10min_TSNTSratio.pdf', UvrD10min_p)
