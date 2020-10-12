library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)


#read all bed file
xlist <- list.files(pattern = "*.txt")
for(i in xlist) {
  x <- read.table((i))
  assign(i,x)
}
y <- tbl_df(xlist)
#name all data
for(i in xlist) {
  i.tmp <- get(i)
  names(i.tmp) <- c("base", "position", "number") 
  assign(i, i.tmp)
}

for(i in xlist) {
  i.tmp <- get(i)
  i.tmp <- i.tmp[ !i.tmp$base=="N", ]
  assign(i, i.tmp)
}

for(i in xlist) {
  i.tmp <- get(i)
  i.tmp <- mutate(i.tmp, Position=position + 1)
  assign(i, i.tmp)
}

for(i in xlist) {
  i.tmp <- get(i)
  i.tmp <- select(i.tmp, base, Position, number)
  assign(i, i.tmp)
}

rm(i.tmp, x, y)
rm(i, xlist)


#nuc
MsmegWT5minBed_12nuc <- summarise(group_by(MsmegWT5minBed_12nuc.txt, Position), sum(number))
MsmegWT5minBed_13nuc <- summarise(group_by(MsmegWT5minBed_13nuc.txt, Position), sum(number))
MsmegWT5minBed_12nuc.txt <- mutate(MsmegWT5minBed_12nuc.txt, length="12nt", percent = (number/438348)*100, rep="Rep1")
MsmegWT5minBed_13nuc.txt <- mutate(MsmegWT5minBed_13nuc.txt, length="13nt", percent = (number/389155)*100, rep="Rep1")
                           MsmegWT5minBed_29nuc.txt, MsmegWT5minBed_30nuc.txt)

MsmegWT5minBed_nuc3 <- rbind(MsmegWT5minBed_12nuc.txt, MsmegWT5minBed_13nuc.txt) 
MsmegWT5minBed_nuc3$base <- factor(MsmegWT5minBed_nuc3$base, levels = c("G", "A", "C", "T"), ordered = TRUE)
MsmegWT5minBed_nuc3$length <- factor(MsmegWT5minBed_nuc3$length, levels = c("12nt", "13nt"), ordered = TRUE)
MsmegWT5minBed_nuc3 <- arrange(MsmegWT5minBed_nuc3, base)
MsmegWT5minBed_nuc3_p <- ggplot(MsmegWT5minBed_nuc3, aes(Position, percent)) + geom_bar(stat='identity', aes(fill = base)) + 
  scale_fill_manual(values = c("G" = "purple4", "C" = "dodgerblue4", "A" = "green4", "T" = "orange")) +
  facet_grid(length ~ .) + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) + ylab("Nucleotide frequency (% of Total)") +
  theme(axis.text = element_text(size = 15), axis.title.x = (element_text(size = 20)), axis.title.y = (element_text(size = 20))) +
  theme(legend.title=element_text(size=13), legend.text=element_text(size=12), strip.text = element_text(size=18))
  #ggtitle("wild type CPD 5min")
MsmegWT5minBed_nuc3_p
ggsave('MsmegWT5minBed_nuc3_p.pdf', MsmegWT5minBed_nuc3_p)

MsmegUvrD10minBed_12nuc <- summarise(group_by(MsmegUvrD10minBed_12nuc.txt, Position), sum(number))
MsmegUvrD10minBed_13nuc <- summarise(group_by(MsmegUvrD10minBed_13nuc.txt, Position), sum(number))
MsmegUvrD10minBed_12nuc.txt <- mutate(MsmegUvrD10minBed_12nuc.txt, length="12nt", percent = (number/541756)*100, rep="Rep1")
MsmegUvrD10minBed_13nuc.txt <- mutate(MsmegUvrD10minBed_13nuc.txt, length="13nt", percent = (number/790517)*100, rep="Rep1")
MsmegUvrD10minBed_nuc3 <- rbind(MsmegUvrD10minBed_12nuc.txt, MsmegUvrD10minBed_13nuc.txt) 
MsmegUvrD10minBed_nuc3$base <- factor(MsmegUvrD10minBed_nuc3$base, levels = c("G", "A", "C", "T"), ordered = TRUE)
MsmegUvrD10minBed_nuc3$length <- factor(MsmegUvrD10minBed_nuc3$length, levels = c("12nt", "13nt"), ordered = TRUE)
MsmegUvrD10minBed_nuc3 <- arrange(MsmegUvrD10minBed_nuc3, base)
MsmegUvrD10minBed_nuc3_p <- ggplot(MsmegUvrD10minBed_nuc3, aes(Position, percent)) + geom_bar(stat='identity', aes(fill = base)) + 
  scale_fill_manual(values = c("G" = "purple4", "C" = "dodgerblue4", "A" = "green4", "T" = "orange")) +
  facet_grid(length ~ .) + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) + ylab("Nucleotide frequency (% of Total)") +
  theme(axis.text = element_text(size = 15), axis.title.x = (element_text(size = 20)), axis.title.y = (element_text(size = 20))) +
  theme(legend.title=element_text(size=13), legend.text=element_text(size=12), strip.text = element_text(size=18))
  #ggtitle("uvrD CPD 5min")
MsmegUvrD10minBed_nuc3_p
ggsave('MsmegUvrD10minBed_nuc3_p.pdf', MsmegUvrD10minBed_nuc3_p)

#length
#length
wt5min_len <- read.table("MsmegWT5minBed_lent_dis.txt", header = T)
uvrD10min_len <- read.table("MsmegUvrD10minBed_lent_dis.txt", header = T)

sum(wt5min_len$Reads)
wt5min_len <- mutate(wt5min_len, percent=(Reads/2309284)*100, rep="Rep1")
wt5min_len_p <- ggplot(data=wt5min_len, aes(Length, percent)) + geom_bar(stat='identity', position = 'identity', color='blue', width=0.5) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text = element_text(size = 20), axis.title.x = (element_text(size = 20, family = "Helvetica")), axis.title.y = (element_text(size = 20, family = "Helvetica"))) + ylab("Read count (% of Total)") +
  theme(plot.title = element_text(hjust = 0.5, size = 15)) +xlab("Read length (nt)") +
  scale_y_continuous(limits = c(0, 40))
wt5min_len_p
ggsave('wt5min_len_p.pdf', wt5min_len_p)

sum(uvrD10min_len$Reads)
uvrD10min_len <- mutate(uvrD10min_len, percent=(Reads/2060476)*100, rep="Rep1")
uvrD10min_len_p <- ggplot(data=uvrD10min_len, aes(Length, percent)) + geom_bar(stat='identity', position = 'identity', color='blue', width=0.5) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("Read count (% of Total)") + xlab("Read length (nt)") +
  theme(axis.text = element_text(size = 20), axis.title.x = (element_text(size = 20, family = "Helvetica")), axis.title.y = (element_text(size = 20, family = "Helvetica"))) + 
  theme(plot.title = element_text(hjust = 0.5, size = 15)) +
  scale_y_continuous(limits = c(0, 40))
uvrD10min_len_p
ggsave('uvrD10min_len_p.pdf', uvrD10min_len_p)



