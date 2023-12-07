library(dada2)
library(ggplot2)
setwd("D:/Radka/NGS/Sekvenace Násilí_ATB + Narko/Sekvence Násilí+ATB/Orez_k_analyze_mysi/trimmednew_0523")

#List of forward and reverse reads
LIST<-list.files()
F_reads<-LIST[grep("_trus-trimmed-pair1.fastq.gz",LIST)]
R_reads<-LIST[grep("_trus-trimmed-pair2.fastq.gz",LIST)]
F_reads_TO<-paste0("S",F_reads)
R_reads_TO<-paste0("S",R_reads)

file.rename(from = F_reads,to = F_reads_TO)
file.rename(from = R_reads,to=R_reads_TO)

LIST<-list.files()

F_reads<-LIST[grep("_trus-trimmed-pair1.fastq.gz",LIST)]
R_reads<-LIST[grep("_trus-trimmed-pair2.fastq.gz",LIST)]

#graphical representation of quality profiles
QP.f<-plotQualityProfile(F_reads[1:10],aggregate = TRUE)+ggtitle("Forward reads")
QP.2<-plotQualityProfile(R_reads[1:10],aggregate = TRUE)+ggtitle("Rewerse reads")


# QP.f
# QP.2

# ggsave(QP.f,filename = "/media/kreising/DATA/data/Radka_Janet/Forward_reads.pdf")
# ggsave(QP.2,filename = "/media/kreising/DATA/data/Radka_Janet/Reverse_reads.pdf")

sample.names<-gsub("_trus-trimmed-pair1.fastq.gz","",F_reads)
sample.names<-gsub("-assigned-","",sample.names)
filtFs <- paste0(sample.names, "_READ1_filt.fastq.gz")
filtRs <- paste0(sample.names, "_READ2_filt.fastq.gz")

#Quality filtering
for(x in 1:length(F_reads)) {
  print(sample.names[x])
  fastqPairedFilter(c(F_reads[x], R_reads[x]), c(filtFs[x], filtRs[x]),
                    maxN=0, maxEE=2, minQ=2,truncQ=2,
                    compress=TRUE, verbose=TRUE, matchIDs = TRUE,
                    minLen = c(270,190),truncLen = c(270,190))
}


###############################################
#DADA DENOISING################################
###############################################

#These commands denoise quality-filtered fastq files and build abundance matrix,
#(samples in rows, ASVs in columns)

list.files("D:/Radka/NGS/Sekvenace Násilí_ATB + Narko/Sekvence Násilí+ATB/Orez_k_analyze_mysi/trimmednew_0523")

#List of quality filtered fastq files
fns <- list.files("D:/Radka/NGS/Sekvenace Násilí_ATB + Narko/Sekvence Násilí+ATB/Orez_k_analyze_mysi/trimmednew_0523")
fns<-fns[grep("filt",fns)]
setwd("D:/Radka/NGS/Sekvenace Násilí_ATB + Narko/Sekvence Násilí+ATB/Orez_k_analyze_mysi/trimmednew_0523")

#fns<-fns[1:283] TOHLE VYBIRA SUMSET VZORKU PRO DADU - ZRUSIT
fastqs <- fns[grepl("_filt.fastq", fns)]
# fastqs <- sort(fastqs) 


fnFs <- fastqs[grepl("READ1_filt.fastq.gz", fastqs)] 
fnRs <- fastqs[grepl("READ2_filt.fastq.gz", fastqs)] 
sample.names <- gsub("_READ1_filt.fastq.gz","",fnFs)

#fastq dereplication
derepFs <- derepFastq(fnFs,verbose=T)
derepRs <- derepFastq(fnRs, verbose=T)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

#deoising
dadaFs <- dada(derepFs, selfConsist = TRUE,MAX_CONSIST=20,multithread = TRUE)
dadaRs <- dada(derepRs, selfConsist = TRUE,MAX_CONSIST=20,multithread = TRUE)

#merge denoised forward and reverse ASVs
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, 
                      minOverlap = 10,maxMismatch=0,justConcatenate=F)

#abundance matrix
seqtab <- makeSequenceTable(mergers)

save(seqtab,file = "D:/Radka/NGS/Sekvenace Násilí_ATB + Narko/Sekvence Násilí+ATB/Orez_k_analyze_mysi/seqtab.R")
write.table(seqtab,file = "D:/Radka/NGS/Sekvenace Násilí_ATB + Narko/Sekvence Násilí+ATB/Orez_k_analyze_mysi/seqtab.R.txt",sep="\t",quote = F)
sum(seqtab)
rowSums(seqtab)
colSums(seqtab)
colnames(seqtab)
