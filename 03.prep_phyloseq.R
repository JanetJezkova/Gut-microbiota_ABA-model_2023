library(Rcpp)
library(dada2)
library(phyloseq)
library(ShortRead)

load("D:/Radka/NGS/Sekvenace Násilí_ATB + Narko/Sekvence Násilí+ATB/Orez_k_analyze_mysi/seqtab.R")
setwd("D:/Radka/NGS/Sekvenace Násilí_ATB + Narko/Sekvence Násilí+ATB/Orez_k_analyze_mysi/")
#########################################
#Chimeric sequences######################
#########################################

#extraxt ASVs fasta from abundance matrix
FASTA<-DNAStringSet(colnames(seqtab))
names(FASTA)<-colnames(seqtab)
writeFasta(FASTA,"haplo.fasta")

# #elimination of chimeric sequences by uchime (Terminal command)
# system("usearch8.0.1517_i86linux32 -uchime_ref haplo.fasta -db ~/DB/gold.fasta -nonchimeras haplo.uchime.fasta -strand plus")


seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

sum(seqtab.nochim)/sum(seqtab)

############################################
#TAXONOMY###################################
############################################
# databaze se da stahnout tady: https://benjjneb.github.io/dada2/training.html

FASTA<-readDNAStringSet("haplo.fasta")
taxa <- assignTaxonomy(seqtab.nochim, 
                       refFasta="D:/Radka/NGS/Taxonomy reference/silva_nr99_v138.1_train_set.fa.gz", 
                       multithread=8,minBoot = 80)

head(unname(taxa))
############################################
#Create phyloseq object#####################
############################################

#OTU TABLE
seqtab<-otu_table(seqtab.nochim,taxa_are_rows = F)

#HAPLO
HAPLO<-readDNAStringSet("haplo.fasta")

#TAXO
TAXO<-tax_table(taxa)

PHYLOSEQ<-merge_phyloseq(seqtab,TAXO,HAPLO)
PHYLOSEQ_dupl<-PHYLOSEQ
sample_names(PHYLOSEQ)
write.table(seqtab,file = "D:/Radka/NGS/Sekvenace Násilí_ATB + Narko/Sekvence Násilí+ATB/Orez_k_analyze_mysi/merge.txt",sep="\t",quote = F)

#Provizorni sample metadata
list.files()

sample_names(PHYLOSEQ)<-gsub("S[0-9]*_","S",sample_names(PHYLOSEQ))

# META<-read.delim("metadata_mysi_bez H2O11.txt")
# SN<-paste0("S",META$ID_individual)
# META1<-META[duplicated(META$ID_individual)==FALSE,]
# META2<-META[duplicated(META$ID_individual)==TRUE,]
# META1$ID_individual<-paste0(META1$ID_individual,"_F_trus1R")
# META2$ID_individual<-paste0(META2$ID_individual,"_F_trus2R")
# META12<-rbind(META1,META2)
# SN<-paste0("S",META12$ID_individual)
# 
# sample_names(PHYLOSEQ)%in%SN
# PHYLOSEQ_dupl<-merge_phyloseq(PHYLOSEQ_dupl,META)
# 
# PHYLOSEQ_dupl<-merge_phyloseq(PHYLOSEQ_dupl,META)
# save(PHYLOSEQ_dupl,file = "D:/Radka/NGS/Sekvenace Násilí_ATB + Narko/Sekvence Násilí+ATB/Orez_k_analyze_mysi/phyloseq_dupl.R")

# Meta data podle Janet
META_J<-read.delim("metadata_mysi_bez H2O11.txt")
END<-sapply(strsplit(META_J$Sample,"_"),function(x) x[4], simplify = T)
NEW_NAME<-paste0("SS",META_J$ID_individual,"_F_",END,paste0(""))

# Ukáže nám jestli jsou vzorky ve phyloseq obsaženy v metadatech, ukáže nám ty co nejsou (!)
NEW_NAME%in%sample_names(PHYLOSEQ_dupl)
META_J<-sample_data(META_J)
sample_names(META_J)
sample_names(META_J)<-NEW_NAME

PHYLOSEQ_dupl<-merge_phyloseq(PHYLOSEQ_dupl,META_J)
save(PHYLOSEQ_dupl,file = "D:/Radka/NGS/Sekvenace Násilí_ATB + Narko/Sekvence Násilí+ATB/Orez_k_analyze_mysi/phyloseq_dupl_J.R")

# load("/media/kreising/DATA/data/Radka_Janet/DADA_PHYLOSEQ/PHYLOSEQ_dupl.R")
# write.table