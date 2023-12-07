library(phyloseq)
library(vegan)
library(ggplot2)
install.packages("FSA")
library(FSA)

load("D:/Radka/NGS/Sekvenace Násilí_ATB + Narko/Sekvence Násilí+ATB/Orez_k_analyze_mysi/merged_J.R")
sort(sample_sums(MERGED))

PHYLOSEQ.final<-MERGED
sample_data(PHYLOSEQ.final)
# rarefaction - znormalizuje data, tak že vezme náhondě např.2000 sekvencí (podle počtu nejnižšího vzorku) a otestuje přítomnost OTU
PHYLOSEQ.final.rare<-rarefy_even_depth(PHYLOSEQ.final)
# proporce - taky normalizace, tím, že se vytvoří proporce asi sekvencí
PHYLOSEQ.final.prop<-transform_sample_counts(PHYLOSEQ.final,function(x) x/sum(x))
# spočítání alpha diverzity
RICH<-estimate_richness(PHYLOSEQ.final.rare)
write.table(RICH,file="D:/Radka/NGS/Sekvenace Násilí_ATB + Narko/Sekvence Násilí+ATB/Orez_k_analyze_mysi/Alpha_diversity.txt", sep="\t",quote = F)
# přidání alpha diverzity k sample datům
RICH<-data.frame(RICH,sample_data(PHYLOSEQ.final.rare))
ggplot(RICH,aes(x=Group,y=Shannon))+geom_boxplot(outlier.shape = NA)+geom_jitter()
ggplot(RICH,aes(x=Group,y=Observed))+geom_boxplot(outlier.shape = NA)+geom_jitter()

# alpha diversita
# můžeme se podívat histogramem na rozložení dat
hist(RICH$Observed)
{qqnorm(resid(model))
  qqline(resid(model))}

# Shannon
model<-lm(Observed~Group,data=RICH)
{qqnorm(resid(model))
  qqline(resid(model))}
summary(model)
anova(model)

model_log<-(lm(log10(Shannon)~Group,data=RICH))
{qqnorm(resid(model_log))
  qqline(resid(model_log))}

summary(model_log)
anova(model_log)

kruskal.test(Observed~Group,data=RICH)
dunnTest(Observed~Group,data=RICH)

# beta diverzita
BC<-vegdist(otu_table(PHYLOSEQ.final.prop))
JA<-vegdist(data.frame(otu_table(PHYLOSEQ.final.rare)),method = "jaccard",binary = T)

ord<-ordinate(PHYLOSEQ.final.prop,method = "PCoA", BC)
plot_ordination(PHYLOSEQ.final.prop,ord,color="Group")
plot_ordination(PHYLOSEQ.final.prop,ord,color="Group")+facet_wrap(.~Group)

ord<-ordinate(PHYLOSEQ.final.prop,method = "PCoA", JA)
plot_ordination(PHYLOSEQ.final.prop,ord,color="Group")
plot_ordination(PHYLOSEQ.final.prop,ord,color="Group")+facet_wrap(.~Group)

SD<-data.frame(sample_data(PHYLOSEQ.final.prop))
adonis2(BC~Group,data=SD)
adonis2(JA~Group,data=SD)

# pairwiseAdonis
library(pairwiseAdonis)

PAIRWISE.TAB<-function(x){
  LIST<-list()
  for(i in 2:length(x)){
    act<-x[[i]]
    LIST[[i-1]]<-act[1,]
  }
  DF<-do.call("rbind",LIST)
  DF<-data.frame(Contrast=names(x)[2:length(x)],
                 DF)
  DF
}

pp<-pairwise.adonis2(BC~Group,data=SD)
pp<-PAIRWISE.TAB(pp)
pp2<-pp[,c(1,5,6,4)]
names(pp2)<-c("Contrast","F","p","R2")  
pp2

pp<-pairwise.adonis2(JA~Group,data=SD)
pp<-PAIRWISE.TAB(pp)
pp2<-pp[,c(1,5,6,4)]
names(pp2)<-c("Contrast","F","p","R2")  
pp2

# korelace BC a JA
library(corrplot)
library(Hmisc)

plot(BC,JA)
cor.test(BC,JA)

```