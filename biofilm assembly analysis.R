setwd("~/Desktop/Dissertation/aim2/analysis/")

##--------------------------------------------------------------------
#load pyloseq
load("~/Desktop/Dissertation/aim2/bacteria/phyloseq.RData")
load("~/Desktop/Dissertation/aim2/fungi/phyloseq.RData")

##--------------------------------------------------------------------
# load libraries
library(ape)
library(phyloseq)
library(ggplot2)
library(vegan)
library(tibble)
#library(ggdendro)
library(plyr)
library(dplyr)
#library(gridExtra)
#library(pals)
library(reshape)
library(ggsignif)
library(pairwiseAdonis)
#library(indicspecies)
#library(dendextend)
library(cowplot)
library(scales)
library(RColorBrewer)
#library(randomcoloR)
library(ALDEx2)
library(pulsar)
library(picante)
library(ComplexHeatmap)
library(funrar)
##--------------------------------------------------------------------
# alpha diversity

#ps <- rarefy_even_depth(ps)
a = estimate_richness(ps, measures=c("Shannon","Observed"))
a$code = rownames(a)
a$code = gsub("X","",a$code)

samp = as.matrix(ps@otu_table)
tree = phy_tree(ps)
comm.pd = pd(samp,tree,include.root=F)
comm.pd$code = rownames(comm.pd)

df = inner_join(comm.pd,ps@sam_data,by="code")
df = left_join(a,df,by="code")
df$day <- factor(df$day, levels=c("0","whey","1","4","7"))

p1 = ggplot(df, aes(x=day, y=Observed,color=day)) + 
  geom_boxplot() + theme_bw() + 
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_color_manual(values = c("#E41A1C","#FF7F00","#377EB8","#4DAF4A","#984EA3")) +
  geom_signif(comparisons = list(
    c("1", "7"),
    c("1", "4"),
    c("4", "7")), 
    map_signif_level=TRUE,step_increase = 0.1,color = "black",
    size = 0.3,textsize = 2)
p1
p2 = ggplot(df, aes(x=day, y=Shannon,color=day)) + 
  geom_boxplot() + theme_bw() + 
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_color_manual(values = c("#E41A1C","#FF7F00","#377EB8","#4DAF4A","#984EA3")) +
  geom_signif(comparisons = list(
    c("1", "7"),
    c("1", "4"),
    c("4", "7")), 
    map_signif_level=TRUE,step_increase = 0.1,color = "black",
    size = 0.3,textsize = 2)
p2
p3 = ggplot(df, aes(x=day, y=PD,color=day)) + 
  geom_boxplot() + theme_bw() + 
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  scale_color_manual(values = c("#E41A1C","#FF7F00","#377EB8","#4DAF4A","#984EA3")) +
  geom_signif(comparisons = list(
    c("1", "7"),
    c("1", "4"),
    c("4", "7")), 
    map_signif_level=TRUE,step_increase = 0.1,color = "black",
    size = 0.3,textsize = 2)
p3

fig2 <- plot_grid(p1, p2, p3,
                  #labels = c("A","", "B",""),rel_widths = c(1,0.5),
                  ncol = 3, nrow = 1,label_size = 12)

ggsave(fig2, file="~/Desktop/fig2.eps",height = 3, width = 8,units = "in",dpi = 300)

##--------------------------------------------------------------------
# beta diversity clustering

#include day 0
prop.ps <- transform_sample_counts(ps, function(otu) otu/sum(otu))

ord.pcoa.bray <- ordinate(prop.ps, method="PCoA", distance="bray")
p = plot_ordination(prop.ps, ord.pcoa.bray, color="day")
p$data$day <- as.character(p$data$day)
p$data$day <- factor(p$data$day, levels=c("0","whey","1","4","7"))
p1 = p + theme_bw() + 
  theme(legend.title=element_blank(),
        legend.text = element_text(size=10)) +
  scale_color_manual(values = c("#E41A1C","#FF7F00","#377EB8","#4DAF4A","#984EA3")) + ggtitle("Bray-Curtis")
p1

ord.pcoa.bray <- ordinate(prop.ps, method="PCoA", distance="wunifrac")
p = plot_ordination(prop.ps, ord.pcoa.bray, color="day")
p$data$day <- as.character(p$data$day)
p$data$day <- factor(p$data$day, levels=c("0","whey","1","4","7"))
p2 = p + theme_bw() + 
  theme(legend.title=element_blank(),
        legend.text = element_text(size=10)) +
  scale_color_manual(values = c("#E41A1C","#FF7F00","#377EB8","#4DAF4A","#984EA3")) + ggtitle("Weighted UniFrac")
p2

fig3 <- plot_grid(p1, p2,
                  #labels = c("A", "B"),
                  ncol = 2, nrow = 1,label_size = 12)

ggsave(fig3,file="~/Desktop/fig3.eps",height = 3, width = 8,units = "in",dpi = 300)



##--------------------------------------------------------------------
# test if biofilm developed under different time same or not
ps.prop = transform_sample_counts(ps, function(otu) otu/sum(otu))
beta <- phyloseq::distance(ps.prop, method = "wunifrac")
# make a data frame from the sample_data
sampledf <- data.frame(sample_data(ps.prop))
# Adonis test
adonis2(beta ~ day, data = sampledf)
pairwise.adonis(otu_table(ps.prop),sampledf$day)

sampledf$day[sampledf$day != 0] <- "other"
adonis2(beta ~ day, data = sampledf)
pairwise.adonis(otu_table(ps.prop),sampledf$day)

##--------------------------------------------------------------------
# beta diversity distance to original whey
prop.ps <- transform_sample_counts(ps, function(otu) otu/sum(otu))
bray <- phyloseq::distance(prop.ps, method = "bray")
bray.m = melt(as.matrix(bray))

# remove self-comparisons
bray.m = bray.m %>%
  filter(as.character(X1) != as.character(X2)) %>%
  mutate_if(is.factor, as.character)

# get sample data (S4 error OK and expected)
sd = sample_data(prop.ps) 
class(sd) = 'data.frame'
  
sd = sd %>% select(code, day,time) %>%
  mutate_if(is.factor,as.character)

# combined distances with sample data
colnames(sd) = c("X1", "day1","time1")
bray.sd = left_join(bray.m, sd, by = "X1")

colnames(sd) = c("X2", "day2","time2")
bray.sd = left_join(bray.sd, sd, by = "X2")


#only keep comparisons to whey
bray.f = subset(bray.sd,bray.sd$day2=="whey")
bray.f = subset(bray.f,bray.f$day1 !="whey")

p3 = ggplot(bray.f, aes(x = day1, y = value, color=day1)) +
  theme_bw() +
  geom_boxplot() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#E41A1C","#377EB8","#4DAF4A","#984EA3")) +
  ylab("Bray-curtis distance") +
  xlab("Disimilarity to the whey")

#only keep comparisons to previous day of sampling
bray.f1 = subset(bray.sd,bray.sd$day1=="whey" & bray.sd$day2 == "1")
bray.f2 = subset(bray.sd,bray.sd$day1=="1" & bray.sd$day2 == "4")
bray.f3 = subset(bray.sd,bray.sd$day1=="4" & bray.sd$day2 == "7")
bray.f = bind_rows(bray.f1, bray.f2)
bray.f = bind_rows(bray.f, bray.f3)

p4 = ggplot(bray.f, aes(x = day2, y = value, color=day2)) +
  theme_bw() +
  geom_boxplot() +
  theme(legend.position = "none") +
  ylab("Bray-curtis distance") +
  scale_color_manual(values = c("#377EB8","#4DAF4A","#984EA3")) +
  xlab("Disimilarity to previous sampling time")

fig3  = plot_grid(p1, p2, p3, p4, labels = c('A','','B',''),
                  ncol = 2, nrow = 2,label_size = 12)


ggsave(fig3,file="~/Desktop/fig3.eps",height = 5, width = 6,units = "in",dpi = 300)

##--------------------------------------------------------------------
# proportion of asv and microbiome to the original whey
prop.ps <- transform_sample_counts(ps, function(otu) otu/sum(otu))

sd = sample_data(prop.ps) 
class(sd) = 'data.frame'
sd = sd %>% select(code, day) %>%
  mutate_if(is.factor,as.character)

asv = as.data.frame(t(otu_table(prop.ps)))

######bacteria
#microbiome prop
asv.whey = as.data.frame(t(asv[(asv['whey_F_filt.fastq.gz']>0),]))
s = rowSums(asv.whey)
sd$sum = s
#asv observation prop
ob = apply(asv.whey,1,function(x) sum(x > 0))
sd$ob = ob
ob.all = apply(asv,2,function(x) sum(x > 0))
sd$ob.all = ob.all
sd$ob.prop = ob/ob.all
sd = sd[1:24,]

######fungi
#microbiome prop
asv.whey = as.data.frame(t(asv[(asv['LSsumwinITS81921C03_F_filt.fastq.gz']>0),]))
s = rowSums(asv.whey)
sd$sum = s
#asv observation prop
ob = apply(asv.whey,1,function(x) sum(x > 0))
sd$ob = ob
ob.all = apply(asv,2,function(x) sum(x > 0))
sd$ob.all = ob.all
sd$ob.prop = ob/ob.all
sd = sd[1:22,]


p3 = ggplot(sd, aes(x = day, y = ob.prop*100, color=day)) +
  theme_bw() +
  geom_boxplot() +
  theme(legend.position = "none") +
  ylab("% ASV shared with whey") +
  scale_color_manual(values = c("#E41A1C","#377EB8","#4DAF4A","#984EA3")) +
  xlab("Day")
p3
p4 = ggplot(sd, aes(x = day, y = sum*100, color=day)) +
  theme_bw() +
  geom_boxplot() +
  theme(legend.position = "none") +
  ylab("% reads shared with whey") +
  scale_color_manual(values = c("#E41A1C","#377EB8","#4DAF4A","#984EA3")) +
  xlab("Day")
p4


fig4  = plot_grid(p1, p2, p3, p4, labels = c('A','','B',''),
                  ncol = 2, nrow = 2,label_size = 12)


ggsave(fig4,file="~/Desktop/fig4.eps",height = 5.5, width = 6.5,units = "in",dpi = 300)

## --------------------------------------------------------------------
# bar plot
merge = merge_samples(ps,"day")
prop.ps = transform_sample_counts(merge, function(otu) otu/sum(otu))
melt = psmelt(prop.ps)

melt1<-melt

#for bacteria
melt1$genus = paste(melt1$Genus,melt1$OTU)
melt1$genus[melt1$Abundance < 0.01] <- "<1% abun."
species = read.csv("bac_species.csv",header = TRUE)
melt1 = left_join(melt1,species,by="genus")
melt1$species <- reorder(melt1$species, -melt1$Abundance)
#for fungi
melt1$Genus <- gsub('g__', '', melt1$Genus)
melt1$Species <- gsub('s__', '', melt1$Species)
melt1$species = paste(melt1$Genus,melt1$Species)
melt1$species[is.na(melt1$Genus)] <- "unclassified"
melt1$species[melt1$Abundance < 0.01] <- "<1% abun."
melt1$species <- reorder(melt1$species, -melt1$Abundance)


bac_colors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 12))
fun_colors = c(brewer.pal(name="Set2", n = 8), brewer.pal(name="Paired", n = 10))

melt1$Sample <- as.character(melt1$Sample)
melt1$Sample <- factor(melt1$Sample, levels=c("0","whey","1","4","7"))

p1 = ggplot(melt1, aes(x=Sample, y = Abundance, fill = species)) +
  geom_bar(stat = "identity") +
  theme(axis.title.x = element_blank(),
        legend.text=element_text(size=10), 
        legend.title = element_blank()) +
  scale_fill_manual(values = c(bac_colors)) +
  ylab("Relative abundance")


p2 = ggplot(melt1, aes(x=Sample, y = Abundance, fill = species)) +
  geom_bar(stat = "identity") +
  theme(axis.title.x = element_blank(),
        legend.text=element_text(size=10), 
        legend.title = element_blank()) +
  scale_fill_manual(values = c(fun_colors)) +
  ylab("Relative abundance")


fig5  = plot_grid(p1, p2,labels = c('A','B'),
                  ncol = 1, nrow = 2,label_size = 12)


ggsave(fig5,file="~/Desktop/fig5.eps",height = 10, width = 10,units = "in",dpi = 300)
## --------------------------------------------------------------------
#sloan's neutral community model
#https://github.com/Russel88/MicEco
library(Hmisc)
library(bbmle)
library(optimx)

taxa = as.data.frame(tax_table(ps))
taxa$asv = rownames(taxa)

ps.d1 <- subset_samples(ps, day == "1")
ps.d1 = filter_taxa(ps.d1, function(x) mean(x) > 0, TRUE)
#ps.d1 = subset_samples(ps.d1, time == "june")
#ps.d1 = subset_samples(ps.d1, time == "july")
d1<-otu_table(ps.d1)
pre.d1 <- neutral.fit(d1)

pre.d1[[2]]$effect[pre.d1[[2]]$freq < pre.d1[[2]]$Lower ] <- "low"
pre.d1[[2]]$effect[pre.d1[[2]]$freq > pre.d1[[2]]$Upper ] <- "high"
pre.d1[[2]]$effect[is.na(pre.d1[[2]]$effect)] <- "neutral"

p6.1 = ggplot(pre.d1[[2]],aes(log10(p),freq)) + 
  geom_point(aes(color = effect)) + 
  theme_bw() +
  scale_color_manual(values = c("red","blue","grey")) +
  geom_line(data=pre.d1[[2]],aes(log10(p),Lower),linetype = "dashed") + 
  geom_line(data=pre.d1[[2]],aes(log10(p),Upper),linetype = "dashed") + 
  xlab("Average relative abundance (Log10)") + ylab("Frequency") +
  geom_label(x = -2, y = 0, label = "R2=0.018") +
  ggtitle("d1")

p6.4 = ggplot(pre.d1[[2]],aes(log10(p),freq)) + 
  geom_point(aes(color = effect)) + 
  theme_bw() +
  scale_color_manual(values = c("red","blue","grey")) +
  geom_line(data=pre.d1[[2]],aes(log10(p),Lower),linetype = "dashed") + 
  geom_line(data=pre.d1[[2]],aes(log10(p),Upper),linetype = "dashed") + 
  xlab("Average relative abundance (Log10)") + ylab("Frequency") +
  geom_label(x = -2, y = 0, label = "R2=0.16") +
  ggtitle("d1")


ps.d4 <- subset_samples(ps, day == "4")
ps.d4 = filter_taxa(ps.d4, function(x) mean(x) > 0, TRUE)
d4<-otu_table(ps.d4)
pre.d4 <- neutral.fit(d4)

pre.d4[[2]]$effect[pre.d4[[2]]$freq < pre.d4[[2]]$Lower ] <- "low"
pre.d4[[2]]$effect[pre.d4[[2]]$freq > pre.d4[[2]]$Upper ] <- "high"
pre.d4[[2]]$effect[is.na(pre.d4[[2]]$effect)] <- "neutral"

p6.2 = ggplot(pre.d4[[2]],aes(log10(p),freq)) + 
  geom_point(aes(color = effect)) + 
  theme_bw() +
  scale_color_manual(values = c("red","blue","grey")) +
  geom_line(data=pre.d4[[2]],aes(log10(p),Lower),linetype = "dashed") + 
  geom_line(data=pre.d4[[2]],aes(log10(p),Upper),linetype = "dashed") + 
  xlab("Average relative abundance (Log10)") + ylab("Frequency") +
  geom_label(x = -2, y = 0, label = "R2=-0.07") +
  ggtitle("d4")

p6.5 = ggplot(pre.d4[[2]],aes(log10(p),freq)) + 
  geom_point(aes(color = effect)) + 
  theme_bw() +
  scale_color_manual(values = c("red","blue","grey")) +
  geom_line(data=pre.d4[[2]],aes(log10(p),Lower),linetype = "dashed") + 
  geom_line(data=pre.d4[[2]],aes(log10(p),Upper),linetype = "dashed") + 
  xlab("Average relative abundance (Log10)") + ylab("Frequency") +
  geom_label(x = -2, y = 0, label = "R2=-0.07") +
  ggtitle("d4")

ps.d7 <- subset_samples(ps, day == "7")
ps.d7 = filter_taxa(ps.d7, function(x) mean(x) > 0, TRUE)
d7<-otu_table(ps.d7)
pre.d7 <- neutral.fit(d7)

pre.d7[[2]]$effect[pre.d7[[2]]$freq < pre.d7[[2]]$Lower ] <- "low"
pre.d7[[2]]$effect[pre.d7[[2]]$freq > pre.d7[[2]]$Upper ] <- "high"
pre.d7[[2]]$effect[is.na(pre.d7[[2]]$effect)] <- "neutral"

p6.3 = ggplot(pre.d7[[2]],aes(log10(p),freq)) + 
  geom_point(aes(color = effect)) + 
  theme_bw() +
  scale_color_manual(values = c("red","blue","grey")) +
  geom_line(data=pre.d7[[2]],aes(log10(p),Lower),linetype = "dashed") + 
  geom_line(data=pre.d7[[2]],aes(log10(p),Upper),linetype = "dashed") + 
  xlab("Average relative abundance (Log10)") + ylab("Frequency") +
  geom_label(x = -2, y = 0, label = "R2=-0.21") +
  ggtitle("d7")

p6.6 = ggplot(pre.d7[[2]],aes(log10(p),freq)) + 
  geom_point(aes(color = effect)) + 
  theme_bw() +
  scale_color_manual(values = c("red","grey")) +
  geom_line(data=pre.d7[[2]],aes(log10(p),Lower),linetype = "dashed") + 
  geom_line(data=pre.d7[[2]],aes(log10(p),Upper),linetype = "dashed") + 
  xlab("Average relative abundance (Log10)") + ylab("Frequency") +
  geom_label(x = -2, y = 0, label = "R2=0.01") +
  ggtitle("d7")

fig6 = plot_grid(p6.1, p6.2, p6.3, p6.4, p6.5, p6.6, 
               labels = c('A','','','B','',''), 
               ncol = 3,nrow = 2,label_size = 12)


ggsave(fig6,file="~/Desktop/fig6.eps",height = 6, width = 12,units = "in",dpi = 300)


## --------------------------------------------------------------------
#bNTI null model
library(picante)
ps.d = subset_samples(ps,day=="0")
otu = as.data.frame(otu_table(ps.d))
phylo = phy_tree(ps.d)

match.phylo.otu = match.phylo.data(phylo, t(otu))

## calculate empirical betaMNTD
beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T));

## calculate randomized betaMNTD
beta.reps = 999
rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.otu$data),ncol(match.phylo.otu$data),beta.reps))

for (rep in 1:beta.reps) {
  
rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.otu$data),taxaShuffle(cophenetic(match.phylo.otu$phy)),abundance.weighted=T,exclude.conspecifics = F))
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(match.phylo.otu$data),ncol=ncol(match.phylo.otu$data))

for (columns in 1:(ncol(match.phylo.otu$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.otu$data)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
}

rownames(weighted.bNTI) = colnames(match.phylo.otu$data);
colnames(weighted.bNTI) = colnames(match.phylo.otu$data);

hist(weighted.bNTI)

bac = read.csv("./bNTI.csv")
bac = melt(bac)
bac$variable = gsub("d","",bac$variable)
bac$variable = as.numeric(bac$variable)
ggplot(bac, aes(x=variable, y=value)) + 
  geom_point(color="grey") + 
  geom_smooth(method=lm) +
  theme_bw() +
  geom_hline(yintercept=2, linetype="solid") +
  geom_hline(yintercept=-2, linetype="solid") +
  ylim(-2,2)

## --------------------------------------------------------------------
#Raup_Crick null model




## --------------------------------------------------------------------
#SpiecEasi: The stability of biofilm microbial community will be determined based on microbial network structure.

ps.d1 <- subset_samples(ps, day == "1")
ps.d1 <- filter_taxa(ps.d1, function(x) sum(x) > 0, TRUE) 
d1 = data.matrix(otu_table(ps.d1))

sparcc.d1 <- sparcc(d1)
## Define arbitrary threshold for SparCC correlation matrix for the graph
sparcc.graph <- abs(sparcc.d1$Cor) >= 0.3
diag(sparcc.graph) <- 0
library(Matrix)
sparcc.graph <- Matrix(sparcc.graph, sparse=TRUE)
## Create igraph objects
ig.sparcc <- adj2igraph(sparcc.graph)
edge_density(ig.sparcc, loops = FALSE)

library(igraph)
## set size of vertex proportional to clr-mean
vsize <- rowMeans(clr(d1, 1))+6

plot(ig.sparcc, vertex.size=vsize, vertex.label=NA, main="sparcc")

## --------------------------------------------------------------------
#extract data for picrust
# whole community
ps.biofilm = subset_samples(ps,day!="whey")
ps.f = filter_taxa(ps.biofilm, function(x) mean(x) > 0, TRUE)

otu = as.data.frame(otu_table(ps.f))
seq = as.data.frame(refseq(ps.f))
write.csv(otu,file="~/Desktop/Dissertation/aim2/analysis/picrust/otu.csv")
write.csv(seq,file="~/Desktop/Dissertation/aim2/analysis/picrust/seq.fna")

#partitions by non-neutral
ps.f= subset_samples(ps,day!="0")
ps.f= subset_samples(ps.f,day!="whey")
ps.f = filter_taxa(ps.f, function(x) mean(x) > 0, TRUE)

otu = as.data.frame(otu_table(ps.f))
seq = as.data.frame(refseq(ps.f))
write.csv(otu,file="~/Desktop/Dissertation/aim2/analysis/picrust/otu.csv")
write.csv(seq,file="~/Desktop/Dissertation/aim2/analysis/picrust/seq.fna")

#manual combining samples by partition
tb = read.table("/Users/langsun/Desktop/Dissertation/aim2/analysis/picrust/from\ sloan\ model/otu_relative\ abundance\ by\ partitions.txt")

rd_tb = tb*2000
rd_tb = round_df(rd_tb, 0)
write.table(rd_tb,file="~/Desktop/Dissertation/aim2/analysis/picrust/otu.txt",sep = "\t")


## --------------------------------------------------------------------
# https://huttenhower.sph.harvard.edu/galaxy/
#run categorize_by_function.r,change "the string separate by ;" to [1] to get level 3 category

dvd =read.table(file =  "/Users/langsun/Desktop/Dissertation/aim2/analysis/picrust/total biofilm microbial community/predicted.tsv",header=TRUE, sep="\t", row.names=1)

kegg_brite_map = read.table("/Users/langsun/Desktop/Dissertation/aim2/analysis/picrust/picrust1_KO_BRITE_map.tsv",header=TRUE, sep="\t", quote = "",stringsAsFactors = FALSE,comment.char="", row.names=1)

bacteria_ko_L3 = categorize_by_function_l1(dvd, kegg_brite_map)
ko_sorted <- bacteria_ko_L3[rownames(bacteria_ko_L3), ]

ko_sorted = as.matrix(ko_sorted)
ko_sorted = make_relative(ko_sorted)
ko_sorted = as.data.frame(ko_sorted)
ko_sorted$pathways = rownames(ko_sorted)

g1 = read.csv(file = "~/Desktop/Dissertation/aim2/analysis/picrust/total biofilm microbial community/ko_3_hm/1_4.csv")
g2 = read.csv(file = "~/Desktop/Dissertation/aim2/analysis/picrust/total biofilm microbial community/ko_3_hm/1_7.csv")
g3 = read.csv(file = "~/Desktop/Dissertation/aim2/analysis/picrust/total biofilm microbial community/ko_3_hm/4_7.csv")

g = full_join(g,g3,by="pathways")

pathway = read.csv(file = "~/Desktop/Dissertation/aim2/analysis/picrust/total biofilm microbial community/ko_3_hm/pathway.csv")


## --------------------------------------------------------------------
# aldex2 on picrust comparison

round_df <- function(x, digits) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}

rd_df = round_df(ko_sorted, 0)

#apply aldex
conds <- c(rep("d4",6),rep("d7",6))
aldex.d1 <- aldex(rd_df, conds, mc.samples=128, test="t", effect=TRUE, 
                 include.sample.summary=FALSE, denom="all", verbose=FALSE)


x.d = aldex.d1
x.d$effect.sig[abs(x.d$effect) <= 1] <- "non.sig"
x.d$effect.sig[x.d$effect > 1] <- "sig.pos"
x.d$effect.sig[x.d$effect < -1] <- "sig.neg"

# add ko level 1 
x.d$pathways = rownames(x.d)
x.d.ko1 = left_join(x.d,pathway,by="pathways")
x.d.ko1$Ko_1[is.na(x.d.ko1$Ko_1)] <- "Other Eukaryotic Pathways"
x.d.ko1 = subset(x.d.ko1,x.d.ko1$Ko_1!="Other Eukaryotic Pathways")

hm = read.csv("~/Desktop/Dissertation/aim2/analysis/picrust/total biofilm microbial community/ko_3_hm/ko_hm.csv")

hm_matrix = hm[,4:6]
rownames(hm_matrix) = hm$pathways

library(circlize)
col_fun = colorRamp2(c(-18,-6, 6), c("blue","white","red"))
my_group <- as.numeric(as.factor(hm$Ko_1))

Heatmap(hm_matrix,row_split = my_group,col = col_fun,
        row_gap = unit(5, "mm"),
        row_names_gp = gpar(fontsize = c(6)),
        row_title_gp = gpar(fill = c("red", "blue", "green","yellow"),fontsize = c(10)))


hm_matrix=as.matrix(hm_matrix)
my_group <- as.numeric(as.factor(hm$Ko_1))
rowSide <- brewer.pal(9, "Set1")[my_group]
heatmap(hm_matrix,RowSideColors=rowSide)

p = ggplot(x.d.ko1, aes(x = effect,y = Ko_1)) +
  geom_jitter(aes(color = effect.sig),position = position_jitter(width=0.1, height=0)) +
  scale_color_manual(values = c("non.sig" = "grey", "sig.pos" = "red","sig.neg" = "blue")) +
  xlim(-15, 15) 

## ------------------------------------------------------------------
# PCoA on KO 
dvd =read.table(file =  "/Users/langsun/Desktop/Dissertation/aim2/analysis/picrust/total biofilm microbial community/ko_3_hm/bio_dev/predicted.tsv",header=TRUE, sep="\t", row.names=1)

metadata=read.csv("picrust/total biofilm microbial community/ko_3_hm/bio_dev/metadata.csv",header = T,row.names=1)

ko.ps <- phyloseq(otu_table(bacteria_ko_L3,taxa_are_rows=TRUE), 
               sample_data(metadata))

ko.prop <- transform_sample_counts(ko.ps, function(otu) otu/sum(otu))
ko.ord.pcoa.bray <- ordinate(ko.prop, method="PCoA", distance="bray")

p = plot_ordination(ko.prop, ko.ord.pcoa.bray, color="day")

p$data$day <- as.character(p$data$day)
p$data$day <- factor(p$data$day, levels=c("0","1","4","7"))

p1 = p + theme_bw() + 
  theme(legend.title=element_blank(),
        legend.text = element_text(size=10)) +
  scale_color_manual(values = c("#E41A1C","#377EB8","#4DAF4A","#984EA3"))

## ------------------------------------------------------------------
# enriched pathways
enrich=read.csv("picrust/total biofilm microbial community/ko_3_hm/enriched.csv",header = T)

ra_ko.prop = as.matrix(otu_table(ko.prop))
ra_ko.prop = prop.table(ra_ko.prop, 1)

ra_ko = phyloseq(otu_table(ra_ko.prop,taxa_are_rows=TRUE), 
         sample_data(metadata))

melt.kegg = psmelt(ra_ko)
colnames(melt.kegg) = c("Sample","pathways","Abundance","sample","day","time")

melt.kegg1 = inner_join(melt.kegg, enrich,by="pathways")

melt.kegg1$pathways = factor(melt.kegg1$pathways,levels=unique(enrich$pathways))

p1 = ggplot(melt.kegg1, aes(x=pathways, y = Abundance, fill = sample)) +
  geom_bar(stat = "identity") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size=6),
        legend.text=element_text(size=10), 
        legend.title = element_blank(),
        legend.position = "bottom") +
  xlab("Relative abundance") +
  scale_fill_manual(values = c("#456355","#FCD16B","#D3DDDC")) +
  facet_grid(enriched~.,switch = "both", 
             scales = "free", space = "free") + 
  coord_flip()

ggsave(p1,filename = "~/Desktop/picrust.eps",height = 10,width = 6)

leg <- ggplot(melt.kegg1, aes(x = 0, y = pathways)) + 
  geom_point(aes(color = Ko_1), shape = 15, size = 3, show.legend = F) + 
  theme_classic()+
  theme(axis.title = element_blank(), axis.line = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(),legend.position = "bottom",
        plot.margin = unit(c(0,0,0,0), "cm"),
        strip.text = element_blank())+
  scale_color_manual(values = c("#E41A1C","#FF7F00","#377EB8","#4DAF4A","#984EA3")) +
  facet_grid(enriched~.,switch = "both", 
             scales = "free", space = "free")

ggsave(leg,filename = "~/Desktop/anno.eps",height = 10,width = 20)

p1 + annotation_custom(ggplotGrob(leg),ymin = 1, ymax = 1.07) 

