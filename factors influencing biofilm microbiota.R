setwd("~/Desktop/Dissertation/aim3/analysis/")

##--------------------------------------------------------------------
#load pyloseq
load("~/Desktop/Dissertation/aim3/bacteria/phyloseq.RData")
load("~/Desktop/Dissertation/aim3/fungi/phyloseq.RData")

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
library(ggpubr)
library(car)
library(rcompanion)
library(lme4)
library(lmerTest)
library(FSA)
library(tidyr)
library(ggforce)
##--------------------------------------------------------------------
# sample filtering
ps@sam_data$group2 = interaction(ps@sam_data$condition,
                                 ps@sam_data$milktyp, 
                                 sep=  "_", lex.order = TRUE)
ps@sam_data$trial <- as.character(ps@sam_data$trial)
ps@sam_data$trial <- factor(ps@sam_data$trial, levels=c("1","4","7"))
ps@sam_data$production2 <- interaction(ps@sam_data$trial,
                                       ps@sam_data$production, 
                                       sep=  "_", lex.order = TRUE)
ps@sam_data$production3 <- interaction(ps@sam_data$condition,
                                       ps@sam_data$trial, 
                                       sep=  "_", lex.order = TRUE)



ps.f = filter_taxa(ps, function(x) mean(x) > 0, TRUE)
ps.f = subset_samples(ps.f,group=="biofilm")
##--------------------------------------------------------------------
# alpha diversity

###########biofilm
p = plot_richness(ps.f, x="trial", 
                  measures=c("Shannon"),
                  color = "production")
f1a = p + geom_point(aes(shape = group2)) + theme_bw() + 
  theme(axis.title=element_blank(),
    #legend.title = element_blank(),
    panel.background = element_rect(fill = 'white', 
                                    colour = 'black')) +
  geom_signif(comparisons = list(c("1", "4"),c("4","7"),c("1","7")),
              map_signif_level=TRUE,step_increase = 0.1,color = "black",
              size = 0.3,textsize = 2) 
f1a
ggsave(f1a, file="~/Desktop/f1a.eps",height = 3, width = 4,units = "in",dpi = 300)

p = plot_richness(ps.f, x="milktype", 
                  measures=c("Shannon"),
                  color = "group2")
f1a = p + geom_point() + theme_bw() + 
  theme(axis.title=element_blank(),
        #legend.title = element_blank(),
        panel.background = element_rect(fill = 'white', 
                                        colour = 'black')) +
  geom_signif(comparisons = list(c("raw", "heated")),
              map_signif_level=TRUE,step_increase = 0.1,color = "black",
              size = 0.3,textsize = 2) +
  facet_grid(.~group)
f1a
ggsave(f1a, file="~/Desktop/f1a.eps",height = 3, width = 6,units = "in",dpi = 300)

###########cheese
p = plot_richness(ps.f, x="trial", 
                  measures=c("Shannon"),
                  color = "production")
f3a = p + geom_point(aes(shape = group2)) + theme_bw() + 
  theme(axis.title=element_blank(),
        #legend.title = element_blank(),
        panel.background = element_rect(fill = 'white', 
                                    colour = 'black')) +
  geom_signif(comparisons = list(c("1", "4"),c("4","7"),c("1","7")),
              map_signif_level=TRUE,step_increase = 0.1,color = "black",
              size = 0.3,textsize = 2) +
  facet_grid(.~condition,switch = "both", 
             scales = "free", space = "free")
f3a
ggsave(f3a, file="~/Desktop/f3a.eps",height = 3, width = 4,units = "in",dpi = 300)
##--------------------------------------------------------------------
# beta diversity clustering
prop.ps <- transform_sample_counts(ps.f, function(otu) otu/sum(otu))
ord.pcoa.bray <- ordinate(prop.ps, method="PCoA", distance="bray")

#f1b
f1b = plot_ordination(prop.ps, ord.pcoa.bray,
                     color="production",shape = "group2") + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'black')) +
  stat_ellipse(aes(color=production,group=production),type = "t")
f1b
ggsave(f1b, file="~/Desktop/f1b.eps",
       height = 3, width = 5,units = "in",dpi = 300)


#f1c
col1 = c("#FCBBA1","#FCBBA1","#FB6A4A","#CB181D","#C6DBEF","#C6DBEF","#6BAED6","#08519C")
col2 = c("#FCBBA1","#FB6A4A","#CB181D","#C6DBEF","#6BAED6","#08519C")
f1c = plot_ordination(prop.ps, ord.pcoa.bray,
                      color="production3",shape = "milktype") + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'black')) +
  scale_color_manual(values = col1) +
  stat_ellipse(aes(group=condition,color=condition),type = "t")
f1c
ggsave(f1c, file="~/Desktop/f1c.eps",
       height = 4, width = 6,units = "in",dpi = 300)

#f3b
prop.ps = subset_samples(prop.ps,trial!="1")
ord.pcoa.bray <- ordinate(prop.ps, method="PCoA", distance="bray")
col2 = c("#FCBBA1","#FB6A4A","#CB181D","#C6DBEF","#6BAED6","#08519C")
f3b = plot_ordination(prop.ps, ord.pcoa.bray,
                      color="group2",shape = "trial") + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'black')) +
  scale_color_manual(values = col2) +
  stat_ellipse(aes(group=condition,color=condition),type = "t")
f3b                   

#s1
fs1 = plot_ordination(prop.ps, ord.pcoa.bray,
                      color="milktype") + 
  theme(panel.background = element_rect(fill = 'white', 
                                        colour = 'black'))
fs1
ggsave(fs1, file="~/Desktop/fs1.eps",
       height = 3, width = 5,units = "in",dpi = 300)

##--------------------------------------------------------------------
# beta diversity 
# Calculate bray curtis distance matrix
prop.ps <- transform_sample_counts(ps.f, function(otu) otu/sum(otu))
beta <- phyloseq::distance(prop.ps, method = "bray")
# make a data frame from the sample_data
sampledf <- data.frame(sample_data(prop.ps))

# Adonis test
adonis2(beta ~ condition, data = sampledf)
adonis2(beta ~ milktype, data = sampledf)
adonis2(beta ~ trial, data = sampledf)

sampledf$production = sampledf$trial
sampledf$production = gsub("4","after", sampledf$production)
sampledf$production = gsub("7","after", sampledf$production)
adonis2(beta ~ production, data = sampledf)

#pairwise for biofilm
adonis2(beta ~ trial*condition, data = sampledf)
sampledf$con = paste(sampledf$trial,sampledf$condition)
pairwise.adonis(otu_table(prop.ps),sampledf$con)

adonis2(beta ~ milktype*condition, data = sampledf)
sampledf$con = paste(sampledf$milktype,sampledf$condition)
pairwise.adonis(otu_table(prop.ps),sampledf$con)


## --------------------------------------------------------------------
#bar plot
ps.f@sam_data$group3 = interaction(ps.f@sam_data$group2,
                                   ps.f@sam_data$trial, 
                                   ps.f@sam_data$group,
                                   sep=  "_", lex.order = TRUE)
ps.f@sam_data$group3= as.character(ps.f@sam_data$group3) 
merge = merge_samples(ps.f,"group3")
prop.ps <- transform_sample_counts(merge, function(otu) otu/sum(otu))
melt <- psmelt(prop.ps)

melt1<-melt
melt1 = separate(data = melt1, col = Sample, into = c("condition","milk", "trial","group"),sep = "_")
melt2 = subset(melt1,melt1$group!="milk")


#######################for bacteria
melt2$genus = paste(melt2$Genus,melt2$OTU)
melt2$genus[melt2$Abundance < 0.01] <- "<1% abun."
species = read.csv("bac_species.csv",header = TRUE)
melt2 = left_join(melt2,species,by="genus")
#melt2$species <- reorder(melt2$species, -melt2$Abundance)

bac.col = c("#E41A1C","#FC4E2A",
            "#FEB24C","#FD8D3C","grey",
            "#9E9AC8","#6A51A3",
            "#E5F5E0","#C7E9C0", "#A1D99B", "#74C476",
            "#238B45","#00441B", #two dominant lactobacillus
            "#6BAED6","#377EB8",
            "#FFEDA0", "#FED976",
            "#41AB5D",
            "#969696")

p.bac = ggplot(melt2, aes(x=trial, y = Abundance, fill = updated_species)) +
  geom_bar(stat = "identity") +
  theme(#axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.text=element_text(size=6), 
    legend.key.size = unit(0.13, "in"),
    legend.title = element_blank(),
    strip.text.x = element_text(size = 6)) +
  facet_grid(group~condition+milk,switch = "both", 
             scales = "free", space = "free") +
  scale_fill_manual(values = bac.col) +
  ylab("Relative Abundance")
p.bac
ggsave(p.bac, file="~/Desktop/p.bac.eps",height = 4, width = 6.8,units = "in",dpi = 300)



################for fungi
melt1$Genus <- gsub('g__', '', melt1$Genus)
melt1$Species <- gsub('s__', '', melt1$Species)
melt1$species = paste(melt1$Genus,melt1$Species)
melt1$species[melt1$Abundance < 0.01] <- "<1% abun."

fun.col = c(brewer.pal(name="Set3", n = 12))

p.fun = ggplot(melt1, aes(x=trial, y = Abundance, fill = species)) +
  geom_bar(stat = "identity") +
  theme(#axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.text=element_text(size=6), 
    legend.key.size = unit(0.13, "in"),
    legend.title = element_blank(),
    strip.text.x = element_text(size = 6)) +
  facet_grid(group~condition+milk,switch = "both", 
             scales = "free", space = "free") +
  scale_fill_manual(values = fun.col) +
  ylab("Relative Abundance")

ggsave(p.fun, file="~/Desktop/fig6.eps",height = 4, width = 6.8,units = "in",dpi = 300)


## --------------------------------------------------------------------
#statistics differentially abundant taxa
####bacteria
prop.ps = transform_sample_counts(ps.f, function(otu) otu/sum(otu))
ps.melt = subset_samples(prop.ps,trial=="7")
ps.melt = subset_samples(prop.ps,production=="after")
ps.melt = filter_taxa(ps.melt, function(x) sum(x) > 0.01, TRUE)
ps.melt
melt2 = psmelt(ps.melt)
melt2$transform = asin(sqrt(melt2$Abundance))
melt2$genus = paste(melt2$Genus,melt2$OTU)
melt2 = left_join(melt2,species,by="genus")


###fungi
prop.ps = transform_sample_counts(ps.f, function(otu) otu/sum(otu))
ps.melt = tax_glom(prop.ps,taxrank = "Species")
ps.melt = subset_samples(ps.melt,production=="after")
ps.melt = filter_taxa(ps.melt, function(x) sum(x) > 0.01, TRUE)
ps.melt
melt2 = psmelt(ps.melt)
melt2$transform = asin(sqrt(melt2$Abundance))
melt2$genus = paste(melt2$Genus,melt2$Species)


variable = as.data.frame(unique(melt2$genus))
colnames(variable) <- c("v")

for (row in 1:19) {
  tryCatch({
    v = variable[row,"v"] 
    print(v)
    
    plotNormalHistogram((subset(melt2,melt2$genus==v))$transform,main = v)
    test = t.test(Abundance ~ condition,
                   data = subset(melt2,melt2$genus==v))
    print(test)
    
    test1 = t.test(transform ~ condition,
                   data = subset(melt2,melt2$genus==v))
    print(test1)
    
    test2 = wilcox.test(transform ~ condition,
                        data = subset(melt2,melt2$genus==v))
    print(test2)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

random_tree = rtree(ntaxa(ps.f), rooted=TRUE, tip.label=taxa_names(ps.f))
plot(random_tree)

## --------------------------------------------------------------------
#statistics differentially abundant taxa using deseq2
###for fungi
ps.deseq = tax_glom(ps.f,taxrank = "Species")
ps.deseq = subset_samples(ps.deseq,production=="after")
ps.deseq
otu_table(ps.deseq) <- otu_table(ps.deseq) + 1

###for bacteria
ps.deseq = ps.f
ps.deseq = subset_samples(ps.f,production=="after")
ps.deseq = subset_samples(ps.deseq,milktype=="raw")
ps.deseq
otu_table(ps.deseq) <- otu_table(ps.deseq) + 1

#run deseq2
deseq = phyloseq_to_deseq2(ps.deseq, ~ condition)
deseq = DESeq(deseq, test="Wald", fitType="parametric")

#result
res = results(deseq, cooksCutoff = FALSE)
res.df = as.data.frame(res@listData)

sigtab = cbind(as(res, "data.frame"), as(tax_table(ps.deseq)[rownames(res), ], "matrix")) 
sigtab$OTU = rownames(sigtab)
sigtab$genus = paste(sigtab$Genus,sigtab$OTU)
sigtab = left_join(sigtab,species,by="genus")

## --------------------------------------------------------------------
#SpiecEasi: The stability of biofilm microbial community will be determined based on microbial network structure.
library(SpiecEasi)
ps.bio <- subset_samples(ps, group == "biofilm")
ps.summer <- subset_samples(ps.bio, condition == "summer")
ps.winter <- subset_samples(ps.bio, condition == "winter")

summer <- data.matrix(otu_table(ps.summer))
summer.mb <- spiec.easi(summer, method='mb', lambda.min.ratio=1e-2,
                          nlambda=20, pulsar.params=list(rep.num=50))

winter.mb <- spiec.easi(ps.winter, method='mb', lambda.min.ratio=1e-2,
                        nlambda=20, pulsar.params=list(rep.num=50))

## Create igraph objects
library(igraph)
ig.summer <- adj2igraph(getRefit(summer.mb))
ig.winter <- adj2igraph(getRefit(winter.mb))

edge_density(ig.mb, loops = FALSE)
edge_density(ig.mb, loops = FALSE)

## set size of vertex proportional to clr-mean
vsize <- rowMeans(clr(d1, 1))+6

plot(ig.mb, vertex.size=vsize, vertex.label=NA, main="sparcc")

