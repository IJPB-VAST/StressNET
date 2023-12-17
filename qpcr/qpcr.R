library("reshape2")
library("agricolae")
library("ggplot2")
library("ggrepel")
library('gtools')
library('circlize')
library("ComplexHeatmap")

library("tibble")
library("stringr")
library("ggpubr")
library("kableExtra")
library("gridExtra")
library("tidyr")
library("dplyr")
library("patchwork")
library("plyr")
library('cowplot')


counts<- read.table("RNA.counts.scaled.log2.tsv",check.names = F)
counts<- as.matrix(counts)
counts<- counts[,!(colnames(counts) %in% c("Bur.0_W1N1_P3ID34_b"))]
ymelt <- reshape2::melt(counts) %>% separate(Var2, c("Genotype","condition","phenoscope","rep"), sep = "_") 



ymelt <- ymelt %>%
  mutate(Genotype = case_when(
    Genotype == "Col.0" ~ "Col-0",
    Genotype == "Cvi.0" ~ "Cvi-0",
    Genotype == "Sha" ~ "Sha",
    Genotype == "Tsu.0" ~ "Tsu-0",
    Genotype == "Bur.0" ~ "Bur-0",
  )) %>% 
  mutate(condition = case_when(
    condition == "W0N0" ~ "W-N-",
    condition == "W0N1" ~ "W-N+",
    condition == "W1N0" ~ "W+N-",
    condition == "W1N1" ~ "W+N+",
  ))
ymelt$condition <- factor(ymelt$condition,levels= c("W+N+","W-N+","W+N-","W-N-")) 
ymelt$Genotype <- factor(ymelt$Genotype,levels= c("Col-0","Cvi-0","Sha","Tsu-0","Bur-0"))  
ymelt$geneID <- factor(ymelt$Var1)  
ymelt$counts <-  ymelt$value  



bur0 <- readr::read_csv('Bur0_ALL.csv')%>% mutate(Geno = 'Bur-0')
col0 <- readr::read_csv('Col0_ALL.csv')%>% mutate(Geno = 'Col-0')
Cvi0 <- readr::read_csv('Cvi0_ALL.csv')%>% mutate(Geno = 'Cvi-0')
sha <- readr::read_csv('Sha_ALL.csv')%>% mutate(Geno = 'Sha')
Tsu0 <- readr::read_csv('Tsu0_ALL.csv')%>% mutate(Geno = 'Tsu-0')


dataframe <- rbind(bur0,Cvi0,col0,sha,Tsu0) %>% as.data.frame()
colnames(dataframe) <- gsub(' ','',colnames(dataframe) )

dataframe <- dataframe %>%
  mutate(geneID = case_when(
    Target == "act8" ~ "AT1G49240",
    Target == "gapd" ~ "AT1G16300",
    Target == "bor5" ~ "AT1G74810",
    Target == "cepd2" ~ "AT2G47880",
    Target == "cer1" ~ "AT1G02205",
    Target == "hho3" ~ "AT1G25550",
    Target == "nit4" ~ "AT5G22300",
    Target == "nrt2.5" ~ "AT1G12940",
    Target == "nudx" ~ "AT5G19470",
    Target == "pyl6" ~ "AT2G40330",
    Target == "roxy13" ~ "AT4G15680",
    Target == "roxy16" ~ "AT1G03020",
    Target == "roxy2" ~ "AT5G14070",
    Target == "roxy4" ~ "AT3G62950",
    Target == "roxy8" ~ "AT3G62960",
    Target == "tar" ~ "AT1G34060",
    TRUE ~ NA_character_   
  )) %>% 
  mutate(Cond = case_when(
    BiologicalGroup == "W0N0" ~ "W-N-",
    BiologicalGroup == "W0N1" ~ "W-N+",
    BiologicalGroup == "W1N0" ~ "W+N-",
    BiologicalGroup == "W1N1" ~ "W+N+",
  ))

# meanExpression 

# newqpcr <- dataframe[,c('Target','BiologicalGroup','MeanCq','CqSD','CqSEM','Geno','geneID','Cond')] %>%
#   group_by(Geno,BiologicalGroup) %>% do(mutate(., deltaCq= MeanCq- (MeanCq[Target == "act8"] + MeanCq[Target == "gapd"])/2  ))%>% ungroup()%>%
#   group_by(Geno,Target) %>% do(mutate(., deltaCon = deltaCq[BiologicalGroup=='W1N1']  ))  %>% ungroup()%>%
#   mutate(meanSEM = CqSEM  )%>%
#   mutate(reference = "ref")%>%
#   mutate(meanExpression = 2^-(deltaCq - deltaCon)  ) 
  
newqpcr <- dataframe[,c('Target','BiologicalGroup','Expression','ExpressionSEM','Geno','geneID','Cond')] %>%
           mutate(meanExpression = .$Expression  )%>%
           mutate(meanSEM = .$ExpressionSEM  )

  
newqpcr$Cond <- factor(newqpcr$Cond,levels= c("W+N+","W-N+","W+N-","W-N-")) 
newqpcr$Geno <- factor(newqpcr$Geno,levels= c("Col-0","Cvi-0","Sha","Tsu-0","Bur-0"))  

gene_aliases <- read.table( "gene_aliases_20181231_summary.tsv",header = T)
newqpcr$geneSymbol <- gene_aliases$symbols[match(newqpcr$geneID ,gene_aliases$name)]


gene = "AT1G02205"
title = na.omit(newqpcr[newqpcr$geneID== gene ,])$geneSymbol[1]
p1 <- ggbarplot(  na.omit(newqpcr[newqpcr$geneID== gene ,]) , x = "Cond", y = "meanExpression",  color = "Cond",palette = c("blue","orange","purple","red"),   lab.vjust = -1.6) +  facet_wrap(. ~ Geno ,ncol = 5,scales = "fixed",drop = FALSE) +theme_bw()  +   geom_errorbar( aes(x=Cond, ymin=meanExpression-meanSEM, ymax=meanExpression + meanSEM), width=0.4, colour="black", alpha=0.9, linewidth=0.3) + theme_bw() + ylab("meanExpression ± meanSEM") + ggtitle(title) + theme(axis.text.x = element_text(angle=90))
p2 <-   ggboxplot(ymelt[ymelt$geneID==gene ,], x = "condition", y = "counts", color = "condition", add = "jitter",palette = c("blue","orange","purple","red"))  +  facet_wrap(. ~ Genotype,ncol = 5,scales = "fixed",drop = FALSE) +  rotate_x_text(90)  + ylab(" Normalized Counts") + stat_compare_means(comparisons = list(c("W+N+","W-N+"),c("W+N+","W+N-"),c("W+N+","W-N-")), label = "p.signif",method = "t.test") + theme_bw()  + ggtitle(title)+ theme(axis.text.x = element_text(angle=90))

gene = "AT1G03020"
title = na.omit(newqpcr[newqpcr$geneID== gene ,])$geneSymbol[1]

p3 <- ggbarplot(  na.omit(newqpcr[newqpcr$geneID== gene ,]), x = "Cond", y = "meanExpression",  color = "Cond",palette = c("blue","orange","purple","red"),   lab.vjust = -1.6) +  facet_wrap(. ~ Geno  ,ncol = 5,scales = "fixed",drop = FALSE) +theme_bw()  +   geom_errorbar( aes(x=Cond, ymin=meanExpression-meanSEM, ymax=meanExpression + meanSEM), width=0.4, colour="black", alpha=0.9, linewidth=0.3) + theme_bw() + ylab("meanExpression ± meanSEM") + ggtitle(title) + theme(axis.text.x = element_text(angle=90))
p4 <-   ggboxplot(ymelt[ymelt$geneID==gene ,], x = "condition", y = "counts", color = "condition", add = "jitter",palette = c("blue","orange","purple","red"))  +  facet_wrap(. ~ Genotype,ncol = 5,scales = "fixed",drop = FALSE) +  rotate_x_text(90)  + ylab(" Normalized Counts") + stat_compare_means(comparisons = list(c("W+N+","W-N+"),c("W+N+","W+N-"),c("W+N+","W-N-")), label = "p.signif",method = "t.test") + theme_bw() +   ggtitle(title)+ theme(axis.text.x = element_text(angle=90))

gene = "AT1G12940"
title = na.omit(newqpcr[newqpcr$geneID== gene ,])$geneSymbol[1]

p5 <- ggbarplot(  na.omit(newqpcr[newqpcr$geneID== gene ,]), x = "Cond", y = "meanExpression",  color = "Cond",palette = c("blue","orange","purple","red"),   lab.vjust = -1.6) +  facet_wrap(. ~ Geno  ,ncol = 5,scales = "fixed",drop = FALSE) +theme_bw()  +   geom_errorbar( aes(x=Cond, ymin=meanExpression-meanSEM, ymax=meanExpression + meanSEM), width=0.4, colour="black", alpha=0.9, linewidth=0.3) + theme_bw() + ylab("meanExpression ± meanSEM") + ggtitle(title) + theme(axis.text.x = element_text(angle=90))
p6 <-  ggboxplot(ymelt[ymelt$geneID==gene ,], x = "condition", y = "counts", color = "condition", add = "jitter",palette = c("blue","orange","purple","red"))  +  facet_wrap(. ~ Genotype,ncol = 5,scales = "fixed",drop = FALSE) +  rotate_x_text(90)  + ylab(" Normalized Counts") + stat_compare_means(comparisons = list(c("W+N+","W-N+"),c("W+N+","W+N-"),c("W+N+","W-N-")), label = "p.signif",method = "t.test") + theme_bw() +   ggtitle(title)+ theme(axis.text.x = element_text(angle=90))


gene = "AT1G25550"
title = na.omit(newqpcr[newqpcr$geneID== gene ,])$geneSymbol[1]

p7 <-  ggbarplot(  na.omit(newqpcr[newqpcr$geneID== gene ,]), x = "Cond", y = "meanExpression",  color = "Cond",palette = c("blue","orange","purple","red"),   lab.vjust = -1.6) +  facet_wrap(. ~ Geno  ,ncol = 5,scales = "fixed",drop = FALSE) +theme_bw()  +   geom_errorbar( aes(x=Cond, ymin=meanExpression-meanSEM, ymax=meanExpression + meanSEM), width=0.4, colour="black", alpha=0.9, linewidth=0.3) + theme_bw() +  ggtitle(title) + theme(axis.text.x = element_text(angle=90))
p8 <-   ggboxplot(ymelt[ymelt$geneID==gene ,], x = "condition", y = "counts", color = "condition", add = "jitter",palette = c("blue","orange","purple","red"))  +  facet_wrap(. ~ Genotype,ncol = 5,scales = "fixed",drop = FALSE) +  rotate_x_text(90)  + ylab(" Normalized Counts") + stat_compare_means(comparisons = list(c("W+N+","W-N+"),c("W+N+","W+N-"),c("W+N+","W-N-")), label = "p.signif",method = "t.test") + theme_bw() +   ggtitle(title)+ theme(axis.text.x = element_text(angle=90))

gene = "AT1G74810"
title = na.omit(newqpcr[newqpcr$geneID== gene ,])$geneSymbol[1]

p9 <- ggbarplot(  na.omit(newqpcr[newqpcr$geneID== gene ,]), x = "Cond", y = "meanExpression",  color = "Cond",palette = c("blue","orange","purple","red"),   lab.vjust = -1.6) +  facet_wrap(. ~ Geno  ,ncol = 5,scales = "fixed",drop = FALSE) +theme_bw()  +   geom_errorbar( aes(x=Cond, ymin=meanExpression-meanSEM, ymax=meanExpression + meanSEM), width=0.4, colour="black", alpha=0.9, linewidth=0.3) + theme_bw() + ylab("meanExpression ± meanSEM") + ggtitle(title) + theme(axis.text.x = element_text(angle=90))
p10 <-   ggboxplot(ymelt[ymelt$geneID==gene ,], x = "condition", y = "counts", color = "condition", add = "jitter",palette = c("blue","orange","purple","red"))  +  facet_wrap(. ~ Genotype,ncol = 5,scales = "fixed",drop = FALSE) +  rotate_x_text(90)  + ylab(" Normalized Counts") + stat_compare_means(comparisons = list(c("W+N+","W-N+"),c("W+N+","W+N-"),c("W+N+","W-N-")), label = "p.signif",method = "t.test") + theme_bw() +  ggtitle(title)+ theme(axis.text.x = element_text(angle=90))

gene = "AT2G47880"
title = na.omit(newqpcr[newqpcr$geneID== gene ,])$geneSymbol[1]

p11 <- ggbarplot(  na.omit(newqpcr[newqpcr$geneID== gene ,]), x = "Cond", y = "meanExpression",  color = "Cond",palette = c("blue","orange","purple","red"),   lab.vjust = -1.6) +  facet_wrap(. ~ Geno  ,ncol = 5,scales = "fixed",drop = FALSE) +theme_bw()  +   geom_errorbar( aes(x=Cond, ymin=meanExpression-meanSEM, ymax=meanExpression + meanSEM), width=0.4, colour="black", alpha=0.9, linewidth=0.3) + theme_bw() + ggtitle(title) + theme(axis.text.x = element_text(angle=90))
p12 <-  ggboxplot(ymelt[ymelt$geneID==gene ,], x = "condition", y = "counts", color = "condition", add = "jitter",palette = c("blue","orange","purple","red"))  +  facet_wrap(. ~ Genotype,ncol = 5,scales = "fixed",drop = FALSE) +  rotate_x_text(90)  + ylab(" Normalized Counts") + stat_compare_means(comparisons = list(c("W+N+","W-N+"),c("W+N+","W+N-"),c("W+N+","W-N-")), label = "p.signif",method = "t.test") + theme_bw() +   ggtitle(title)+ theme(axis.text.x = element_text(angle=90))

gene = "AT5G22300"
title = na.omit(newqpcr[newqpcr$geneID== gene ,])$geneSymbol[1]


p13 <- ggbarplot(  na.omit(newqpcr[newqpcr$geneID== gene ,]), x = "Cond", y = "meanExpression",  color = "Cond",palette = c("blue","orange","purple","red"),   lab.vjust = -1.6) +  facet_wrap(. ~ Geno  ,ncol = 5,scales = "fixed",drop = FALSE) +theme_bw()  +   geom_errorbar( aes(x=Cond, ymin=meanExpression-meanSEM, ymax=meanExpression + meanSEM), width=0.4, colour="black", alpha=0.9, linewidth=0.3) + theme_bw()  + ggtitle(title) + theme(axis.text.x = element_text(angle=90))
p14 <-   ggboxplot(ymelt[ymelt$geneID==gene ,], x = "condition", y = "counts", color = "condition", add = "jitter",palette = c("blue","orange","purple","red"))  +  facet_wrap(. ~ Genotype,ncol = 5,scales = "fixed",drop = FALSE) +  rotate_x_text(90)  + ylab(" Normalized Counts") + stat_compare_means(comparisons = list(c("W+N+","W-N+"),c("W+N+","W+N-"),c("W+N+","W-N-")), label = "p.signif",method = "t.test") + theme_bw()   + ggtitle(title)+ theme(axis.text.x = element_text(angle=90))


gene = "AT3G62960"
title = na.omit(newqpcr[newqpcr$geneID== gene ,])$geneSymbol[1]

p15 <- ggbarplot(  na.omit(newqpcr[newqpcr$geneID== gene ,]), x = "Cond", y = "meanExpression",  color = "Cond",palette = c("blue","orange","purple","red"),   lab.vjust = -1.6) +  facet_wrap(. ~ Geno  ,ncol = 5,scales = "fixed",drop = FALSE) +theme_bw()  +   geom_errorbar( aes(x=Cond, ymin=meanExpression-meanSEM, ymax=meanExpression + meanSEM), width=0.4, colour="black", alpha=0.9, linewidth=0.3) + theme_bw() + ylab("meanExpression ± meanSEM") + ggtitle(title) + theme(axis.text.x = element_text(angle=90))
p16 <-   ggboxplot(ymelt[ymelt$geneID==gene ,], x = "condition", y = "counts", color = "condition", add = "jitter",palette = c("blue","orange","purple","red"))  +  facet_wrap(. ~ Genotype,ncol = 5,scales = "fixed",drop = FALSE) +  rotate_x_text(90)  + ylab(" Normalized Counts") + stat_compare_means(comparisons = list(c("W+N+","W-N+"),c("W+N+","W+N-"),c("W+N+","W-N-")), label = "p.signif",method = "t.test") + theme_bw() + ylab("counts ± meanSEM")  + ggtitle(title)+ theme(axis.text.x = element_text(angle=90))


gene = "AT4G15680"
title = na.omit(newqpcr[newqpcr$geneID== gene ,])$geneSymbol[1]

p17 <-  ggbarplot(  na.omit(newqpcr[newqpcr$geneID== gene ,]), x = "Cond", y = "meanExpression",  color = "Cond",palette = c("blue","orange","purple","red"),   lab.vjust = -1.6) +  facet_wrap(. ~ Geno  ,ncol = 5,scales = "fixed",drop = FALSE) +theme_bw()  +   geom_errorbar( aes(x=Cond, ymin=meanExpression-meanSEM, ymax=meanExpression + meanSEM), width=0.4, colour="black", alpha=0.9, linewidth=0.3) + theme_bw() + ylab("meanExpression ± meanSEM") + ggtitle(title) + theme(axis.text.x = element_text(angle=90))
p18 <-  ggboxplot(ymelt[ymelt$geneID==gene ,], x = "condition", y = "counts", color = "condition", add = "jitter",palette = c("blue","orange","purple","red"))  +  facet_wrap(. ~ Genotype,ncol = 5,scales = "fixed",drop = FALSE) +  rotate_x_text(90)  + ylab(" Normalized Counts") + stat_compare_means(comparisons = list(c("W+N+","W-N+"),c("W+N+","W+N-"),c("W+N+","W-N-")), label = "p.signif",method = "t.test") + theme_bw() + ylab("counts ± meanSEM")  + ggtitle(title)+ theme(axis.text.x = element_text(angle=90))

gene = "AT2G40330"
title = na.omit(newqpcr[newqpcr$geneID== gene ,])$geneSymbol[1]

p19 <- ggbarplot(  na.omit(newqpcr[newqpcr$geneID== gene ,]), x = "Cond", y = "meanExpression",  color = "Cond",palette = c("blue","orange","purple","red"),   lab.vjust = -1.6) +  facet_wrap(. ~ Geno  ,ncol = 5,scales = "fixed",drop = FALSE) +theme_bw()  +   geom_errorbar( aes(x=Cond, ymin=meanExpression-meanSEM, ymax=meanExpression + meanSEM), width=0.4, colour="black", alpha=0.9, linewidth=0.3) + theme_bw() + ylab("meanExpression ± meanSEM") + ggtitle(title) + theme(axis.text.x = element_text(angle=90))
p20 <-  ggboxplot(ymelt[ymelt$geneID==gene ,], x = "condition", y = "counts", color = "condition", add = "jitter",palette = c("blue","orange","purple","red"))  +  facet_wrap(. ~ Genotype,ncol = 5,scales = "fixed",drop = FALSE) +  rotate_x_text(90)  + ylab(" Normalized Counts") + stat_compare_means(comparisons = list(c("W+N+","W-N+"),c("W+N+","W+N-"),c("W+N+","W-N-")), label = "p.signif",method = "t.test") + theme_bw() + ylab("counts ± meanSEM")  + ggtitle(title)+ theme(axis.text.x = element_text(angle=90))

gene = "AT3G62950"
title = na.omit(newqpcr[newqpcr$geneID== gene ,])$geneSymbol[1]


p21 <- ggbarplot(  na.omit(newqpcr[newqpcr$geneID== gene ,]), x = "Cond", y = "meanExpression",  color = "Cond",palette = c("blue","orange","purple","red"),   lab.vjust = -1.6) +  facet_wrap(. ~ Geno  ,ncol = 5,scales = "fixed",drop = FALSE) +theme_bw()  +   geom_errorbar( aes(x=Cond, ymin=meanExpression-meanSEM, ymax=meanExpression + meanSEM), width=0.4, colour="black", alpha=0.9, linewidth=0.3) + theme_bw() + ylab("meanExpression ± meanSEM") + ggtitle(title) + theme(axis.text.x = element_text(angle=90))
p22 <-   ggboxplot(ymelt[ymelt$geneID==gene ,], x = "condition", y = "counts", color = "condition", add = "jitter",palette = c("blue","orange","purple","red"))  +  facet_wrap(. ~ Genotype,ncol = 5,scales = "fixed",drop = FALSE) +  rotate_x_text(90)  + ylab(" Normalized Counts") + stat_compare_means(comparisons = list(c("W+N+","W-N+"),c("W+N+","W+N-"),c("W+N+","W-N-")), label = "p.signif",method = "t.test") + theme_bw() + ylab("counts ± meanSEM")  + ggtitle(title)+ theme(axis.text.x = element_text(angle=90))
# 

gene = "AT5G19470"
title = na.omit(newqpcr[newqpcr$geneID== gene ,])$geneSymbol[1]


p23 <- ggbarplot(  na.omit(newqpcr[newqpcr$geneID== gene ,]), x = "Cond", y = "meanExpression",  color = "Cond",palette = c("blue","orange","purple","red"),   lab.vjust = -1.6) +  facet_wrap(. ~ Geno  ,ncol = 5,scales = "fixed",drop = FALSE) +theme_bw()  +   geom_errorbar( aes(x=Cond, ymin=meanExpression-meanSEM, ymax=meanExpression + meanSEM), width=0.4, colour="black", alpha=0.9, linewidth=0.3) + theme_bw() + ylab("meanExpression ± meanSEM") + ggtitle(title) + theme(axis.text.x = element_text(angle=90))
p24 <-  ggboxplot(ymelt[ymelt$geneID==gene ,], x = "condition", y = "counts", color = "condition", add = "jitter",palette = c("blue","orange","purple","red"))  +  facet_wrap(. ~ Genotype,ncol = 5,scales = "fixed",drop = FALSE) +  rotate_x_text(90)  + ylab(" Normalized Counts") + stat_compare_means(comparisons = list(c("W+N+","W-N+"),c("W+N+","W+N-"),c("W+N+","W-N-")), label = "p.signif",method = "t.test") + theme_bw() + ylab("counts ± meanSEM")  + ggtitle(title)+ theme(axis.text.x = element_text(angle=90))

gene = "AT1G34060"
title = na.omit(newqpcr[newqpcr$geneID== gene ,])$geneSymbol[1]


p25 <- ggbarplot(  na.omit(newqpcr[newqpcr$geneID== gene ,]), x = "Cond", y = "meanExpression",  color = "Cond",palette = c("blue","orange","purple","red"),   lab.vjust = -1.6) +  facet_wrap(. ~ Geno  ,ncol = 5,scales = "fixed",drop = FALSE) +theme_bw()  +   geom_errorbar( aes(x=Cond, ymin=meanExpression-meanSEM, ymax=meanExpression + meanSEM), width=0.4, colour="black", alpha=0.9, linewidth=0.3) + theme_bw() + ylab("meanExpression ± meanSEM") + ggtitle(title) + theme(axis.text.x = element_text(angle=90))
p26 <-  ggboxplot(ymelt[ymelt$geneID==gene ,], x = "condition", y = "counts", color = "condition", add = "jitter",palette = c("blue","orange","purple","red"))  +  facet_wrap(. ~ Genotype,ncol = 5,scales = "fixed",drop = FALSE) +  rotate_x_text(90)  + ylab(" Normalized Counts") + stat_compare_means(comparisons = list(c("W+N+","W-N+"),c("W+N+","W+N-"),c("W+N+","W-N-")), label = "p.signif",method = "t.test") + theme_bw() + ylab("counts ± meanSEM")  + ggtitle(title)+ theme(axis.text.x = element_text(angle=90))

gene = "AT5G14070"
title = na.omit(newqpcr[newqpcr$geneID== gene ,])$geneSymbol[1]

p27 <- ggbarplot(  na.omit(newqpcr[newqpcr$geneID== gene ,]), x = "Cond", y = "meanExpression",  color = "Cond",palette = c("blue","orange","purple","red"),   lab.vjust = -1.6) +  facet_wrap(. ~ Geno  ,ncol = 5,scales = "fixed",drop = FALSE) +theme_bw()  +   geom_errorbar( aes(x=Cond, ymin=meanExpression-meanSEM, ymax=meanExpression + meanSEM), width=0.4, colour="black", alpha=0.9, linewidth=0.3) + theme_bw() + ylab("meanExpression ± meanSEM") + ggtitle(title) + theme(axis.text.x = element_text(angle=90))
p28 <-  ggboxplot(ymelt[ymelt$geneID==gene ,], x = "condition", y = "counts", color = "condition", add = "jitter",palette = c("blue","orange","purple","red"))  +  facet_wrap(. ~ Genotype,ncol = 5,scales = "fixed",drop = FALSE) +  rotate_x_text(90)  + ylab(" Normalized Counts") + stat_compare_means(comparisons = list(c("W+N+","W-N+"),c("W+N+","W+N-"),c("W+N+","W-N-")), label = "p.signif",method = "t.test") + theme_bw() + ylab("counts ± meanSEM")  + ggtitle(title)+ theme(axis.text.x = element_text(angle=90))


(p1|p2)/(p3|p4)/(p5|p6)/(p7|p8)/(p9|p10) 
ggsave("figures/Figure_S3-1.jpg", width = 30, height = 60, units = c("cm"), dpi = 600)
ggsave("figures/Figure_S3-1.svg", width = 30, height = 60, units = c("cm"), dpi = 600)

(p11|p12)/(p13|p14)/ (p15|p16)/(p17|p18)/ (p19|p20) 
ggsave("figures/Figure_S3-2.jpg", width = 30, height = 60, units = c("cm"), dpi = 600)
ggsave("figures/Figure_S3-2.svg", width = 30, height = 60, units = c("cm"), dpi = 600)

(p21|p22)/(p23|p24)/(p25|p26)/(p27|p28) 
ggsave("figures/Figure_S3-3.jpg", width = 30, height = 60, units = c("cm"), dpi = 600)
ggsave("figures/Figure_S3-3.svg", width = 30, height = 60, units = c("cm"), dpi = 600)
