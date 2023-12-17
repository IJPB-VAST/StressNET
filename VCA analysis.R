### Anova 

counts<- read.table(paste0("data/Counts/merge_counts_scaled_log2.tsv"),check.names = F)

adjvalue = 0.05
alpha = 1
DEG=NULL

for (ecotype in c("Col-0","Sha","Cvi-0","Tsu-0","Bur-0"))  {
  for  (cond in c("W-N+","W+N-","W-N-"))  {
    res <- read.table(paste0("data/Counts/",ecotype,'_',cond,"_refer_whole.tsv"),  header=T)
    res<- subset(res,padj< adjvalue & abs(log2FoldChange) > alpha)
    DEG<- c(DEG,as.character(res$row))
    DEG<- DEG[!duplicated(DEG)] 
  }
}


data <- reshape::melt(counts)
data<- separate(data, "variable", c("G","Cond","Exp","rep"), sep =  "_",remove=FALSE)
data$W <- gsub("N[0-9]","",data$Cond)
data$N <- gsub("W[0-9]","",data$Cond)

mlist<- NULL
for (i in levels(factor(DEG)) ){
  subdata <- na.omit(data[data$GeneID == i,])
  if ( length(levels(factor(subdata$G))) > 1 ){
    pW <- t(anovaVCA(value ~ G*W*N*Exp, subdata)$aov.tab[-1,5])
    q <-  data.frame(i,pW)
    colnames(q)<- c("name",colnames(q[,-1]))
    mlist<- rbind(mlist,q)
    colnames(mlist)<- colnames(q)
  }
}

mlist$residual <- mlist$G.W + mlist$G.N + mlist$G.Exp + mlist$W.Exp + mlist$N.Exp+ mlist$G.W.Exp + mlist$G.N.Exp + mlist$W.N.Exp + mlist$G.W.N.Exp + mlist$error

mlist <- mlist[,(colnames(mlist) %in% c("name","G","W","N","W.N","G.W.N", "Exp","residual"))]

write.csv(mlist,'data/vca.csv')
write.csv(DEG,'data/deg_genelist.csv')



ST002201_AN003604 <- read.delim("data/ST002201_AN003604.txt", header=FALSE)
metabolite_matrix<- ST002201_AN003604[2:137,]
label <- ST002201_AN003604[139:nrow(ST002201_AN003604),]
colnames(metabolite_matrix) <- ST002201_AN003604[1,]
colnames(label) <-  ST002201_AN003604[138,]
label$class <- gsub("^.*_","",label$metabolite_fullname)

data <- reshape::melt(metabolite_matrix,id = "metabolite_name")
data <- data %>% separate(variable,c("Geno","Cond","Exp","rep"),sep = "_")
Meta= levels(as.factor(data$metabolite_name))
Genotype=levels(as.factor(data$Geno))
data$value <- as.numeric(data$value)
Pval= NULL

for (j in Meta ){
  for (i in Genotype ){
    Metaj <- subset(data, metabolite_name == j & data$Geno ==i)
    PvaljW <- summary(aov(lm(value ~ Cond, data = Metaj[Metaj$Cond %in% c("W+N+","W-N+"),])))[[1]][1,5]
    PvaljN <- summary(aov(lm(value ~ Cond, data = Metaj[Metaj$Cond %in% c("W+N+","W+N-"),])))[[1]][1,5]
    PvaljNW <- summary(aov(lm(value ~ Cond, data = Metaj[Metaj$Cond %in% c("W+N+","W-N-"),])))[[1]][1,5]
    Pvalj <- cbind(PvaljW,PvaljN,PvaljNW)
    colnames(Pvalj) <- c("Wstress","Nstress","Combined_stress")
    Labeli <- data.frame(j,i,Pvalj)
    Pval <- rbind(Pval,Labeli)
  }
} 

names(Pval) <- c("metabolite","Geno","Wstress","Nstress","combined_stress")
Pval$Geno.padj<- p.adjust(Pval$Geno, method = "BH")
Pval$Wstress.padj<- p.adjust(Pval$Wstress, method = "BH")
Pval$Nstress.padj<- p.adjust(Pval$Nstress, method = "BH")
Pval$combined_stress.padj<- p.adjust(Pval$combined_stress, method = "BH")
write.csv(Pval,'data/Pval_metabolite.csv')



### metabolite vca
ST002201_AN003604 <- read.delim("data/ST002201_AN003604.txt", header=FALSE)
metabolite_matrix<- ST002201_AN003604[2:137,]
label <- ST002201_AN003604[139:nrow(ST002201_AN003604),]
colnames(metabolite_matrix) <- ST002201_AN003604[1,]
colnames(label) <-  ST002201_AN003604[138,]
label$class <- gsub("^.*_","",label$metabolite_fullname)

data <- reshape::melt(metabolite_matrix,id = "metabolite_name")
data <- data %>% separate(variable,c("G","Cond","Exp","rep"),sep = "_")
data$Cond <- gsub('\\+',1,data$Cond)
data$Cond <- gsub('\\-',0,data$Cond)
data$W <- gsub("N[0-9]","",data$Cond)
data$N <- gsub("W[0-9]","",data$Cond)
data$value <- as.numeric(data$value)

mlist<- NULL
for (i in levels(factor(data$metabolite_name)) ){
  subdata <- na.omit(data[data$metabolite_name == i,])
  if ( length(levels(factor(subdata$G))) > 1 ){
    pW <- t(anovaVCA(value  ~ G*W*N*Exp, subdata)$aov.tab[-1,5])
    q <-  data.frame(i,pW)
    colnames(q)<- c("name",colnames(q[,-1]))
    mlist<- rbind(mlist,q)
    colnames(mlist)<- colnames(q)
  }
}

mlist$residual <- mlist$G.W + mlist$G.N + mlist$G.Exp + mlist$W.Exp + mlist$N.Exp+ mlist$G.W.Exp + mlist$G.N.Exp + mlist$W.N.Exp + mlist$G.W.N.Exp + mlist$error

mlist <- mlist[,(colnames(mlist) %in% c("name","G","W","N","W.N","G.W.N", "Exp","residual"))]

write.csv(mlist,'data/vca_metabo.csv')

