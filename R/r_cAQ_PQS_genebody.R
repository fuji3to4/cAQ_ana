library(patchwork)
library(ggridges)
#source("r_cAQ_PQS_ChIPseeker.R")



#genebody-----------------------------------------------------
data<-peakAnnoList[["gencode.v29"]]@anno %>% as.data.frame()
ids.genebody <- grep("Promoter",data$annotation,invert = T)
ids.notna<-which(!is.na(data$annotation))
ids.genebody<-intersect(ids.genebody,ids.notna)
#data[ids.genebody,] %>% head(.,100)

g4.count.genebody<-data[ids.genebody,"geneId"] %>% table

gencode.genes<-genes(gencode.v29.txdb)$gene_id
notin.genes<-gencode.genes[!(gencode.genes %in% names(g4.count.genebody))]
x<-rep(0,length(notin.genes))
names(x)<-notin.genes
g4.count.genebody<-c(g4.count.genebody,x)


#ttt<-ChIPseeker:::getGene(gencode.v29.txdb,by="gene")
gene<-genes(gencode.v29.txdb)
genelength<-width(gene)
names(genelength)<-names(gene)
genelength<-as.data.frame(genelength)

g4.count.gb.genelen<-genelength[names(g4.count.genebody),]

g4.count.genebody.df.list<-list("all"=g4.count.genebody,"Width"=g4.count.gb.genelen)
g4.count.genebody.df<-data.frame(Category="all",Gene=names(g4.count.genebody),Count=as.numeric(g4.count.genebody),Width=as.numeric(g4.count.gb.genelen))
for(name in names(toptags.genes.ensembl)){
  ids<-names(g4.count.genebody) %in% toptags.genes.ensembl[[name]]
  tmp.count<-g4.count.genebody[ids]
  tmp.genelen<-g4.count.gb.genelen[ids]
  g4.count.genebody.df.list[[name]]<-tmp.count
  tmp.df<-data.frame(Category=name,Gene=names(tmp.count),Count=as.numeric(tmp.count),Width=as.numeric(tmp.genelen))
  g4.count.genebody.df<-rbind(g4.count.genebody.df,tmp.df)
}
boxplot(log10(Count+1)~Category,data=g4.count.genebody.df)
boxplot(Count/Width~Category,data=g4.count.genebody.df)

group_by(g4.count.genebody.df,Category) %>% summarise(median=median(Count),mean=mean(Count))
group_by(g4.count.genebody.df,Category) %>% summarise(median=median(Count/Width),mean=mean(Count/Width))

plot(Count~Width,data=g4.count.genebody.df,cex=.1)


p.riges <- ggplot(g4.count.genebody.df,aes(x=log10(Count+1),y=Category,fill=Category)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2, alpha = .7,from=0,
                      vline_color = "red",vline_size = 0.5)+
  scale_x_continuous(name="log10(PQS Counts)") +
  scale_y_discrete(limits=c("all",names(toptags.genes.ensembl))) +
  scale_fill_cyclical(values =c("gray",rep(c("#8080ff","#ff8080"),4)))+
  theme_bw()

p.riges
#ggsave(file = "./out/geom_density_PQScounts_genebody.pdf", plot = p.riges, width = 11, height = 8.5)

p.riges <- ggplot(g4.count.genebody.df,aes(x=Count/Width,y=Category,fill=Category)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2, vline_color = "red",alpha = .7,from=0,to=0.002)+
  scale_x_continuous(name="PQS Counts / Gene Length") +
  scale_y_discrete(limits=c("all",names(toptags.genes.ensembl))) +
  scale_fill_cyclical(values =c("gray",rep(c("#8080ff","#ff8080"),4)))+
  theme_bw()

p.riges
#ggsave(file = "./out/geom_density_PQScounts_div_genelen_genebody.pdf", plot = p.riges, width = 11, height = 8.5)




sapply(g4.count.genebody.df.list,length)
g4.count.genebody.dummy<-data.frame(gene=names(g4.count.genebody.df.list[["all"]]),Count=as.numeric(g4.count.genebody.df.list[["all"]]),Width=as.numeric(g4.count.gb.genelen))
for(name in names(toptags.genes.ensembl)){
  g4.count.genebody.dummy[,name]<-g4.count.genebody.dummy$gene %in% names(g4.count.genebody.df.list[[name]]) %>% as.integer() %>% as.factor() 
}



glm.genebody.g4<-list()
for(name in names(toptags.genes.ensembl)){
  glm.genebody.g4[[name]] <- glm(formula = as.formula(paste0("Count~",name)),data=g4.count.genebody.dummy,family = poisson())
}
lapply(glm.genebody.g4, summary)

glm.offset.genebody.g4<-list()
for(name in names(toptags.genes.ensembl)){
  glm.offset.genebody.g4[[name]] <- glm(formula = as.formula(paste0("Count~",name)),offset = log(Width),data=g4.count.genebody.dummy,family = poisson())
}
lapply(glm.offset.genebody.g4, summary)


wilcox.genebody.g4<-list()
for(name in names(toptags.genes.ensembl)){
  wilcox.genebody.g4[[name]] <- wilcox.test(formula = as.formula(paste0("Count~",name)),data=g4.count.genebody.dummy,family = poisson())
}
wilcox.genebody.g4

wilcox.genebody.offset.g4<-list()
for(name in names(toptags.genes.ensembl)){
  wilcox.genebody.offset.g4[[name]] <- wilcox.test(formula = as.formula(paste0("Count/Width~",name)),data=g4.count.genebody.dummy,family = poisson())
}
wilcox.genebody.offset.g4

g4.count.genebody.reg<-data.frame(Category=NULL,gene=NULL,Count=NULL,Width=NULL,reg=NULL)
for(name in names(toptags.genes.ensembl)){
  tmp<-data.frame(Category=name,gene=g4.count.genebody.dummy$gene,Count=g4.count.genebody.dummy$Count,Width=g4.count.genebody.dummy$Width,reg=g4.count.genebody.dummy[,name])
  g4.count.genebody.reg<-rbind(g4.count.genebody.reg,tmp)
}

#write.csv(g4.count.genebody.reg,file="./out/g4.count.genebody.reg.csv",row.names = F,quote = F)


g4.count.genebody.reg$Category<-as.factor(g4.count.genebody.reg$Category)

#glm.offset.genebody.sum<-glm(Count~Category+reg,offset = log(Width),data=g4.count.genebody.reg,family = poisson())
#summary(glm.offset.genebody.sum)

#Intron--------------------------------------
data<-peakAnnoList[["gencode.v29"]]@anno %>% as.data.frame()
#peakAnnoList[["gencode.v29"]]@annoStat$Feature
ids.intron<-grep("intron",data$annotation)
#ids.1stintron<-grep("intron 1 of",data$annotation)
#ids.otherintron<-setdiff(ids.intron,ids.1stintron)

#ChPseeker:::getGenomicAnnotation.internal
genomicRegion<-intronList <- intronsByTranscript(gencode.v29.txdb)
GRegion <- unlist(genomicRegion)
GRegionLen <- elementNROWS(genomicRegion)
names(GRegionLen) <- names(genomicRegion)
GRegion$gene_id <- rep(names(genomicRegion), times = GRegionLen)
gr2 <- GRegion[!duplicated(GRegion$gene_id)]
strd <- as.character(strand(gr2))
len <- GRegionLen[GRegionLen != 0]
GRegion$intron_rank <- lapply(seq_along(strd), function(i) {
  rank <- seq(1, len[i])
  if (strd[i] == "-") 
    rank <- rev(rank)
  return(rank)
}) %>% unlist

txidinfo <- transcripts(gencode.v29.txdb, columns = c("tx_id", "tx_name", "gene_id"))
idx <- which(sapply(txidinfo$gene_id, length) == 0)
txidinfo[idx, ]$gene_id <- txidinfo[idx, ]$tx_name
txid2geneid <- paste(mcols(txidinfo)[["tx_name"]], mcols(txidinfo)[["gene_id"]], sep = "/")
txid2geneid <- sub("/NA", "", txid2geneid)
names(txid2geneid) <- mcols(txidinfo)[["tx_id"]]

nn<-txid2geneid[names(genomicRegion)]
#nn <- ChIPseeker:::TXID2EG(names(genomicRegion))
names(GRegionLen) <- nn
GRegion$gene_id <- rep(nn, times = GRegionLen)
anno <- paste("Intron", " (", GRegion$gene_id, ", intron ", GRegion$intron_rank, " of ", GRegionLen[GRegion$gene_id], ")", sep = "")
g4.intron<-data[ids.intron,]
which(anno %in% g4.intron$annotation[10])

intron.width<-width(GRegion)
names(intron.width)<-anno


g4.count.each.intron<-g4.intron$annotation %>% table %>% as.data.frame()

notin.genes<-anno[!(anno %in% g4.count.each.intron[,1])]
tmp<-data.frame("."=notin.genes,"Freq"=rep(0,length(notin.genes)))
g4.count.each.intron<-rbind(g4.count.each.intron,tmp)

row.names(g4.count.each.intron)<-g4.count.each.intron[,1]
g4.count.each.intron[,1]<-NULL

g4.count.each.intron$width<-intron.width[row.names(g4.count.each.intron)]

g4.count.each.intron$gene_id<-sub("^Intron \\(ENST[0-9]+\\.[0-9]+/(ENSG[0-9]+\\.[0-9]+),.+","\\1",row.names(g4.count.each.intron),perl = T)
g4.count.each.intron$tx_name<-sub("^Intron \\((ENST[0-9]+\\.[0-9]+)/ENSG[0-9]+\\.[0-9]+,.+","\\1",row.names(g4.count.each.intron),perl = T)
g4.count.each.intron<-g4.count.each.intron[,c(3,4,1,2)]
colnames(g4.count.each.intron)[3]<-"Count"

g4.count.intron.dummy<-group_by(g4.count.each.intron,gene_id) %>% summarise(Count=sum(Count),Width=sum(width))

for(name in names(toptags.genes.ensembl)){
  g4.count.intron.dummy[,name]<-g4.count.intron.dummy$gene_id %in% toptags.genes.ensembl[[name]] %>% as.integer() %>% as.factor()
}
g4.count.intron.dummy

g4.count.intron.df<-g4.count.intron.dummy[,1:3]
g4.count.intron.df$Category <-"all"
for(name in names(toptags.genes.ensembl)){
  tmp<-g4.count.intron.dummy[which(g4.count.intron.dummy[,name]==1),1:3]
  tmp$Category <-name
  g4.count.intron.df<-rbind(g4.count.intron.df,tmp)
}
boxplot(log10(Count+1)~Category,data=g4.count.intron.df)
boxplot(Count/Width~Category,data=g4.count.intron.df,ylim=c(0,0.002))

group_by(g4.count.intron.df,Category) %>% summarise(median=median(Count),mean=mean(Count))
group_by(g4.count.intron.df,Category) %>% summarise(median=median(Count/Width),mean=mean(Count/Width))


library(ggridges)
p.riges <- ggplot(g4.count.intron.df,aes(x=log10(Count+1),y=Category,fill=Category)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2, vline_color = "red",alpha = .7,from=0)+
  scale_x_continuous(name="log10(PQS Counts)") +
  scale_y_discrete(limits=c("all",names(toptags.genes.ensembl))) +
  scale_fill_cyclical(values =c("gray",rep(c("#8080ff","#ff8080"),4)))+
  theme_bw()

p.riges
ggsave(file = "./out/geom_density_PQScounts_intron.pdf", plot = p.riges, width = 11, height = 8.5)

p.riges <- ggplot(g4.count.intron.df,aes(x=Count/Width,y=Category,fill=Category)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2, vline_color = "red",alpha = .7,from=0,to=0.002)+
  scale_x_continuous(name="PQS Counts/Width") +
  scale_y_discrete(limits=c("all",names(toptags.genes.ensembl))) +
  scale_fill_cyclical(values =c("gray",rep(c("#8080ff","#ff8080"),4)))+
  theme_bw()

p.riges
ggsave(file = "./out/geom_density_PQScounts_div_genelen_intron.pdf", plot = p.riges, width = 11, height = 8.5)



glm.intron.g4<-list()
for(name in names(toptags.genes.ensembl)){
  glm.intron.g4[[name]] <- glm(formula = as.formula(paste0("Count~",name)),data=g4.count.intron.dummy,family = poisson())
}
lapply(glm.intron.g4, summary)

glm.offset.intron.g4<-list()
for(name in names(toptags.genes.ensembl)){
  glm.offset.intron.g4[[name]] <- glm(formula = as.formula(paste0("Count~",name)),offset = log(Width),data=g4.count.intron.dummy,family = poisson())
}
lapply(glm.offset.intron.g4, summary)

g4.count.intron.reg<-data.frame(Category=NULL,gene=NULL,Count=NULL,Width=NULL,reg=NULL)
for(name in names(toptags.genes.ensembl)){
  tmp<-data.frame(Category=name,gene=g4.count.intron.dummy$gene_id,Count=g4.count.intron.dummy$Count,Width=g4.count.intron.dummy$Width,reg=g4.count.intron.dummy[,name])
  colnames(tmp)[5]<-"reg"
  g4.count.intron.reg<-rbind(g4.count.intron.reg,tmp)
}

#write.csv(g4.count.intron.reg,file="./out/g4.count.intron.reg.csv",row.names = F,quote = F)


#1st intron
ids.1stintron<-grep("intron 1 of",row.names(g4.count.each.intron))
ids.otherintron<-grep("intron 1 of",row.names(g4.count.each.intron),invert = T)
g4.count.each.intron[ids.1stintron,]

g4.count.1stintron.dummy<-group_by(g4.count.each.intron[ids.1stintron,],gene_id) %>% summarise(Count=sum(Count),Width=sum(width))

for(name in names(toptags.genes.ensembl)){
  g4.count.1stintron.dummy[,name]<-g4.count.1stintron.dummy$gene_id %in% toptags.genes.ensembl[[name]] %>% as.integer() %>% as.factor()
}
g4.count.1stintron.dummy

g4.count.1stintron.df<-g4.count.1stintron.dummy[,1:3]
g4.count.1stintron.df$Category <-"all"
for(name in names(toptags.genes.ensembl)){
  tmp<-g4.count.1stintron.dummy[which(g4.count.1stintron.dummy[,name]==1),1:3]
  tmp$Category <-name
  g4.count.1stintron.df<-rbind(g4.count.1stintron.df,tmp)
}

group_by(g4.count.1stintron.df,Category) %>% summarise(median=median(Count),mean=mean(Count))
group_by(g4.count.1stintron.df,Category) %>% summarise(median=median(Count/Width),mean=mean(Count/Width))

boxplot(log10(Count+1)~Category,data=g4.count.1stintron.df)
boxplot(Count/Width~Category,data=g4.count.1stintron.df,ylim=c(0,0.002))

library(ggridges)
p.riges <- ggplot(g4.count.1stintron.df,aes(x=log10(Count+1),y=Category,fill=Category)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2, vline_color = "red",alpha = .7,from=0)+
  scale_x_continuous(name="log10(PQS Counts)") +
  scale_y_discrete(limits=c("all",names(toptags.genes.ensembl))) +
  scale_fill_cyclical(values =c("gray",rep(c("#8080ff","#ff8080"),4)))+
  theme_bw()

p.riges
ggsave(file = "./out/geom_density_PQScounts_1stintron.pdf", plot = p.riges, width = 11, height = 8.5)

p.riges <- ggplot(g4.count.1stintron.df,aes(x=Count/Width,y=Category,fill=Category)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2, vline_color = "red",alpha = .7,from=0,to=0.001)+
  scale_x_continuous(name="PQS Counts/Width") +
  scale_y_discrete(limits=c("all",names(toptags.genes.ensembl))) +
  scale_fill_cyclical(values =c("gray",rep(c("#8080ff","#ff8080"),4)))+
  theme_bw()

p.riges
ggsave(file = "./out/geom_density_PQScounts_div_genelen_1stintron.pdf", plot = p.riges, width = 11, height = 8.5)

glm.1stintron.g4<-list()
for(name in names(toptags.genes.ensembl)){
  glm.1stintron.g4[[name]] <- glm(formula = as.formula(paste0("Count~",name)),data=g4.count.1stintron.dummy,family = poisson())
}
lapply(glm.1stintron.g4, summary)

glm.offset.1stintron.g4<-list()
for(name in names(toptags.genes.ensembl)){
  glm.offset.1stintron.g4[[name]] <- glm(formula = as.formula(paste0("Count~",name)),offset = log(Width),data=g4.count.1stintron.dummy,family = poisson())
}
lapply(glm.offset.1stintron.g4, summary)


g4.count.1stintron.reg<-data.frame(Category=NULL,gene=NULL,Count=NULL,Width=NULL,reg=NULL)
for(name in names(toptags.genes.ensembl)){
  tmp<-data.frame(Category=name,gene=g4.count.1stintron.dummy$gene_id,Count=g4.count.1stintron.dummy$Count,Width=g4.count.1stintron.dummy$Width,reg=g4.count.1stintron.dummy[,name])
  colnames(tmp)[5]<-"reg"
  g4.count.1stintron.reg<-rbind(g4.count.1stintron.reg,tmp)
}

#write.csv(g4.count.1stintron.reg,file="./out/g4.count.1stintron.reg.csv",row.names = F,quote = F)

#other intron
ids.otherintron<-grep("intron 1 of",row.names(g4.count.each.intron),invert = T)

g4.count.otherintron.dummy<-group_by(g4.count.each.intron[ids.otherintron,],gene_id) %>% summarise(Count=sum(Count),Width=sum(width))

for(name in names(toptags.genes.ensembl)){
  g4.count.otherintron.dummy[,name]<-g4.count.otherintron.dummy$gene_id %in% toptags.genes.ensembl[[name]] %>% as.integer() %>% as.factor()
}
g4.count.otherintron.dummy

g4.count.otherintron.df<-g4.count.otherintron.dummy[,1:3]
g4.count.otherintron.df$Category <-"all"
for(name in names(toptags.genes.ensembl)){
  tmp<-g4.count.otherintron.dummy[which(g4.count.otherintron.dummy[,name]==1),1:3]
  tmp$Category <-name
  g4.count.otherintron.df<-rbind(g4.count.otherintron.df,tmp)
}

group_by(g4.count.otherintron.df,Category) %>% summarise(median=median(Count),mean=mean(Count))
group_by(g4.count.otherintron.df,Category) %>% summarise(median=median(Count/Width),mean=mean(Count/Width))

boxplot(log10(Count+1)~Category,data=g4.count.otherintron.df)
boxplot(Count/Width~Category,data=g4.count.otherintron.df,ylim=c(0,0.002))

library(ggridges)
p.riges <- ggplot(g4.count.otherintron.df,aes(x=log10(Count+1),y=Category,fill=Category)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2, vline_color = "red",alpha = .7)+
  scale_x_continuous(name="log10(PQS Counts)") +
  scale_y_discrete(limits=c("all",names(toptags.genes.ensembl))) +
  scale_fill_cyclical(values =c("gray",rep(c("#8080ff","#ff8080"),4)))+
  theme_bw()

p.riges
ggsave(file = "./out/geom_density_PQScounts_otherintron.pdf", plot = p.riges, width = 11, height = 8.5)

p.riges <- ggplot(g4.count.otherintron.df,aes(x=Count/Width,y=Category,fill=Category)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2, vline_color = "red",alpha = .7,from=0,to=0.002)+
  scale_x_continuous(name="PQS Counts/Width") +
  scale_y_discrete(limits=c("all",names(toptags.genes.ensembl))) +
  scale_fill_cyclical(values =c("gray",rep(c("#8080ff","#ff8080"),4)))+
  theme_bw()

p.riges
ggsave(file = "./out/geom_density_PQScounts_div_genelen_otherintron.pdf", plot = p.riges, width = 11, height = 8.5)


glm.otherintron.g4<-list()
for(name in names(toptags.genes.ensembl)){
  glm.otherintron.g4[[name]] <- glm(formula = as.formula(paste0("Count~",name)),data=g4.count.otherintron.dummy,family = poisson())
}
lapply(glm.otherintron.g4, summary)

glm.offset.otherintron.g4<-list()
for(name in names(toptags.genes.ensembl)){
  glm.offset.otherintron.g4[[name]] <- glm(formula = as.formula(paste0("Count~",name)),offset = log(Width),data=g4.count.otherintron.dummy,family = poisson())
}
lapply(glm.offset.otherintron.g4, summary)

g4.count.otherintron.reg<-data.frame(Category=NULL,gene=NULL,Count=NULL,Width=NULL,reg=NULL)
for(name in names(toptags.genes.ensembl)){
  tmp<-data.frame(Category=name,gene=g4.count.otherintron.dummy$gene_id,Count=g4.count.otherintron.dummy$Count,Width=g4.count.otherintron.dummy$Width,reg=g4.count.otherintron.dummy[,name])
  colnames(tmp)[5]<-"reg"
  g4.count.otherintron.reg<-rbind(g4.count.otherintron.reg,tmp)
}

#write.csv(g4.count.otherintron.reg,file="./out/g4.count.otherintron.reg.csv",row.names = F,quote = F)

#Exon------------------------------------------------------------
data<-peakAnnoList[["gencode.v29"]]@anno %>% as.data.frame()
#peakAnnoList[["gencode.v29"]]@annoStat$Feature
ids.exon<-grep("exon",data$annotation)
#ids.1stexon<-grep("exon 1 of",data$annotation)
#ids.otherexon<-setdiff(ids.exon,ids.1stexon)

#ChPseeker:::getGenomicAnnotation.internal
genomicRegion<-exonList <- exonsBy(gencode.v29.txdb)
GRegion <- unlist(genomicRegion)
GRegionLen <- elementNROWS(genomicRegion)
names(GRegionLen) <- names(genomicRegion)
GRegion$gene_id <- rep(names(genomicRegion), times = GRegionLen)
gr2 <- GRegion[!duplicated(GRegion$gene_id)]
strd <- as.character(strand(gr2))
len <- GRegionLen[GRegionLen != 0]
GRegion$exon_rank <- lapply(seq_along(strd), function(i) {
  rank <- seq(1, len[i])
  if (strd[i] == "-") 
    rank <- rev(rank)
  return(rank)
}) %>% unlist

txidinfo <- transcripts(gencode.v29.txdb, columns = c("tx_id", "tx_name", "gene_id"))
idx <- which(sapply(txidinfo$gene_id, length) == 0)
txidinfo[idx, ]$gene_id <- txidinfo[idx, ]$tx_name
txid2geneid <- paste(mcols(txidinfo)[["tx_name"]], mcols(txidinfo)[["gene_id"]], sep = "/")
txid2geneid <- sub("/NA", "", txid2geneid)
names(txid2geneid) <- mcols(txidinfo)[["tx_id"]]
nn<-txid2geneid[names(genomicRegion)]
#nn <- ChIPseeker:::TXID2EG(names(genomicRegion))
names(GRegionLen) <- nn
GRegion$gene_id <- rep(nn, times = GRegionLen)
anno <- paste("Exon", " (", GRegion$gene_id, ", exon ", GRegion$exon_rank, " of ", GRegionLen[GRegion$gene_id], ")", sep = "")

g4.exon<-data[ids.exon,]
which(anno %in% g4.exon$annotation[2])

exon.width<-width(GRegion)
names(exon.width)<-anno

g4.count.each.exon<-g4.exon$annotation %>% table %>% as.data.frame()

notin.genes<-anno[!(anno %in% g4.count.each.exon[,1])]
tmp<-data.frame("."=notin.genes,"Freq"=rep(0,length(notin.genes)))
g4.count.each.exon<-rbind(g4.count.each.exon,tmp)

row.names(g4.count.each.exon)<-g4.count.each.exon[,1]
g4.count.each.exon[,1]<-NULL
g4.count.each.exon$width<-exon.width[row.names(g4.count.each.exon)]

g4.count.each.exon$gene_id<-sub("^Exon \\(ENST[0-9]+\\.[0-9]+/(ENSG[0-9]+\\.[0-9]+),.+","\\1",row.names(g4.count.each.exon),perl = T)
g4.count.each.exon$tx_name<-sub("^Exon \\((ENST[0-9]+\\.[0-9]+)/ENSG[0-9]+\\.[0-9]+,.+","\\1",row.names(g4.count.each.exon),perl = T)
g4.count.each.exon<-g4.count.each.exon[,c(3,4,1,2)]
colnames(g4.count.each.exon)[3]<-"Count"

g4.count.exon.dummy<-group_by(g4.count.each.exon,gene_id) %>% summarise(Count=sum(Count),Width=sum(width))

for(name in names(toptags.genes.ensembl)){
  g4.count.exon.dummy[,name]<-g4.count.exon.dummy$gene_id %in% toptags.genes.ensembl[[name]] %>% as.integer() %>% as.factor()
}
g4.count.exon.dummy

g4.count.exon.df<-g4.count.exon.dummy[,1:3]
g4.count.exon.df$Category <-"all"
for(name in names(toptags.genes.ensembl)){
  tmp<-g4.count.exon.dummy[which(g4.count.exon.dummy[,name]==1),1:3]
  tmp$Category <-name
  g4.count.exon.df<-rbind(g4.count.exon.df,tmp)
}
boxplot(log10(Count+1)~Category,data=g4.count.exon.df)
boxplot(Count/Width~Category,data=g4.count.exon.df,ylim=c(0,0.001))

group_by(g4.count.exon.df,Category) %>% summarise(median=median(Count),mean=mean(Count))
group_by(g4.count.exon.df,Category) %>% summarise(median=median(Count/Width),mean=mean(Count/Width))



library(ggridges)
p.riges <- ggplot(g4.count.exon.df,aes(x=log10(Count+1),y=Category,fill=Category)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2, vline_color = "red",alpha = .7,from=0)+
  scale_x_continuous(name="log10(PQS Counts)") +
  scale_y_discrete(limits=c("all",names(toptags.genes.ensembl))) +
  scale_fill_cyclical(values =c("gray",rep(c("#8080ff","#ff8080"),4)))+
  theme_bw()

p.riges
ggsave(file = "./out/geom_density_PQScounts_exon.pdf", plot = p.riges, width = 11, height = 8.5)

p.riges <- ggplot(g4.count.exon.df,aes(x=Count/Width,y=Category,fill=Category)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2, vline_color = "red",alpha = .7,from=0,to=0.001)+
  scale_x_continuous(name="PQS Counts/Width") +
  scale_y_discrete(limits=c("all",names(toptags.genes.ensembl))) +
  scale_fill_cyclical(values =c("gray",rep(c("#8080ff","#ff8080"),4)))+
  theme_bw()

p.riges
ggsave(file = "./out/geom_density_PQScounts_div_genelen_exon.pdf", plot = p.riges, width = 11, height = 8.5)


glm.exon.g4<-list()
for(name in names(toptags.genes.ensembl)){
  glm.exon.g4[[name]] <- glm(formula = as.formula(paste0("Count~",name)),data=g4.count.exon.dummy,family = poisson())
}
lapply(glm.exon.g4, summary)

glm.offset.exon.g4<-list()
for(name in names(toptags.genes.ensembl)){
  glm.offset.exon.g4[[name]] <- glm(formula = as.formula(paste0("Count~",name)),offset = log(Width),data=g4.count.exon.dummy,family = poisson())
}
lapply(glm.offset.exon.g4, summary)


g4.count.exon.reg<-data.frame(Category=NULL,gene=NULL,Count=NULL,Width=NULL,reg=NULL)
for(name in names(toptags.genes.ensembl)){
  tmp<-data.frame(Category=name,gene=g4.count.exon.dummy$gene_id,Count=g4.count.exon.dummy$Count,Width=g4.count.exon.dummy$Width,reg=g4.count.exon.dummy[,name])
  colnames(tmp)[5]<-"reg"
  g4.count.exon.reg<-rbind(g4.count.exon.reg,tmp)
}

#write.csv(g4.count.exon.reg,file="./out/g4.count.exon.reg.csv",row.names = F,quote = F)

#1st exon
ids.1stexon<-grep("exon 1 of",row.names(g4.count.each.exon))
ids.otherexon<-grep("exon 1 of",row.names(g4.count.each.exon),invert = T)
g4.count.each.exon[ids.1stexon,]

g4.count.1stexon.dummy<-group_by(g4.count.each.exon[ids.1stexon,],gene_id) %>% summarise(Count=sum(Count),Width=sum(width))

for(name in names(toptags.genes.ensembl)){
  g4.count.1stexon.dummy[,name]<-g4.count.1stexon.dummy$gene_id %in% toptags.genes.ensembl[[name]] %>% as.integer() %>% as.factor()
}
g4.count.1stexon.dummy

g4.count.1stexon.df<-g4.count.1stexon.dummy[,1:3]
g4.count.1stexon.df$Category <-"all"
for(name in names(toptags.genes.ensembl)){
  tmp<-g4.count.1stexon.dummy[which(g4.count.1stexon.dummy[,name]==1),1:3]
  tmp$Category <-name
  g4.count.1stexon.df<-rbind(g4.count.1stexon.df,tmp)
}

group_by(g4.count.1stexon.df,Category) %>% summarise(median=median(Count),mean=mean(Count))
group_by(g4.count.1stexon.df,Category) %>% summarise(median=median(Count/Width),mean=mean(Count/Width))

boxplot(log10(Count+1)~Category,data=g4.count.1stexon.df)
boxplot(Count/Width~Category,data=g4.count.exon.df,ylim=c(0,0.001))

library(ggridges)
p.riges <- ggplot(g4.count.1stexon.df,aes(x=log10(Count+1),y=Category,fill=Category)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2, vline_color = "red",alpha = .7,from=0)+
  scale_x_continuous(name="log10(PQS Counts)") +
  scale_y_discrete(limits=c("all",names(toptags.genes.ensembl))) +
  scale_fill_cyclical(values =c("gray",rep(c("#8080ff","#ff8080"),4)))+
  theme_bw()

p.riges
ggsave(file = "./out/geom_density_PQScounts_1stexon.pdf", plot = p.riges, width = 11, height = 8.5)

p.riges <- ggplot(g4.count.1stexon.df,aes(x=Count/Width,y=Category,fill=Category)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2, vline_color = "red",alpha = .7,from=0,to=0.001)+
  scale_x_continuous(name="PQS Counts/Width") +
  scale_y_discrete(limits=c("all",names(toptags.genes.ensembl))) +
  scale_fill_cyclical(values =c("gray",rep(c("#8080ff","#ff8080"),4)))+
  theme_bw()

p.riges
ggsave(file = "./out/geom_density_PQScounts_div_genelen_1stexon.pdf", plot = p.riges, width = 11, height = 8.5)

glm.1stexon.g4<-list()
for(name in names(toptags.genes.ensembl)){
  glm.1stexon.g4[[name]] <- glm(formula = as.formula(paste0("Count~",name)),data=g4.count.1stexon.dummy,family = poisson())
}
lapply(glm.1stexon.g4, summary)

glm.offset.1stexon.g4<-list()
for(name in names(toptags.genes.ensembl)){
  glm.offset.1stexon.g4[[name]] <- glm(formula = as.formula(paste0("Count~",name)),offset = log(Width),data=g4.count.1stexon.dummy,family = poisson())
}
lapply(glm.offset.1stexon.g4, summary)


g4.count.1stexon.reg<-data.frame(Category=NULL,gene=NULL,Count=NULL,Width=NULL,reg=NULL)
for(name in names(toptags.genes.ensembl)){
  tmp<-data.frame(Category=name,gene=g4.count.1stexon.dummy$gene_id,Count=g4.count.1stexon.dummy$Count,Width=g4.count.1stexon.dummy$Width,reg=g4.count.1stexon.dummy[,name])
  colnames(tmp)[5]<-"reg"
  g4.count.1stexon.reg<-rbind(g4.count.1stexon.reg,tmp)
}

#write.csv(g4.count.1stexon.reg,file="./out/g4.count.1stexon.reg.csv",row.names = F,quote = F)

#other exon
ids.otherexon<-grep("exon 1 of",row.names(g4.count.each.exon),invert = T)

g4.count.otherexon.dummy<-group_by(g4.count.each.exon[ids.otherexon,],gene_id) %>% summarise(Count=sum(Count),Width=sum(width))

for(name in names(toptags.genes.ensembl)){
  g4.count.otherexon.dummy[,name]<-g4.count.otherexon.dummy$gene_id %in% toptags.genes.ensembl[[name]] %>% as.integer() %>% as.factor()
}
g4.count.otherexon.dummy

g4.count.otherexon.df<-g4.count.otherexon.dummy[,1:3]
g4.count.otherexon.df$Category <-"all"
for(name in names(toptags.genes.ensembl)){
  tmp<-g4.count.otherexon.dummy[which(g4.count.otherexon.dummy[,name]==1),1:3]
  tmp$Category <-name
  g4.count.otherexon.df<-rbind(g4.count.otherexon.df,tmp)
}

group_by(g4.count.otherexon.df,Category) %>% summarise(median=median(Count),mean=mean(Count))
group_by(g4.count.otherexon.df,Category) %>% summarise(median=median(Count/Width),mean=mean(Count/Width))

boxplot(log10(Count+1)~Category,data=g4.count.otherexon.df)
boxplot(Count/Width~Category,data=g4.count.otherexon.df,ylim=c(0,0.001))

library(ggridges)
p.riges <- ggplot(g4.count.otherexon.df,aes(x=log10(Count+1),y=Category,fill=Category)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2, vline_color = "red",alpha = .7,from=0)+
  scale_x_continuous(name="PQS Counts") +
  scale_y_discrete(limits=c("all",names(toptags.genes.ensembl))) +
  scale_fill_cyclical(values =c("gray",rep(c("#8080ff","#ff8080"),4)))+
  theme_bw()

p.riges
ggsave(file = "./out/geom_density_PQScounts_otherexon.pdf", plot = p.riges, width = 11, height = 8.5)

p.riges <- ggplot(g4.count.otherexon.df,aes(x=Count/Width,y=Category,fill=Category)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2, vline_color = "red",alpha = .7,from=0,to=0.001)+
  scale_x_continuous(name="PQS Counts/Width") +
  scale_y_discrete(limits=c("all",names(toptags.genes.ensembl))) +
  scale_fill_cyclical(values =c("gray",rep(c("#8080ff","#ff8080"),4)))+
  theme_bw()

p.riges
ggsave(file = "./out/geom_density_PQScounts_div_genelen_otherexon.pdf", plot = p.riges, width = 11, height = 8.5)


glm.otherexon.g4<-list()
for(name in names(toptags.genes.ensembl)){
  glm.otherexon.g4[[name]] <- glm(formula = as.formula(paste0("Count~",name)),data=g4.count.otherexon.dummy,family = poisson())
}
lapply(glm.otherexon.g4, summary)

glm.offset.otherexon.g4<-list()
for(name in names(toptags.genes.ensembl)){
  glm.offset.otherexon.g4[[name]] <- glm(formula = as.formula(paste0("Count~",name)),offset = log(Width),data=g4.count.otherexon.dummy,family = poisson())
}
lapply(glm.offset.otherexon.g4, summary)

g4.count.otherexon.reg<-data.frame(Category=NULL,gene=NULL,Count=NULL,Width=NULL,reg=NULL)
for(name in names(toptags.genes.ensembl)){
  tmp<-data.frame(Category=name,gene=g4.count.otherexon.dummy$gene_id,Count=g4.count.otherexon.dummy$Count,Width=g4.count.otherexon.dummy$Width,reg=g4.count.otherexon.dummy[,name])
  colnames(tmp)[5]<-"reg"
  g4.count.otherexon.reg<-rbind(g4.count.otherexon.reg,tmp)
}

#write.csv(g4.count.otherexon.reg,file="./out/g4.count.otherexon.reg.csv",row.names = F,quote = F)

#5'UTR-----------------------------------------------------------------------
data<-peakAnnoList[["gencode.v29"]]@anno %>% as.data.frame()
ids.5utr <- grep("5' UTR",data$annotation)
#data[ids.5utr,] %>% head(.,100)

#ChPseeker:::getGenomicAnnotation.internal
genomicRegion<-fiveUTRList <- fiveUTRsByTranscript(gencode.v29.txdb)
GRegion <- unlist(genomicRegion)
GRegionLen <- elementNROWS(genomicRegion)
# names(GRegionLen) <- names(genomicRegion)
# GRegion$gene_id <- rep(names(genomicRegion), times = GRegionLen)

txidinfo <- transcripts(gencode.v29.txdb, columns = c("tx_id", "tx_name", "gene_id"))
idx <- which(sapply(txidinfo$gene_id, length) == 0)
txidinfo[idx, ]$gene_id <- txidinfo[idx, ]$tx_name
txid2geneid <- unlist(mcols(txidinfo)[["gene_id"]])
txid2geneid <- sub("NA", "", txid2geneid)
names(txid2geneid) <- mcols(txidinfo)[["tx_id"]]
nn<-txid2geneid[names(genomicRegion)]
#nn <- ChIPseeker:::TXID2EGID(names(genomicRegion))
names(GRegionLen) <- nn
GRegion$gene_id <- rep(nn, times = GRegionLen)


fiveUTR.peak<-peakAnnoList[["gencode.v29"]]@anno[ids.5utr]
ol<-findOverlaps(fiveUTR.peak,GRegion)

##ChIPseekerのアルゴリズムではpeakに対して最初にHitした遺伝子領域にしぼっている。
qh <- queryHits(ol)
hit.idx <- which(!duplicated(qh))
ol <- ol[hit.idx]

queryIndex <- queryHits(ol)
subjectIndex <- subjectHits(ol)

data.frame(fiveUTR.peak[queryIndex]$geneId,GRegion[subjectIndex]$gene_id)
#example
#ENSG00000182533.6, ENSG00000180914.10
# #hitしたUTRと出力されたgene_idと同じものを抜き出そうとしたが、異
# peak.geneid<-fiveUTR.peak$geneId
# GR.geneid<-GRegion$gene_id
# ol.macth<-c()
# for(i in 1:20){
#   print(fiveUTR.peak[queryIndex[i]])
#   print(GRegion[subjectIndex[i]])
#   #ol.macth<-c(ol.macth,peak.geneid[queryIndex[i]] == GR.geneid[subjectIndex[i]])
# }
# 
# #この遺伝子には5'UTRがない
# GRegion[GRegion$gene_id == "ENSG00000255054.3"]


fiveUTR.peak$width<-width(GRegion[subjectIndex])
mcols(fiveUTR.peak)$gene_id.corr<-as.character(GRegion[subjectIndex]$gene_id)

g4.count.5utr<-fiveUTR.peak$gene_id.corr %>% table

gencode.genes<-genes(gencode.v29.txdb)$gene_id
notin.genes<-gencode.genes[!(gencode.genes %in% names(g4.count.5utr))]
x<-rep(0,length(notin.genes))
names(x)<-notin.genes
g4.count.5utr<-c(g4.count.5utr,x)

g4.count.5utr.genelen<-fiveUTR.peak$width
names(g4.count.5utr.genelen)<-fiveUTR.peak$gene_id.corr
g4.count.5utr.genelen<-g4.count.5utr.genelen[!duplicated(names(g4.count.5utr.genelen))]

notinGR<-GRegion[!(GRegion$gene_id %in% names(g4.count.5utr.genelen))]
notinGR<-notinGR[!duplicated(notinGR$gene_id)]
xxx<-width(notinGR)
names(xxx)<-notinGR$gene_id
g4.count.5utr.genelen<-c(g4.count.5utr.genelen,xxx)

#5UTRないものは削除
g4.count.5utr<-g4.count.5utr[names(g4.count.5utr.genelen)]

# notin.genes<-gencode.genes[!(gencode.genes %in% names(g4.count.5utr.genelen))]
# x<-rep(0,length(notin.genes))
# names(x)<-notin.genes
# g4.count.5utr.genelen<-c(g4.count.5utr.genelen,x)

g4.count.5utr.dummy<-data.frame(Count=g4.count.5utr,Width=g4.count.5utr.genelen)
g4.count.5utr.dummy$gene_id<-row.names(g4.count.5utr.dummy)
g4.count.5utr.dummy<-as.tibble(g4.count.5utr.dummy[,c(3,1,2)])

for(name in names(toptags.genes.ensembl)){
  g4.count.5utr.dummy[,name]<-g4.count.5utr.dummy$gene_id %in% toptags.genes.ensembl[[name]] %>% as.integer() %>% as.factor()
}
summary(g4.count.5utr.dummy)

g4.count.5utr.df<-g4.count.5utr.dummy[,1:3]
g4.count.5utr.df$Category <-"all"
for(name in names(toptags.genes.ensembl)){
  tmp<-g4.count.5utr.dummy[which(g4.count.5utr.dummy[,name]==1),1:3]
  tmp$Category <-name
  g4.count.5utr.df<-rbind(g4.count.5utr.df,tmp)
}

boxplot(log10(Count+1)~Category,data=g4.count.5utr.df)
boxplot(Count/Width~Category,data=g4.count.5utr.df,ylim=c(0,0.01))


group_by(g4.count.5utr.df,Category) %>% summarise(median=median(Count),mean=mean(Count))
group_by(g4.count.5utr.df,Category) %>% summarise(median=median(Count/Width),mean=mean(Count/Width))

library(ggridges)
p.riges <- ggplot(g4.count.5utr.df,aes(x=log10(Count+1),y=Category,fill=Category)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2, vline_color = "red",alpha = .7,from=0)+
  scale_x_continuous(name="PQS Counts") +
  scale_y_discrete(limits=c("all",names(toptags.genes.ensembl))) +
  scale_fill_cyclical(values =c("gray",rep(c("#8080ff","#ff8080"),4)))+
  theme_bw()

p.riges
ggsave(file = "./out/geom_density_PQScounts_5utr.pdf", plot = p.riges, width = 11, height = 8.5)

p.riges <- ggplot(g4.count.5utr.df,aes(x=Count/(Width+10),y=Category,fill=Category)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2, vline_color = "red",alpha = .7,from=0,to=0.1)+
  scale_x_continuous(name="PQS Counts/Width") +
  scale_y_discrete(limits=c("all",names(toptags.genes.ensembl))) +
  scale_fill_cyclical(values =c("gray",rep(c("#8080ff","#ff8080"),4)))+
  theme_bw()

p.riges
ggsave(file = "./out/geom_density_PQScounts_div_genelen_5utr.pdf", plot = p.riges, width = 11, height = 8.5)



glm.5utr.g4<-list()
for(name in names(toptags.genes.ensembl)){
  glm.5utr.g4[[name]] <- glm(formula = as.formula(paste0("Count~",name)),data=g4.count.5utr.dummy,family = poisson())
}
lapply(glm.5utr.g4, summary)

glm.offset.5utr.g4<-list()
for(name in names(toptags.genes.ensembl)){
  glm.offset.5utr.g4[[name]] <- glm(formula = as.formula(paste0("Count~",name)),offset = log(Width),data=g4.count.5utr.dummy,family = poisson())
}
lapply(glm.offset.5utr.g4, summary)

g4.count.5utr.reg<-data.frame(Category=NULL,gene=NULL,Count=NULL,Width=NULL,reg=NULL)
for(name in names(toptags.genes.ensembl)){
  tmp<-data.frame(Category=name,gene=g4.count.5utr.dummy$gene_id,Count=g4.count.5utr.dummy$Count,Width=g4.count.5utr.dummy$Width,reg=g4.count.5utr.dummy[,name])
  colnames(tmp)[5]<-"reg"
  g4.count.5utr.reg<-rbind(g4.count.5utr.reg,tmp)
}

#write.csv(g4.count.5utr.reg,file="./out/g4.count.5utr.reg.csv",row.names = F,quote = F)

#3'UTR--------------------------------------------------------
data<-peakAnnoList[["gencode.v29"]]@anno %>% as.data.frame()
ids.3utr <- grep("3' UTR",data$annotation)
#data[ids.3utr,] %>% head(.,100)

#ChPseeker:::getGenomicAnnotation.internal
genomicRegion<-threeUTRList <- threeUTRsByTranscript(gencode.v29.txdb)
GRegion <- unlist(genomicRegion)
GRegionLen <- elementNROWS(genomicRegion)

txidinfo <- transcripts(gencode.v29.txdb, columns = c("tx_id", "tx_name", "gene_id"))
idx <- which(sapply(txidinfo$gene_id, length) == 0)
txidinfo[idx, ]$gene_id <- txidinfo[idx, ]$tx_name
txid2geneid <- unlist(mcols(txidinfo)[["gene_id"]])
txid2geneid <- sub("NA", "", txid2geneid)
names(txid2geneid) <- mcols(txidinfo)[["tx_id"]]
nn<-txid2geneid[names(genomicRegion)]
#nn <- ChIPseeker:::TXID2EGID(names(genomicRegion))
names(GRegionLen) <- nn
GRegion$gene_id <- rep(nn, times = GRegionLen)


threeUTR.peak<-peakAnnoList[["gencode.v29"]]@anno[ids.3utr]
ol<-findOverlaps(threeUTR.peak,GRegion)

##ChIPseekerのアルゴリズムではpeakに対して最初にHitした遺伝子領域にしぼっている。
qh <- queryHits(ol)
hit.idx <- which(!duplicated(qh))
ol <- ol[hit.idx]

queryIndex <- queryHits(ol)
subjectIndex <- subjectHits(ol)

data.frame(threeUTR.peak[queryIndex]$geneId,GRegion[subjectIndex]$gene_id)


threeUTR.peak$width<-width(GRegion[subjectIndex])
mcols(threeUTR.peak)$gene_id.corr<-as.character(GRegion[subjectIndex]$gene_id)

g4.count.3utr<-threeUTR.peak$gene_id.corr %>% table

gencode.genes<-genes(gencode.v29.txdb)$gene_id
notin.genes<-gencode.genes[!(gencode.genes %in% names(g4.count.3utr))]
x<-rep(0,length(notin.genes))
names(x)<-notin.genes
g4.count.3utr<-c(g4.count.3utr,x)

g4.count.3utr.genelen<-threeUTR.peak$width
names(g4.count.3utr.genelen)<-threeUTR.peak$gene_id.corr
g4.count.3utr.genelen<-g4.count.3utr.genelen[!duplicated(names(g4.count.3utr.genelen))]

notinGR<-GRegion[!(GRegion$gene_id %in% names(g4.count.3utr.genelen))]
notinGR<-notinGR[!duplicated(notinGR$gene_id)]
xxx<-width(notinGR)
names(xxx)<-notinGR$gene_id
g4.count.3utr.genelen<-c(g4.count.3utr.genelen,xxx)

#3UTRないものは削除
g4.count.3utr<-g4.count.3utr[names(g4.count.3utr.genelen)]


g4.count.3utr.dummy<-data.frame(Count=g4.count.3utr,Width=g4.count.3utr.genelen)
g4.count.3utr.dummy$gene_id<-row.names(g4.count.3utr.dummy)
g4.count.3utr.dummy<-as.tibble(g4.count.3utr.dummy[,c(3,1,2)])

for(name in names(toptags.genes.ensembl)){
  g4.count.3utr.dummy[,name]<-g4.count.3utr.dummy$gene_id %in% toptags.genes.ensembl[[name]] %>% as.integer() %>% as.factor()
}
summary(g4.count.3utr.dummy)

g4.count.3utr.df<-g4.count.3utr.dummy[,1:3]
g4.count.3utr.df$Category <-"all"
for(name in names(toptags.genes.ensembl)){
  tmp<-g4.count.3utr.dummy[which(g4.count.3utr.dummy[,name]==1),1:3]
  tmp$Category <-name
  g4.count.3utr.df<-rbind(g4.count.3utr.df,tmp)
}

boxplot(log10(Count+1)~Category,data=g4.count.3utr.df)
boxplot(Count/Width~Category,data=g4.count.3utr.df,ylim=c(0,0.01))

group_by(g4.count.3utr.df,Category) %>% summarise(median=median(Count),mean=mean(Count))
group_by(g4.count.3utr.df,Category) %>% summarise(median=median(Count/Width),mean=mean(Count/Width))

library(ggridges)
p.riges <- ggplot(g4.count.3utr.df,aes(x=log10(Count+1),y=Category,fill=Category)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2, vline_color = "red",alpha = .7,from=0)+
  scale_x_continuous(name="PQS Counts") +
  scale_y_discrete(limits=c("all",names(toptags.genes.ensembl))) +
  scale_fill_cyclical(values =c("gray",rep(c("#8080ff","#ff8080"),4)))+
  theme_bw()

p.riges
ggsave(file = "./out/geom_density_PQScounts_3utr.pdf", plot = p.riges, width = 11, height = 8.5)

p.riges <- ggplot(g4.count.3utr.df,aes(x=Count/(Width+0.01),y=Category,fill=Category)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2, vline_color = "red",alpha = .7,from=0,to=0.02)+
  scale_x_continuous(name="PQS Counts/Width") +
  scale_y_discrete(limits=c("all",names(toptags.genes.ensembl))) +
  scale_fill_cyclical(values =c("gray",rep(c("#8080ff","#ff8080"),4)))+
  theme_bw()

p.riges
ggsave(file = "./out/geom_density_PQScounts_div_genelen_3utr.pdf", plot = p.riges, width = 11, height = 8.5)



glm.3utr.g4<-list()
for(name in names(toptags.genes.ensembl)){
  glm.3utr.g4[[name]] <- glm(formula = as.formula(paste0("Count~",name)),data=g4.count.3utr.dummy,family = poisson())
}
lapply(glm.3utr.g4, summary)

glm.offset.3utr.g4<-list()
for(name in names(toptags.genes.ensembl)){
  glm.offset.3utr.g4[[name]] <- glm(formula = as.formula(paste0("Count~",name)),offset = log(Width),data=g4.count.3utr.dummy,family = poisson())
}
lapply(glm.offset.3utr.g4, summary)

g4.count.3utr.reg<-data.frame(Category=NULL,gene=NULL,Count=NULL,Width=NULL,reg=NULL)
for(name in names(toptags.genes.ensembl)){
  tmp<-data.frame(Category=name,gene=g4.count.3utr.dummy$gene_id,Count=g4.count.3utr.dummy$Count,Width=g4.count.3utr.dummy$Width,reg=g4.count.3utr.dummy[,name])
  colnames(tmp)[5]<-"reg"
  g4.count.3utr.reg<-rbind(g4.count.3utr.reg,tmp)
}

#write.csv(g4.count.3utr.reg,file="./out/g4.count.3utr.reg.csv",row.names = F,quote = F)

