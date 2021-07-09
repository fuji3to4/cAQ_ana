library(edgeR)
library(tidyverse)


cAQ_treat<-read_csv("cAQ_treat_fileids.csv")
cAQ_treat$name<-factor(cAQ_treat$name,
                       levels=c("Control.0.6","Control.0.24","AQ_ac.1.6","AQ_ac.1.24","cAQ_mBen.1.6","cAQ_mBen.1.24","cAQ_mBen.5.6","cAQ_mBen.5.24"))
#cAQ_treat$name<-relevel(cAQ_treat$name,ref="Control.0.6")

counts_list<-list()
for(id in cAQ_treat$libraryID){
  counts_list[[id]]<-read_tsv(paste0("./RNA-Seq/HNWHLDSXY_",id,"/Aligned.count.htseq.txt"),col_names = c("Ensembl","Gene_name","Counts"))
}

#save(counts_list,file = "cAQ_RNA-Seq_cout.rda")
#load("cAQ_RNA-Seq_cout.rda")

tmp<-as.data.frame(counts_list)
gene<-tmp[,1:2]
colnames(gene)<-c("Ensembl","Gene_name")

#By EnsDb.Hsapiens.v86
library(EnsDb.Hsapiens.v86)
columns(EnsDb.Hsapiens.v86)

gene$genebiotype<-mapIds(EnsDb.Hsapiens.v86,
                  keys=substr(gene$Ensembl,1,15),
                  column="GENEBIOTYPE",
                  keytype="GENEID",
                  multiVals="first")



cAQ_dge<-DGEList(tmp[,seq(3,ncol(tmp),3)],genes=gene,group = cAQ_treat$name)
colnames(cAQ_dge)<-sub(".Counts","",colnames(cAQ_dge))
rm(tmp)


design <-model.matrix(~0+name,data=cAQ_treat)
colnames(design)<- colnames(design) %>% sub("name","",.)
keep<-filterByExpr(cAQ_dge,design)
cAQ_dge<-cAQ_dge[keep, , keep.lib.sizes=FALSE]

cAQ_dge<-calcNormFactors(cAQ_dge)
cAQ_dge$samples

cAQ_g<-cAQ_treat$name
levels(cAQ_g)<-1:8
plotMDS(cAQ_dge,col=as.integer(cAQ_g))
plotMD(cAQ_dge)

cAQ_g<-as.factor(cAQ_treat$additives)
levels(cAQ_g)<-1:3
plotMDS(cAQ_dge,col=as.integer(cAQ_g),labels =as.character(cAQ_treat$name),cex=.6)

cAQ_dge <- estimateDisp(cAQ_dge, design, robust=TRUE)
cAQ_dge$common.dispersion
plotBCV(cAQ_dge)

#Grouping-----------------------------------------

fit <- glmQLFit(cAQ_dge, design, robust=TRUE)
plotQLDisp(fit)
qlf <- glmQLFTest(fit,coef=1:6)
topTags(qlf)
#is.de <- decideTestsDGE(qlf)
#summary(is.de)
#plotMD(qlf, status=is.de, values=c(1,-1), col=c("red","blue"),legend="topright")

my.contrasts<-makeContrasts(
  AQ_ac.1uM.6h = AQ_ac.1.6 - Control.0.6,
  AQ_ac.1uM.24h = AQ_ac.1.24 - Control.0.24,
  cAQ_mBen.1uM.6h = cAQ_mBen.1.6 - Control.0.6,
  cAQ_mBen.1uM.24h = cAQ_mBen.1.24 - Control.0.24,
  cAQ_mBen.5uM.6h = cAQ_mBen.5.6 - Control.0.6,
  cAQ_mBen.5uM.24h = cAQ_mBen.5.24 - Control.0.24, 
  cAQ_mBenvsAQ_ac.1uM.6h = cAQ_mBen.1.6 - AQ_ac.1.6,
  cAQ_mBenvsAQ_ac.1uM.24h = cAQ_mBen.1.24 - AQ_ac.1.24,
  cAQ_mBen.1uMvs5uM.6h = cAQ_mBen.5.6 - cAQ_mBen.1.6,
  cAQ_mBen.1uMvs5uM.24h = cAQ_mBen.5.24 - cAQ_mBen.1.24,
  Control.6hvs24h = Control.0.6 - Control.0.24,
  AQ_ac.1uM.6hvs24h = (AQ_ac.1.6 - Control.0.6) - (AQ_ac.1.24 - Control.0.24),
  cAQ_mBen.1uM.6hvs24h = (cAQ_mBen.1.6 - Control.0.6) - (cAQ_mBen.1.24 - Control.0.24),
  cAQ_mBen.5uM.6hvs24h = (cAQ_mBen.5.6 - Control.0.6) - (cAQ_mBen.5.24 - Control.0.24),
  cAQ_mBen.1uMvs5uM.6hvs24h = (cAQ_mBen.1.24 - cAQ_mBen.1.6) - (cAQ_mBen.5.24 - cAQ_mBen.5.6),
  levels=design)
  
#options(matprod = "internal")
# res <- glmQLFTest(fit,contrast = my.contrasts[,"cAQ_mBen.5uM.6h"])
# tt<-topTags(res,p.value = 0.05,n=NULL)
# is.de <- decideTestsDGE(res)
# summary(is.de)
# plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"),legend="topright")

pdf("./out/plotMD_cAQ_FDR0.05.pdf")
res<-list()
deg.summary<-list()
for(nn in colnames(my.contrasts)){
  res[[nn]] <- glmQLFTest(fit,contrast = my.contrasts[,nn])
  is.de <- decideTestsDGE(res[[nn]],p.value = 0.05)
  su<-summary(is.de)
  deg.summary[[nn]]<-su
  # if(su[1]+su[3] != 0){
  #   topTags(res[[nn]],p.value = 0.05,n=su[1]+su[3]) %>% as.data.frame() %>% write_csv(.,file=paste0("./out/toptags_",nn,".csv"))
  # }else{
  #   topTags(res[[nn]],n=10) %>% as.data.frame() %>% write_csv(.,file=paste0("./out/toptags_",nn,".csv"))
  # }
  plotMD(res[[nn]], status=is.de, values=c(1,-1), col=c("red","blue"),legend="topright")
}
dev.off()

deg.summary.df<-do.call(cbind,deg.summary) %>% t %>% as.data.frame()
#Reduce(function(a,b)cbind(a,b),deg.summary)
write.csv(deg.summary.df,file = "./out/deg_summary_fdr0.05.csv",)

save(cAQ_treat,counts_list,res,fit,cAQ_dge,file="cAQ_RNA-Seq_de.RData")

#plot CPM----------------------------------------------
plot.id<-function(id,dge,...){
  if(length(grep("^ENS",id))==1){
    n<-which(dge$genes$Ensembl == id)
  }else if(length(grep("^[0-9]",id))==1){
    n<-which(dge$genes$entrez == id)
  }else{
    n<-which(dge$genes$Gene_name == toupper(id))
  }
  
  dge$cpm<-cpm(dge,log = T)
  
  plot(jitter(as.numeric(dge$samples$group)),dge$cpm[n,],xaxt="n", 
       col=as.numeric(factor(sub(".[0-9]*$","",dge$samples$group),levels=c("Control.0","AQ_ac.1","cAQ_mBen.1","cAQ_mBen.5"))),
       main=paste(dge$genes[n,1],dge$genes[n,2],dge$genes[n,3],sep=" / "),xlab="Treat",ylab="logCPM")
  axis(1,at=1:length(levels(dge$samples$group)),labels = levels(dge$samples$group))
}



plot.id("MYC",cAQ_dge)

for(nn in colnames(my.contrasts)){
  tt<-topTags(res[[nn]],n=20)
  pdf(paste0("./out/plot_logCPM_top20_",nn,".pdf"),width=12, height=10)
  for(x in tt$table$Ensembl){
    plot.id(x,cAQ_dge)
  }
  dev.off()
}

#Venn Diagram--------------------------
toptags.genes.id<-list(cAQ_mBen.1uM.6h=topTags(res[["cAQ_mBen.1uM.6h"]],p.value = 0.05,n=NULL)$table$entrez,
                       cAQ_mBen.1uM.24h=topTags(res[["cAQ_mBen.1uM.24h"]],p.value = 0.05,n=NULL)$table$entrez,
                       cAQ_mBen.5uM.6h=topTags(res[["cAQ_mBen.5uM.6h"]],p.value = 0.05,n=NULL)$table$entrez,
                       cAQ_mBen.5uM.24h=topTags(res[["cAQ_mBen.5uM.24h"]],p.value = 0.05,n=NULL)$table$entrez)

toptags.genes.id<-lapply(toptags.genes.id, as.numeric)
toptags.genes.id<-lapply(toptags.genes.id, function(x)x[!is.na(x)])
library(UpSetR)
library(ggvenn)
pdf("./out/DEG/deg_venn_fdr0.05.pdf",width = 10,height = 10)
upset(fromList(toptags.genes.id),order.by = "freq")

ggvenn(
  toptags.genes.id, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
#  fill_color = c("#D55E0050","#0072B250", "#FF000050","#0000FF50"),
  stroke_size = 0.5, set_name_size = 4
)
dev.off()

#combine PQS and toptags----------------------
source("r_cAQ_PQS_ChIPseeker.R")
source("r_cAQ_PQS_promoter.R")
source("r_cAQ_PQS_genebody.R")
g4.count.merge<-merge(g4.count.u1kd1k.dummy[,1:2],g4.count.genebody.dummy[,1:3],by="gene")
colnames(g4.count.merge)<-c("Ensembl","g4.count.u1kd1k","g4.count.genebody","gene_length")

for(name in c("cAQ_mBen.1uM.6h","cAQ_mBen.1uM.24h","cAQ_mBen.5uM.6h","cAQ_mBen.5uM.24h")){
  tmp<-topTags(res[[name]],n=NULL)$table[,c("Ensembl","logFC","logCPM","FDR")]
  colnames(tmp)[2:4]<-paste(colnames(tmp)[2:4],name,sep=".")
  g4.count.merge<-merge(g4.count.merge,tmp,by="Ensembl")
}
tmp<-topTags(res[[name]],n=NULL)$table[,c("Ensembl","Gene_name")]

g4.count.merge<-merge(g4.count.merge,tmp,by="Ensembl")
ncol(g4.count.merge)
#g4.count.merge<-g4.count.merge[,c(1,13,4,2,3,5:12)]


library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
filters = listFilters(ensembl)
filters
attributes = listAttributes(ensembl)
attributes
g4.count.merge$ensembl_gene_id<-substr(g4.count.merge$Ensembl,1,15)
x<-getBM(attributes=c('ensembl_gene_id','description'), 
         filters = 'ensembl_gene_id', 
         values = g4.count.merge$ensembl_gene_id, 
         mart = ensembl)
g4.count.merge.withdesc<-merge(g4.count.merge,x,by="ensembl_gene_id",all.x=TRUE)

write_csv(g4.count.merge.withdesc,file="./out/g4.count.merge.csv")

