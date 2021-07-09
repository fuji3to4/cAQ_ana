library(edgeR)
library(tidyverse)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(org.Hs.eg.db)
library(gplots)
library(ggnewscale)
load("cAQ_RNA-Seq_de.RData")




toptags.genes.id<-list(cAQ_mBen.1uM.6h=topTags(res[["cAQ_mBen.1uM.6h"]],p.value = 0.05,n=NULL)$table$entrez,
                       cAQ_mBen.1uM.24h=topTags(res[["cAQ_mBen.1uM.24h"]],p.value = 0.05,n=NULL)$table$entrez,
                       cAQ_mBen.5uM.6h=topTags(res[["cAQ_mBen.5uM.6h"]],p.value = 0.05,n=NULL)$table$entrez,
                       cAQ_mBen.5uM.24h=topTags(res[["cAQ_mBen.5uM.24h"]],p.value = 0.05,n=NULL)$table$entrez)
                       
toptags.genes.id<-lapply(toptags.genes.id, function(x)x[!is.na(x)])

venn(toptags.genes.id[c(1,2)])
library(ggvenn)
ggvenn(
  toptags.genes.id, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  #  fill_color = c("#D55E0050","#0072B250", "#FF000050","#0000FF50"),
  stroke_size = 0.5, set_name_size = 4
)


##Gene Enrichment Analysis---------


comp.ego.MF <- compareCluster(toptags.genes.id, fun="enrichGO",OrgDb="org.Hs.eg.db",ont="MF",readable=TRUE)
head(comp.ego.MF)
write.table(comp.ego.MF,"./out/comp.ego.MF_table.csv",sep=",",row.names = F)
png("./out/comp.ego.MF_dotplot.png",width =1200 ,height = 1000)
dotplot(comp.ego.MF)
dev.off()

comp.ego.BP <- compareCluster(toptags.genes.id, fun="enrichGO",OrgDb="org.Hs.eg.db",ont="BP",readable=TRUE)
head(comp.ego.BP)
write.table(comp.ego.BP,"./out/comp.ego.BP_table.csv",sep=",",row.names = F)
png("./out/comp.ego.BP_dotplot.png",width =1200 ,height = 1000)
dotplot(comp.ego.BP)
dev.off()

comp.ego.CC <- compareCluster(toptags.genes.id, fun="enrichGO",OrgDb="org.Hs.eg.db",ont="CC",readable=TRUE)
head(comp.ego.CC)
write.table(comp.ego.CC,"./out/comp.ego.CC_table.csv",sep=",",row.names = F)
png("./out/comp.ego.CC_dotplot.png",width =1200 ,height = 1000)
dotplot(comp.ego.CC)
dev.off()

comp.ekegg <- compareCluster(toptags.genes.id, fun="enrichKEGG",organism="hsa", pvalueCutoff=0.05)
head(comp.ekegg)
# comp.ekegg.df<-as.data.frame(comp.ekegg)
# x<-strsplit(comp.ekegg.df$geneID,split = "/")
# xx<-lapply(x, function(a)bitr(a,fromType = "ENTREZID",toType = "SYMBOL",OrgDb="org.Hs.eg.db")$SYMBOL)
# comp.ekegg.df$geneID_SYMBOL<-lapply(xx, function(a)paste(a,collapse = "/")) %>% unlist
comp.ekegg<-setReadable(comp.ekegg, 'org.Hs.eg.db', 'ENTREZID')
write.table(comp.ekegg,"./out/comp.ekegg_table.csv",sep=",",row.names = F)
png("./out/comp.ekegg_dotplot.png",width =1200 ,height = 1000)
dotplot(comp.ekegg)
dev.off()
x2 <- pairwise_termsim(comp.ekegg)
#png("./out/comp.ekegg_emapplot.png",width =1000 ,height = 1000)
pdf("./out/comp.ekegg_emapplot.pdf",width =10 ,height = 10)
emapplot(x2,pie="Count",showCategory = 50)
dev.off()

comp.edo <- compareCluster(toptags.genes.id, fun="enrichDO",readable=TRUE)
head(comp.edo)
write.table(comp.edo,"./out/comp.edo_table.csv",sep=",",row.names = F)
png("./out/comp.edo_dotplot.png",width =1200 ,height = 1000)
dotplot(comp.edo)
dev.off()


##each data--------------------------------------
#cAQ_mBen.1uM.6h
tt<-topTags(res[["cAQ_mBen.1uM.6h"]],p.value = 0.05,n=NULL)
geneList.df<-data.frame(tt$table$entrez,tt$table$logFC)
geneList.df<-geneList.df[!duplicated(geneList.df$tt.table.entrez),]
geneList_cAQ_mBen.1uM.6h<-geneList.df$tt.table.logFC
names(geneList_cAQ_mBen.1uM.6h)<-geneList.df$tt.table.entrez

ekegg_cAQ_mBen.1uM.6h <- enrichKEGG(gene          = names(geneList_cAQ_mBen.1uM.6h),
                                    organism     = 'hsa',
                                    pvalueCutoff = 0.05)
ekegg_cAQ_mBen.1uM.6h
ekeggx_cAQ_mBen.1uM.6h<- setReadable(ekegg_cAQ_mBen.1uM.6h, 'org.Hs.eg.db', 'ENTREZID')
png("./out/ekeggx_cAQ_mBen.1uM.6h_cnetplot.png",width =1000 ,height = 1000)
cnetplot(ekeggx_cAQ_mBen.1uM.6h, categorySize="pvalue", foldChange=geneList_cAQ_mBen.1uM.6h)
dev.off()

#cAQ_mBen.1uM.24h
tt<-topTags(res[["cAQ_mBen.1uM.24h"]],p.value = 0.05,n=NULL)
geneList.df<-data.frame(tt$table$entrez,tt$table$logFC)
geneList.df<-geneList.df[!duplicated(geneList.df$tt.table.entrez),]
geneList_cAQ_mBen.1uM.24h<-geneList.df$tt.table.logFC
names(geneList_cAQ_mBen.1uM.24h)<-geneList.df$tt.table.entrez

ekegg_cAQ_mBen.1uM.24h <- enrichKEGG(gene          = names(geneList_cAQ_mBen.1uM.24h),
                                    organism     = 'hsa',
                                    pvalueCutoff = 0.05)
ekegg_cAQ_mBen.1uM.24h
ekeggx_cAQ_mBen.1uM.24h<- setReadable(ekegg_cAQ_mBen.1uM.24h, 'org.Hs.eg.db', 'ENTREZID')

png("./out/ekeggx_cAQ_mBen.1uM.24h_cnetplot.png",width =1000 ,height = 1000)
cnetplot(ekeggx_cAQ_mBen.1uM.24h, categorySize="pvalue", foldChange=geneList_cAQ_mBen.1uM.24h)
dev.off()

#cAQ_mBen.5uM.6h
tt<-topTags(res[["cAQ_mBen.5uM.6h"]],p.value = 0.05,n=NULL)
geneList.df<-data.frame(tt$table$entrez,tt$table$logFC)
geneList.df<-geneList.df[!duplicated(geneList.df$tt.table.entrez),]
geneList_cAQ_mBen.5uM.6h<-geneList.df$tt.table.logFC
names(geneList_cAQ_mBen.5uM.6h)<-geneList.df$tt.table.entrez

ekegg_cAQ_mBen.5uM.6h <- enrichKEGG(gene          = names(geneList_cAQ_mBen.5uM.6h),
                                    organism     = 'hsa',
                                    pvalueCutoff = 0.05)
ekegg_cAQ_mBen.5uM.6h
ekeggx_cAQ_mBen.5uM.6h<- setReadable(ekegg_cAQ_mBen.5uM.6h, 'org.Hs.eg.db', 'ENTREZID')

png("./out/ekeggx_cAQ_mBen.5uM.6h_cnetplot.png",width =1000 ,height = 1000)
cnetplot(ekeggx_cAQ_mBen.5uM.6h, categorySize="pvalue", foldChange=geneList_cAQ_mBen.5uM.6h)
dev.off()

#cAQ_mBen.5uM.24h
tt<-topTags(res[["cAQ_mBen.5uM.24h"]],p.value = 0.05,n=NULL)
geneList.df<-data.frame(tt$table$entrez,tt$table$logFC)
geneList.df<-geneList.df[!duplicated(geneList.df$tt.table.entrez),]
geneList_cAQ_mBen.5uM.24h<-geneList.df$tt.table.logFC
names(geneList_cAQ_mBen.5uM.24h)<-geneList.df$tt.table.entrez

ekegg_cAQ_mBen.5uM.24h <- enrichKEGG(gene          = names(geneList_cAQ_mBen.5uM.24h),
                                     organism     = 'hsa',
                                     pvalueCutoff = 0.05)
ekegg_cAQ_mBen.5uM.24h
ekeggx_cAQ_mBen.5uM.24h<- setReadable(ekegg_cAQ_mBen.5uM.24h, 'org.Hs.eg.db', 'ENTREZID')

png("./out/ekeggx_cAQ_mBen.5uM.24h_cnetplot.png",width =1000 ,height = 1000)
cnetplot(ekeggx_cAQ_mBen.5uM.24h, categorySize="pvalue", foldChange=geneList_cAQ_mBen.5uM.24h)
dev.off()


#Consider up-or down-regulation----------------------------------------

toptags.genes.id.updown<-list()
for(name in c("cAQ_mBen.1uM.6h","cAQ_mBen.1uM.24h","cAQ_mBen.5uM.6h","cAQ_mBen.5uM.24h")){
  tmp<-topTags(res[[name]],p.value = 0.05,n=NULL)$table
  toptags.genes.id.updown[[paste0(name,".up")]]<-tmp[tmp$logFC >=0,]$entrez
  toptags.genes.id.updown[[paste0(name,".down")]]<-tmp[tmp$logFC <0,]$entrez
}

toptags.genes.id<-lapply(toptags.genes.id, function(x)x[!is.na(x)])
#toptags.genes.id.updown<-lapply(toptags.genes.id.updown, as.numeric)

library(UpSetR)
library(ggvenn)

pdf("./out/DEG/deg_venn_fdr0.05_up.pdf",width = 10,height = 10)
upset(fromList(toptags.genes.id.updown[grep("up",names(toptags.genes.id.updown))]),order.by = "freq")
ggvenn(
  toptags.genes.id.updown[grep("up",names(toptags.genes.id.updown))], 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  #  fill_color = c("#D55E0050","#0072B250", "#FF000050","#0000FF50"),
  stroke_size = 0.5, set_name_size = 4
)
dev.off()

pdf("./out/DEG/deg_venn_fdr0.05_down.pdf",width = 10,height = 10)
upset(fromList(toptags.genes.id.updown[grep("down",names(toptags.genes.id.updown))]),order.by = "freq")
ggvenn(
  toptags.genes.id.updown[grep("down",names(toptags.genes.id.updown))], 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  #  fill_color = c("#D55E0050","#0072B250", "#FF000050","#0000FF50"),
  stroke_size = 0.5, set_name_size = 4
)
dev.off()


##Gene Enrichment Analysis---------

comp.updown.ego.MF <- compareCluster(toptags.genes.id.updown, fun="enrichGO",OrgDb="org.Hs.eg.db",ont="MF",readable=TRUE)
head(comp.updown.ego.MF)
write.table(comp.updown.ego.MF,"./out/comp.updown.ego.MF_table.csv",sep=",",row.names = F)
png("./out/comp.updown.ego.MF_dotplot.png",width =1500 ,height = 1000)
dotplot(comp.updown.ego.MF)
dev.off()

x2 <- pairwise_termsim(comp.updown.ego.MF)
pdf("./out/comp.updown.ego.MF_emapplot.pdf",width =10 ,height = 10)
emapplot(x2,pie="Count",showCategory = 30)+
  scale_fill_manual(values = c("#0072B250","#0072B2","#FF000050","#0000FF50","#FF0000","#0000FF")) #FDR=0.05
  #scale_fill_manual(values = c("#0072B250","#D55E00","#0072B2","#FF000050","#0000FF50","#FF0000","#0000FF")) #FDR=0.01
dev.off()

comp.updown.ego.BP <- compareCluster(toptags.genes.id.updown, fun="enrichGO",OrgDb="org.Hs.eg.db",ont="BP",readable=TRUE)
head(comp.updown.ego.BP)
write.table(comp.updown.ego.BP,"./out/comp.updown.ego.BP_table.csv",sep=",",row.names = F)
png("./out/comp.updown.ego.BP_dotplot.png",width =1800 ,height = 1000)
dotplot(comp.updown.ego.BP)
dev.off()

x2 <- pairwise_termsim(comp.updown.ego.BP)
pdf("./out/comp.updown.ego.BP_emapplot.pdf",width =10 ,height = 10)
emapplot(x2,pie="Count",showCategory = 50)+
  scale_fill_manual(values = c("#D55E0050","#0072B250","#D55E00","#0072B2","#FF000050","#0000FF50","#FF0000","#0000FF")) #FDR=0.05
  #scale_fill_manual(values = c("#D55E0050","#0072B250","#0072B2","#FF000050","#0000FF50","#FF0000","#0000FF")) #FDR=0.01
dev.off()

comp.updown.ego.CC <- compareCluster(toptags.genes.id.updown, fun="enrichGO",OrgDb="org.Hs.eg.db",ont="CC",readable=TRUE)
head(comp.updown.ego.CC)
write.table(comp.updown.ego.CC,"./out/comp.updown.ego.CC_table.csv",sep=",",row.names = F)
png("./out/comp.updown.ego.CC_dotplot.png",width =1500 ,height = 1000)
dotplot(comp.updown.ego.CC)
dev.off()

x2 <- pairwise_termsim(comp.updown.ego.CC)
pdf("./out/comp.updown.ego.CC_emapplot.pdf",width =10 ,height = 10)
emapplot(x2,pie="Count",showCategory = 50)+
  scale_fill_manual(values = c("#0072B250","#D55E00","#0072B2","#FF000050","#0000FF50","#FF0000","#0000FF")) #FDR=0.05,0.01
dev.off()

comp.updown.ekegg <- compareCluster(toptags.genes.id.updown, fun="enrichKEGG",organism="hsa", pvalueCutoff=0.05)

head(comp.updown.ekegg)
# comp.updown.ekegg.df<-as.data.frame(comp.updown.ekegg)
# x<-strsplit(comp.updown.ekegg.df$geneID,split = "/")
# xx<-lapply(x, function(a)bitr(a,fromType = "ENTREZID",toType = "SYMBOL",OrgDb="org.Hs.eg.db")$SYMBOL)
# comp.updown.ekegg.df$geneID_SYMBOL<-lapply(xx, function(a)paste(a,collapse = "/")) %>% unlist
comp.updown.ekegg<-setReadable(comp.updown.ekegg, 'org.Hs.eg.db', 'ENTREZID')
write.table(comp.updown.ekegg,"./out/comp.updown.ekegg_table.csv",sep=",",row.names = F)
png("./out/comp.updown.ekegg_dotplot.png",width =1500 ,height = 1000)
dotplot(comp.updown.ekegg)
dev.off()

x2 <- pairwise_termsim(comp.updown.ekegg)
#png("./out/comp.updown.ekegg_emapplot.png",width =1000 ,height = 1000)
pdf("./out/comp.updown.ekegg_emapplot.pdf",width =10 ,height = 10)
emapplot(x2,pie="Count",showCategory = 50)+
  scale_fill_manual(values = c("#0072B250","#D55E0050","#0072B2","#FF000050","#0000FF50","#FF0000","#0000FF")) #FDR=0.05
  #scale_fill_manual(values = c("#0072B250","#0072B2","#FF000050","#0000FF50","#FF0000","#0000FF")) #FDR=0.01
dev.off()

comp.updown.edo <- compareCluster(toptags.genes.id.updown, fun="enrichDO",readable=TRUE)
head(comp.updown.edo)
write.table(comp.updown.edo,"./out/comp.updown.edo_table.csv",sep=",",row.names = F)
png("./out/comp.updown.edo_dotplot.png",width =1200 ,height = 1000)
dotplot(comp.updown.edo)
dev.off()

x2 <- pairwise_termsim(comp.updown.edo)
pdf("./out/comp.updown.edo_emapplot.pdf",width =10 ,height = 10)
emapplot(x2,pie="Count",showCategory = 50)+
  scale_fill_manual(values = c("#FF000050","#FF0000","#0000FF")) #FDR=0.05
  #scale_fill_manual(values = c("#D55E0050","#0072B250","#D55E00","#FF000050","#FF0000")) #FDR=0.01
dev.off()


##each data-----------------------------
#cAQ_mBen.1uM.6h
tt<-topTags(res[["cAQ_mBen.1uM.6h"]],p.value = 0.05,n=NULL)
geneList.df<-data.frame(tt$table$entrez,tt$table$logFC)
geneList.df<-geneList.df[!duplicated(geneList.df$tt.table.entrez),]
geneList_cAQ_mBen.1uM.6h<-geneList.df$tt.table.logFC
names(geneList_cAQ_mBen.1uM.6h)<-geneList.df$tt.table.entrez

geneList_cAQ_mBen.1uM.6h.up<-geneList_cAQ_mBen.1uM.6h[geneList_cAQ_mBen.1uM.6h> 0]

ekegg_cAQ_mBen.1uM.6h.up <- enrichKEGG(gene          = names(geneList_cAQ_mBen.1uM.6h.up),
                                        organism     = 'hsa',
                                        pvalueCutoff = 0.05)
ekeggx_cAQ_mBen.1uM.6h.up<- setReadable(ekegg_cAQ_mBen.1uM.6h.up, 'org.Hs.eg.db', 'ENTREZID')
ekeggx_cAQ_mBen.1uM.6h.up
#なし
png("./out/ekeggx_cAQ_mBen.1uM.6h.up_cnetplot.png",width =1000 ,height = 1000)
cnetplot(ekeggx_cAQ_mBen.1uM.6h.up, categorySize="pvalue", foldChange=geneList_cAQ_mBen.1uM.6h.up)
dev.off()

geneList_cAQ_mBen.1uM.6h.down<-geneList_cAQ_mBen.1uM.6h[geneList_cAQ_mBen.1uM.6h< 0]

ekegg_cAQ_mBen.1uM.6h.down <- enrichKEGG(gene          = names(geneList_cAQ_mBen.1uM.6h.down),
                                          organism     = 'hsa',
                                          pvalueCutoff = 0.05)
ekeggx_cAQ_mBen.1uM.6h.down<- setReadable(ekegg_cAQ_mBen.1uM.6h.down, 'org.Hs.eg.db', 'ENTREZID')

png("./out/ekeggx_cAQ_mBen.1uM.6h.down_cnetplot.png",width =1000 ,height = 1000)
cnetplot(ekeggx_cAQ_mBen.1uM.6h.down, categorySize="pvalue", foldChange=geneList_cAQ_mBen.1uM.6h.down)
dev.off()

#cAQ_mBen.1uM.24h
tt<-topTags(res[["cAQ_mBen.1uM.24h"]],p.value = 0.05,n=NULL)
geneList.df<-data.frame(tt$table$entrez,tt$table$logFC)
geneList.df<-geneList.df[!duplicated(geneList.df$tt.table.entrez),]
geneList_cAQ_mBen.1uM.24h<-geneList.df$tt.table.logFC
names(geneList_cAQ_mBen.1uM.24h)<-geneList.df$tt.table.entrez

geneList_cAQ_mBen.1uM.24h.up<-geneList_cAQ_mBen.1uM.24h[geneList_cAQ_mBen.1uM.24h> 0]

ekegg_cAQ_mBen.1uM.24h.up <- enrichKEGG(gene          = names(geneList_cAQ_mBen.1uM.24h.up),
                                     organism     = 'hsa',
                                     pvalueCutoff = 0.05)
ekegg_cAQ_mBen.1uM.24h.up
ekeggx_cAQ_mBen.1uM.24h.up<- setReadable(ekegg_cAQ_mBen.1uM.24h.up, 'org.Hs.eg.db', 'ENTREZID')

png("./out/ekeggx_cAQ_mBen.1uM.24h.up_cnetplot.png",width =1000 ,height = 1000)
cnetplot(ekeggx_cAQ_mBen.1uM.24h.up, categorySize="pvalue", foldChange=geneList_cAQ_mBen.1uM.24h.up)
dev.off()

geneList_cAQ_mBen.1uM.24h.down<-geneList_cAQ_mBen.1uM.24h[geneList_cAQ_mBen.1uM.24h< 0]

ekegg_cAQ_mBen.1uM.24h.down <- enrichKEGG(gene          = names(geneList_cAQ_mBen.1uM.24h.down),
                                        organism     = 'hsa',
                                        pvalueCutoff = 0.05)
ekeggx_cAQ_mBen.1uM.24h.down<- setReadable(ekegg_cAQ_mBen.1uM.24h.down, 'org.Hs.eg.db', 'ENTREZID')

png("./out/ekeggx_cAQ_mBen.1uM.24h.down_cnetplot.png",width =1000 ,height = 1000)
cnetplot(ekeggx_cAQ_mBen.1uM.24h.down, categorySize="pvalue", foldChange=geneList_cAQ_mBen.1uM.24h.down)
dev.off()

#cAQ_mBen.5uM.6h
tt<-topTags(res[["cAQ_mBen.5uM.6h"]],p.value = 0.05,n=NULL)
geneList.df<-data.frame(tt$table$entrez,tt$table$logFC)
geneList.df<-geneList.df[!duplicated(geneList.df$tt.table.entrez),]
geneList_cAQ_mBen.5uM.6h<-geneList.df$tt.table.logFC
names(geneList_cAQ_mBen.5uM.6h)<-geneList.df$tt.table.entrez

geneList_cAQ_mBen.5uM.6h.up<-geneList_cAQ_mBen.5uM.6h[geneList_cAQ_mBen.5uM.6h> 0]

ekegg_cAQ_mBen.5uM.6h.up <- enrichKEGG(gene          = names(geneList_cAQ_mBen.5uM.6h.up),
                                       organism     = 'hsa',
                                       pvalueCutoff = 0.05)
ekeggx_cAQ_mBen.5uM.6h.up<- setReadable(ekegg_cAQ_mBen.5uM.6h.up, 'org.Hs.eg.db', 'ENTREZID')
ekeggx_cAQ_mBen.5uM.6h.up

png("./out/ekeggx_cAQ_mBen.5uM.6h.up_cnetplot.png",width =1000 ,height = 1000)
cnetplot(ekeggx_cAQ_mBen.5uM.6h.up, categorySize="pvalue", foldChange=geneList_cAQ_mBen.5uM.6h.up)
dev.off()

geneList_cAQ_mBen.5uM.6h.down<-geneList_cAQ_mBen.5uM.6h[geneList_cAQ_mBen.5uM.6h< 0]

ekegg_cAQ_mBen.5uM.6h.down <- enrichKEGG(gene          = names(geneList_cAQ_mBen.5uM.6h.down),
                                         organism     = 'hsa',
                                         pvalueCutoff = 0.05)
ekeggx_cAQ_mBen.5uM.6h.down<- setReadable(ekegg_cAQ_mBen.5uM.6h.down, 'org.Hs.eg.db', 'ENTREZID')

png("./out/ekeggx_cAQ_mBen.5uM.6h.down_cnetplot.png",width =1000 ,height = 1000)
cnetplot(ekeggx_cAQ_mBen.5uM.6h.down, categorySize="pvalue", foldChange=geneList_cAQ_mBen.5uM.6h.down)
dev.off()


#cAQ_mBen.5uM.24h
tt<-topTags(res[["cAQ_mBen.5uM.24h"]],p.value = 0.05,n=NULL)
geneList.df<-data.frame(tt$table$entrez,tt$table$logFC)
geneList.df<-geneList.df[!duplicated(geneList.df$tt.table.entrez),]
geneList_cAQ_mBen.5uM.24h<-geneList.df$tt.table.logFC
names(geneList_cAQ_mBen.5uM.24h)<-geneList.df$tt.table.entrez

geneList_cAQ_mBen.5uM.24h.up<-geneList_cAQ_mBen.5uM.24h[geneList_cAQ_mBen.5uM.24h> 0]

ekegg_cAQ_mBen.5uM.24h.up <- enrichKEGG(gene          = names(geneList_cAQ_mBen.5uM.24h.up),
                                        organism     = 'hsa',
                                        pvalueCutoff = 0.05)
ekegg_cAQ_mBen.5uM.24h.up
ekeggx_cAQ_mBen.5uM.24h.up<- setReadable(ekegg_cAQ_mBen.5uM.24h.up, 'org.Hs.eg.db', 'ENTREZID')

png("./out/ekeggx_cAQ_mBen.5uM.24h.up_cnetplot.png",width =1000 ,height = 1000)
cnetplot(ekeggx_cAQ_mBen.5uM.24h.up, categorySize="pvalue", foldChange=geneList_cAQ_mBen.5uM.24h.up)
dev.off()

geneList_cAQ_mBen.5uM.24h.down<-geneList_cAQ_mBen.5uM.24h[geneList_cAQ_mBen.5uM.24h< 0]

ekegg_cAQ_mBen.5uM.24h.down <- enrichKEGG(gene          = names(geneList_cAQ_mBen.5uM.24h.down),
                                          organism     = 'hsa',
                                          pvalueCutoff = 0.05)
ekeggx_cAQ_mBen.5uM.24h.down<- setReadable(ekegg_cAQ_mBen.5uM.24h.down, 'org.Hs.eg.db', 'ENTREZID')

png("./out/ekeggx_cAQ_mBen.5uM.24h.down_cnetplot.png",width =1000 ,height = 1000)
cnetplot(ekeggx_cAQ_mBen.5uM.24h.down, categorySize="pvalue", foldChange=geneList_cAQ_mBen.5uM.24h.down)
dev.off()



