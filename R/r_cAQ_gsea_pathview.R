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
                       
toptags.genes.id<-lapply(toptags.genes.id, as.numeric)
#toptags.genes.id<-lapply(toptags.genes.id,unique)

venn(toptags.genes.id[c(1,2)])


#--------------------------------------------------
#Gene Enrichment Analysis
#---------------------------------------------------
comp.ego.MF <- compareCluster(toptags.genes.id, fun="enrichGO",OrgDb="org.Hs.eg.db",ont="MF",readable=TRUE)
head(comp.ego.MF)
#write.table(comp.ego.MF,"./out/comp.ego.MF_table.csv",sep=",",row.names = F)
#png("./out/comp.ego.MF_dotplot.png",width =1200 ,height = 1000)
dotplot(comp.ego.MF)
#dev.off()

comp.ego.BP <- compareCluster(toptags.genes.id, fun="enrichGO",OrgDb="org.Hs.eg.db",ont="BP",readable=TRUE)
head(comp.ego.BP)
#write.table(comp.ego.BP,"./out/comp.ego.BP_table.csv",sep=",",row.names = F)
#png("./out/comp.ego.BP_dotplot.png",width =1200 ,height = 1000)
dotplot(comp.ego.BP)
#dev.off()

comp.ego.CC <- compareCluster(toptags.genes.id, fun="enrichGO",OrgDb="org.Hs.eg.db",ont="CC",readable=TRUE)
head(comp.ego.CC)
#write.table(comp.ego.CC,"./out/comp.ego.CC_table.csv",sep=",",row.names = F)
#png("./out/comp.ego.CC_dotplot.png",width =1200 ,height = 1000)
dotplot(comp.ego.CC)
#dev.off()

comp.ekegg <- compareCluster(toptags.genes.id, fun="enrichKEGG",organism="hsa", pvalueCutoff=0.05)
head(comp.ekegg)
#write.table(comp.ekegg,"./out/comp.ekegg_table.csv",sep=",",row.names = F)
#png("./out/comp.ekegg_dotplot.png",width =1200 ,height = 1000)
dotplot(comp.ekegg)
#dev.off()
x2 <- pairwise_termsim(comp.ekegg)
#png("./out/comp.ekegg_emapplot.png",width =1000 ,height = 1000)
emapplot(x2)
#dev.off()

comp.edo <- compareCluster(toptags.genes.id, fun="enrichDO",readable=TRUE)
head(comp.edo)
#write.table(comp.edo,"./out/comp.edo_table.csv",sep=",",row.names = F)
#png("./out/comp.edo_dotplot.png",width =1200 ,height = 1000)
dotplot(comp.edo)
#dev.off()


##each data
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
#png("./out/ekeggx_cAQ_mBen.1uM.6h_cnetplot.png",width =1000 ,height = 1000)
cnetplot(ekeggx_cAQ_mBen.1uM.6h, categorySize="pvalue", foldChange=geneList_cAQ_mBen.1uM.6h)
#dev.off()

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

#png("./out/ekeggx_cAQ_mBen.1uM.24h_cnetplot.png",width =1000 ,height = 1000)
cnetplot(ekeggx_cAQ_mBen.1uM.24h, categorySize="pvalue", foldChange=geneList_cAQ_mBen.1uM.24h)
#dev.off()

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

#png("./out/ekeggx_cAQ_mBen.5uM.6h_cnetplot.png",width =1000 ,height = 1000)
cnetplot(ekeggx_cAQ_mBen.5uM.6h, categorySize="pvalue", foldChange=geneList_cAQ_mBen.5uM.6h)
#dev.off()

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

#png("./out/ekeggx_cAQ_mBen.5uM.24h_cnetplot.png",width =1000 ,height = 1000)
cnetplot(ekeggx_cAQ_mBen.5uM.24h, categorySize="pvalue", foldChange=geneList_cAQ_mBen.5uM.24h)
#dev.off()

library(pathview)
#Spliceosome
browseKEGG(ekegg_cAQ_mBen.5uM.6h,"hsa03040")
hsa03040 <- pathview(gene.data  = geneList_cAQ_mBen.5uM.6h,
                     pathway.id = "hsa03040",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList_cAQ_mBen.5uM.6h)), cpd=1))

#mTOR signaling pathway
browseKEGG(ekegg_cAQ_mBen.5uM.6h,"hsa04150")
hsa04150 <- pathview(gene.data  = geneList_cAQ_mBen.5uM.6h,
                     pathway.id = "	hsa04150",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList_cAQ_mBen.5uM.6h)), cpd=1))

#AMPK signaling pathway
browseKEGG(ekegg_cAQ_mBen.5uM.6h,"hsa04152")
hsa04152 <- pathview(gene.data  = geneList_cAQ_mBen.5uM.6h,
                     pathway.id = "hsa04152",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList_cAQ_mBen.5uM.6h)), cpd=1))


#AMPK signaling pathway
browseKEGG(ekegg_cAQ_mBen.5uM.6h,"hsa04310")
hsa04310  <- pathview(gene.data  = geneList_cAQ_mBen.5uM.6h,
                     pathway.id = "hsa04310",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList_cAQ_mBen.5uM.6h)), cpd=1))

browseKEGG(ekegg_cAQ_mBen.5uM.6h,"hsa05410")
hsa05410   <- pathview(gene.data  = geneList_cAQ_mBen.5uM.6h,
                      pathway.id = "hsa05410",
                      species    = "hsa",
                      limit      = list(gene=max(abs(geneList_cAQ_mBen.5uM.6h)), cpd=1))
#enrichDO
edo_cAQ_mBen.5uM.6h <- enrichDO(gene          = names(geneList_cAQ_mBen.5uM.6h),
              ont           = "DO",
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = TRUE)

head(edo_cAQ_mBen.5uM.6h)
cnetplot(edo_cAQ_mBen.5uM.6h, categorySize="pvalue", foldChange=geneList_cAQ_mBen.5uM.24h)
