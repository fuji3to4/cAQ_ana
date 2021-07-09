library(tidyverse)
library(edgeR)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(org.Hs.eg.db)
library(GenomicRanges)
library(GenomicFeatures)
library(ChIPseeker)
library(clusterProfiler)
load("cAQ_RNA-Seq_de.RData")


path<-"/home/fujii/work/G4Hunter/"
#chr<-"chr2"
g4<-data.frame()
for(chr in paste0("chr",c(1:22,"X","Y","M"))){

  tmp<-read.table(pipe(paste0("awk '{print $1,$2,$3,$5}' ",path,"Results_",chr,"/",chr,"-Merged.txt")),skip=2)
  colnames(tmp)<-c("start","end","sequece","score")
  tmp$chr<-chr

  g4<-rbind(g4,tmp)
}
#tmp<-read_tsv(paste0(path,"Results_",chr,"/",chr,"-Merged.txt"),
#             skip = 1,col_types = list(col_integer(),col_integer(),col_character(),col_integer(),col_double(),col_integer()))
#なぜか固まるHDDのせいか????
g4.gr<-GRanges(seqnames = g4$chr,IRanges(start = g4$start,end=g4$end),score=g4$score,sequence=DNAStringSet(g4$sequece))
rm(g4)
strand(g4.gr)<-ifelse(g4.gr$score>0,"+","-")
g4.gr$g4ID<-1:length(g4.gr)

tmp<-seqinfo(BSgenome.Hsapiens.NCBI.GRCh38)
tmp<-tmp[seqnames(tmp)[1:25]]
seqnames(tmp)<-paste0("chr",c(1:22,"X","Y","M"))
seqinfo(g4.gr)<-tmp

gencode.v29.gr<-import("~/QNAP2/GENCODEv29/gencode.v29.annotation.gtf")
seqinfo(gencode.v29.gr)<-tmp
gencode.v29.txdb<-makeTxDbFromGRanges(gencode.v29.gr)

#Test-------------------------------------------------
# #covplot(g4.gr,weightCol = "score",chrs=c("chr17"), xlim=c(4.5e7, 4.6e7))
# 
# promoter<-getPromoters(TxDb = gencode.v29.txdb,upstream = 3000,downstream = 3000)
# tagMatrix<-getTagMatrix(g4.gr,windows = promoter)
# # png("tagHeatmap_g4.png")
# # tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")
# # dev.off()
# 
# plotAvgProf(tagMatrix ,xlim=c(-3000,3000),
#             xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
# #            conf=0.95,resample=1000)
# 
# g4.peakAnno <- annotatePeak(g4.gr, tssRegion=c(-3000, 3000),level="gene",
#                          TxDb=gencode.v29.txdb, annoDb="org.Hs.eg.db")
# plotAnnoPie(g4.peakAnno)
# vennpie(g4.peakAnno)
# 
# #upsetplot(g4.peakAnno,vennpie=TRUE)
#
# promoter<-getPromoters(TxDb = gencode.v29.txdb,upstream = 3000,downstream = 3000,by = "gene")
# gene<-genes(gencode.v29.txdb)
# ol.pro<-findOverlaps(g4.gr,promoter)
# ol.gen<-findOverlaps(g4.gr,gene)
# assign(g4.gencode.v29.gr,g4.gr[c(from(ol.pro),from(ol.gen)) %>% unique])
# tagMatrixList[["g4.gencode.v29"]]<-getTagMatrix(g4.gencode.v29.gr,windows = promoter)
# peakAnnoList[["g4.gencode.v29"]]<-annotatePeak(g4.gencode.v29.gr, tssRegion=c(-3000, 3000),level="gene",
#                                    TxDb=gencode.v29.txdb, annoDb="org.Hs.eg.db")
#--------------------------------------------------------------

toptags.genes.ensembl<-list()
for(name in c("cAQ_mBen.1uM.6h","cAQ_mBen.1uM.24h","cAQ_mBen.5uM.6h","cAQ_mBen.5uM.24h")){
  tmp<-topTags(res[[name]],p.value = 0.05,n=NULL)$table
  toptags.genes.ensembl[[paste0(name,".up")]]<-tmp[tmp$logFC >=0,]$Ensembl
  toptags.genes.ensembl[[paste0(name,".down")]]<-tmp[tmp$logFC <0,]$Ensembl
}


for(name in names(toptags.genes.ensembl)){
  assign(paste0(name,".gr"),gencode.v29.gr[gencode.v29.gr$gene_id %in% toptags.genes.ensembl[[name]]])
  assign(paste0(name,".txdb"),makeTxDbFromGRanges(get(paste0(name,".gr"))))
}



# g4.grが大きいため計算に時間がかかってしまう。毎回回すのは非効率なので事前に絞っておく
# tagMatrixList<-list()
# peakAnnoList<-list()
# for(name in c("gencode.v29",names(toptags.genes.ensembl))){
#   promoter<-getPromoters(TxDb = get(paste0(name,".txdb")),upstream = 3000,downstream = 3000,by = "gene")
#   tagMatrixList[[name]]<-getTagMatrix(g4.gr,windows = promoter)
#   peakAnnoList[[name]]<-annotatePeak(g4.gr, tssRegion=c(-3000, 3000),TxDb=get(paste0(name,".txdb")),
#                                      level="gene",sameStrand = TRUE,overlap = "TSS",
#                                      genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream"))
# }




for(name in c("gencode.v29",names(toptags.genes.ensembl))){
  #name<-"cAQ_mBen.5uM.24h.down"
  promoter<-getPromoters(TxDb = get(paste0(name,".txdb")),upstream = 3000,downstream = 3000,by = "gene")
  gene<-genes(get(paste0(name,".txdb")))
  ol.pro<-findOverlaps(g4.gr,promoter)
  ol.gen<-findOverlaps(g4.gr,gene)
  assign(paste0("g4.",name,".gr"),g4.gr[c(from(ol.pro),from(ol.gen)) %>% unique])
}

tagMatrixList<-list()
peakAnnoList<-list()
for(name in c("gencode.v29",names(toptags.genes.ensembl))){
  promoter<-getPromoters(TxDb = get(paste0(name,".txdb")),upstream = 3000,downstream = 3000,by = "gene")
  tagMatrixList[[name]]<-getTagMatrix(get(paste0("g4.",name,".gr")),windows = promoter)
  peakAnnoList[[name]]<-annotatePeak(get(paste0("g4.",name,".gr")), tssRegion=c(-3000, 3000),TxDb=get(paste0(name,".txdb")),
                                     level="gene",sameStrand = TRUE,overlap = "all",
                                     genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream"))
}

anno.genes<-peakAnnoList[["gencode.v29"]]@anno$geneId %>% unique()
gencode.genes<-genes(gencode.v29.txdb)$gene_id

gencode.genes %in% anno.genes
# peakAnnoList.all<- peakAnnoList
# (peakAnnoList.all[["gencode.v29"]]@anno$annotation == peakAnnoList[["gencode.v29"]]@anno$annotation
# peakAnnoList[["gencode.v29"]]@anno %>% as.data.frame -> A
# View(A)
# peakAnnoList.all[["gencode.v29"]]@anno %>% as.data.frame ->B
# View(B)
# sum(A$annotation =="Promoter (<=1kb)",na.rm = T)
# sum(B$annotation =="Promoter (<=1kb)",na.rm = T)
# nrow(A)
# nrow(B)
#save(tagMatrixList,peakAnnoList,file="peakAnnoList_deg.rda")
#load(file="peakAnnoList_deg.rda")



pdf("./out/plot_genomicanno.pdf",width = 11,height = 8.5)
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList)
dev.off()

annobar.df<-data.frame(Category=c(),Feature=c(),Freq=c(),Count=c())
for(i in names(peakAnnoList)){
  annobar.df.tmp<-data.frame(
    Category=i,
    Feature=peakAnnoList[[i]]@annoStat[,1],
    Freq=peakAnnoList[[i]]@annoStat[,2],
    Count=round(peakAnnoList[[i]]@annoStat[,2] /100 * peakAnnoList[[i]]@peakNum)
  )
  annobar.df<-rbind(annobar.df,annobar.df.tmp)
}

p.annobar <- ggplot(annobar.df[annobar.df$Category !="gencode.v29",],aes(x = Category,y = Count,fill = Feature))
p.annobar1 <- p.annobar + geom_bar(stat = "identity") +
  scale_x_discrete(limits=names(toptags.genes.ensembl)) +
  #scale_fill_discrete(limits=tb.nb$labels) +
  scale_fill_brewer(palette = "Paired") +
  coord_flip() +theme_bw()
p.annobar1
# ggsave(file = "./out/plot_AnnoBar_genomicanno.pdf", plot = p.annobar1, width = 11, height = 8.5)
# 
# pdf("./out/plotAvgProf_comp_PQS.pdf",width = 11,height = 8.5)
# plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
# dev.off()
# pdf("./out/plotAvgProf_comp_PQS_row.pdf",width = 10,height = 20)
# plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")
# dev.off()




