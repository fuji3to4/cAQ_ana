library(patchwork)
library(ggridges)
#source("r_cAQ_PQS_ChIPseeker.R")


#-1k~+1kbp-------------------------------------------
data<-peakAnnoList[["gencode.v29"]]@anno %>% as.data.frame()
ids.u1kd1k<-data$distanceToTSS <=1000 & data$distanceToTSS >=-1000 #結局Promoter(<=1kb)と同じ
g4.count.u1kd1k<-data[ids.u1kd1k,"geneId"] %>% table

gencode.genes<-genes(gencode.v29.txdb)$gene_id
notin.genes<-gencode.genes[!(gencode.genes %in% names(g4.count.u1kd1k))]
x<-rep(0,length(notin.genes))
names(x)<-notin.genes
g4.count.u1kd1k<-c(g4.count.u1kd1k,x)

g4.count.u1kd1k.df.list<-list("all"=g4.count.u1kd1k)
g4.count.u1kd1k.df<-data.frame(Category="all",Gene=names(g4.count.u1kd1k),Count=as.numeric(g4.count.u1kd1k))
for(name in names(toptags.genes.ensembl)){
  ids<-names(g4.count.u1kd1k) %in% toptags.genes.ensembl[[name]]
  tmp.count<-g4.count.u1kd1k[ids]
  g4.count.u1kd1k.df.list[[name]]<-tmp.count
  tmp.df<-data.frame(Category=name,Gene=names(tmp.count),Count=as.numeric(tmp.count))
  g4.count.u1kd1k.df<-rbind(g4.count.u1kd1k.df,tmp.df)
}

boxplot(log2(Count+1)~Category,data=g4.count.u1kd1k.df)

sapply(g4.count.u1kd1k.df.list, median)
sapply(g4.count.u1kd1k.df.list, mean)

library(ggridges)

p.riges <- ggplot(g4.count.u1kd1k.df,aes(x=Count,y=Category,fill=Category)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2, vline_color = "red",alpha = .7,to=30)+
  scale_x_continuous(name="PQS Counts") +
  scale_y_discrete(limits=c("all",names(toptags.genes.ensembl))) +
  scale_fill_cyclical(values =c("gray",rep(c("#8080ff","#ff8080"),4)))+
  theme_bw()

p.riges
ggsave(file = "./out/geom_density_PQScounts_u1kd1k.pdf", plot = p.riges, width = 11, height = 8.5)



sapply(g4.count.u1kd1k.df.list,length)
g4.count.u1kd1k.dummy<-data.frame(gene=names(g4.count.u1kd1k.df.list[["all"]]),Count=as.numeric(g4.count.u1kd1k.df.list[["all"]]))

for(name in names(toptags.genes.ensembl)){
  g4.count.u1kd1k.dummy[,name]<-g4.count.u1kd1k.dummy$gene %in% names(g4.count.u1kd1k.df.list[[name]]) %>% as.integer() %>% as.factor() 
}

glm.u1kd1k.g4<-list()
for(name in names(toptags.genes.ensembl)){
  glm.u1kd1k.g4[[name]] <- glm(formula = as.formula(paste0("Count~",name)),data=g4.count.u1kd1k.dummy,family = poisson())
}
lapply(glm.u1kd1k.g4, summary)



wilcox.u1kd1k.g4<-list()
for(name in names(toptags.genes.ensembl)){
  wilcox.u1kd1k.g4[[name]] <- wilcox.test(formula = as.formula(paste0("Count~",name)),data=g4.count.u1kd1k.dummy,family = poisson())
}
wilcox.u1kd1k.g4


g4.count.u1kd1k.reg<-data.frame(Category=NULL,gene=NULL,Count=NULL,reg=NULL)
for(name in names(toptags.genes.ensembl)){
  tmp<-data.frame(Category=name,gene=g4.count.u1kd1k.dummy$gene,Count=g4.count.u1kd1k.dummy$Count,reg=g4.count.u1kd1k.dummy[,name])
  g4.count.u1kd1k.reg<-rbind(g4.count.u1kd1k.reg,tmp)
}

#write.csv(g4.count.u1kd1k.reg,file="./out/g4.count.u1kd1k.reg.csv",row.names = F,quote = F)

#Promoter (<=1kb)-------------------------------------------
# data<-peakAnnoList[["gencode.v29"]]@anno %>% as.data.frame()
# ids.prom1k<-data$annotation =="Promoter (<=1kb)"
# g4.count.prom1k<-data[ids.prom1k,"geneId"] %>% table
# 
# gencode.genes<-genes(gencode.v29.txdb)$gene_id
# notin.genes<-gencode.genes[!(gencode.genes %in% names(g4.count.prom1k))]
# x<-rep(0,length(notin.genes))
# names(x)<-notin.genes
# g4.count.prom1k<-c(g4.count.prom1k,x)
# 
# g4.count.prom1k.df.list<-list("all"=g4.count.prom1k)
# g4.count.prom1k.df<-data.frame(Category="all",Gene=names(g4.count.prom1k),Count=as.numeric(g4.count.prom1k))
# for(name in names(toptags.genes.ensembl)){
#   ids<-names(g4.count.prom1k) %in% toptags.genes.ensembl[[name]]
#   tmp.count<-g4.count.prom1k[ids]
#   g4.count.prom1k.df.list[[name]]<-tmp.count
#   tmp.df<-data.frame(Category=name,Gene=names(tmp.count),Count=as.numeric(tmp.count))
#   g4.count.prom1k.df<-rbind(g4.count.prom1k.df,tmp.df)
# }
# #g4.cout.df$Category<-as.factor(g4.cout.df$Category)
# boxplot(log2(Count+1)~Category,data=g4.count.prom1k.df)
# 
# sapply(g4.count.prom1k.df.list, median)
# sapply(g4.count.prom1k.df.list, mean)
# 
# library(ggridges)
# 
# p.riges <- ggplot(g4.count.prom1k.df,aes(x=Count,y=Category,fill=Category)) +
#   geom_density_ridges(quantile_lines = TRUE, quantiles = 2, vline_color = "red",alpha = .7,to=30)+
#   #geom_density_ridges(stat="binline",alpha = .7,bins=30)+
#   scale_x_continuous(name="PQS Counts") +
#   scale_y_discrete(limits=c("all",names(toptags.genes.ensembl))) +
#   scale_fill_cyclical(values =c("gray",rep(c("#8080ff","#ff8080"),4)))+
#   theme_bw()
# 
# p.riges
# ggsave(file = "./out/geom_density_PQScounts_prom1k.pdf", plot = p.riges, width = 11, height = 8.5)
# 
# #t.test(g4.cout.list[["all"]],g4.cout.list[["cAQ_mBen.1uM.6h.up"]])
# 
# sapply(g4.count.prom1k.df.list,length)
# g4.count.prom1k.dummy<-data.frame(gene=names(g4.count.prom1k.df.list[["all"]]),Count=as.numeric(g4.count.prom1k.df.list[["all"]]))
# #name<-"cAQ_mBen.1uM.6h.up"
# for(name in names(toptags.genes.ensembl)){
#   g4.count.prom1k.dummy[,name]<-g4.count.prom1k.dummy$gene %in% names(g4.count.prom1k.df.list[[name]]) %>% as.integer() %>% as.factor() 
# }
# 
# glm.prom1k.g4<-list()
# #name<-"cAQ_mBen.5uM.6h.down"
# for(name in names(toptags.genes.ensembl)){
#   glm.prom1k.g4[[name]] <- glm(formula = as.formula(paste0("Count~",name)),data=g4.count.prom1k.dummy,family = poisson())
# }
# lapply(glm.prom1k.g4, summary)
# 
# 
# #0~+1kbp-------------------------------------------
# data<-peakAnnoList[["gencode.v29"]]@anno %>% as.data.frame()
# ids.d1k<-data$distanceToTSS <=1000 & data$distanceToTSS >=0 
# g4.count.d1k<-data[ids.d1k,"geneId"] %>% table
# 
# gencode.genes<-genes(gencode.v29.txdb)$gene_id
# notin.genes<-gencode.genes[!(gencode.genes %in% names(g4.count.d1k))]
# x<-rep(0,length(notin.genes))
# names(x)<-notin.genes
# g4.count.d1k<-c(g4.count.d1k,x)
# 
# g4.count.d1k.df.list<-list("all"=g4.count.d1k)
# g4.count.d1k.df<-data.frame(Category="all",Gene=names(g4.count.d1k),Count=as.numeric(g4.count.d1k))
# for(name in names(toptags.genes.ensembl)){
#   ids<-names(g4.count.d1k) %in% toptags.genes.ensembl[[name]]
#   tmp.count<-g4.count.d1k[ids]
#   g4.count.d1k.df.list[[name]]<-tmp.count
#   tmp.df<-data.frame(Category=name,Gene=names(tmp.count),Count=as.numeric(tmp.count))
#   g4.count.d1k.df<-rbind(g4.count.d1k.df,tmp.df)
# }
# 
# boxplot(log2(Count+1)~Category,data=g4.count.d1k.df)
# 
# sapply(g4.count.d1k.df.list, median)
# sapply(g4.count.d1k.df.list, mean)
# 
# library(ggridges)
# 
# p.riges <- ggplot(g4.count.d1k.df,aes(x=Count,y=Category,fill=Category)) +
#   geom_density_ridges(quantile_lines = TRUE, quantiles = 2, vline_color = "red",alpha = .7,to=30)+
#   scale_x_continuous(name="PQS Counts") +
#   scale_y_discrete(limits=c("all",names(toptags.genes.ensembl))) +
#   scale_fill_cyclical(values =c("gray",rep(c("#8080ff","#ff8080"),4)))+
#   theme_bw()
# 
# p.riges
# ggsave(file = "./out/geom_density_PQScounts_d1k.pdf", plot = p.riges, width = 11, height = 8.5)
# 
# 
# 
# sapply(g4.count.d1k.df.list,length)
# g4.count.d1k.dummy<-data.frame(gene=names(g4.count.d1k.df.list[["all"]]),Count=as.numeric(g4.count.d1k.df.list[["all"]]))
# 
# for(name in names(toptags.genes.ensembl)){
#   g4.count.d1k.dummy[,name]<-g4.count.d1k.dummy$gene %in% names(g4.count.d1k.df.list[[name]]) %>% as.integer() %>% as.factor() 
# }
# 
# glm.d1k.g4<-list()
# for(name in names(toptags.genes.ensembl)){
#   glm.d1k.g4[[name]] <- glm(formula = as.formula(paste0("Count~",name)),data=g4.count.d1k.dummy,family = poisson())
# }
# lapply(glm.d1k.g4, summary)
# 
# 
# 
# g4.count.d1k.reg<-data.frame(Category=NULL,gene=NULL,Count=NULL,reg=NULL)
# for(name in names(toptags.genes.ensembl)){
#   tmp<-data.frame(Category=name,gene=g4.count.d1k.dummy$gene,Count=g4.count.d1k.dummy$Count,reg=g4.count.d1k.dummy[,name])
#   g4.count.d1k.reg<-rbind(g4.count.d1k.reg,tmp)
# }
# 
# #write.csv(g4.count.d1k.reg,file="./out/g4.count.d1k.reg.csv",row.names = F,quote = F)
# 
# #-1k~0bp-------------------------------------------
# data<-peakAnnoList[["gencode.v29"]]@anno %>% as.data.frame()
# ids.u1k<-data$distanceToTSS <=0 & data$distanceToTSS >=-1000 #結局Promoter(<=1kb)と同じ
# g4.count.u1k<-data[ids.u1k,"geneId"] %>% table
# 
# gencode.genes<-genes(gencode.v29.txdb)$gene_id
# notin.genes<-gencode.genes[!(gencode.genes %in% names(g4.count.u1k))]
# x<-rep(0,length(notin.genes))
# names(x)<-notin.genes
# g4.count.u1k<-c(g4.count.u1k,x)
# 
# g4.count.u1k.df.list<-list("all"=g4.count.u1k)
# g4.count.u1k.df<-data.frame(Category="all",Gene=names(g4.count.u1k),Count=as.numeric(g4.count.u1k))
# for(name in names(toptags.genes.ensembl)){
#   ids<-names(g4.count.u1k) %in% toptags.genes.ensembl[[name]]
#   tmp.count<-g4.count.u1k[ids]
#   g4.count.u1k.df.list[[name]]<-tmp.count
#   tmp.df<-data.frame(Category=name,Gene=names(tmp.count),Count=as.numeric(tmp.count))
#   g4.count.u1k.df<-rbind(g4.count.u1k.df,tmp.df)
# }
# 
# boxplot(log2(Count+1)~Category,data=g4.count.u1k.df)
# 
# sapply(g4.count.u1k.df.list, median)
# sapply(g4.count.u1k.df.list, mean)
# 
# library(ggridges)
# 
# p.riges <- ggplot(g4.count.u1k.df,aes(x=Count,y=Category,fill=Category)) +
#   geom_density_ridges(quantile_lines = TRUE, quantiles = 2, vline_color = "red",alpha = .7,to=30)+
#   scale_x_continuous(name="PQS Counts") +
#   scale_y_discrete(limits=c("all",names(toptags.genes.ensembl))) +
#   scale_fill_cyclical(values =c("gray",rep(c("#8080ff","#ff8080"),4)))+
#   theme_bw()
# 
# p.riges
# ggsave(file = "./out/geom_density_PQScounts_u1k.pdf", plot = p.riges, width = 11, height = 8.5)
# 
# 
# 
# sapply(g4.count.u1k.df.list,length)
# g4.count.u1k.dummy<-data.frame(gene=names(g4.count.u1k.df.list[["all"]]),Count=as.numeric(g4.count.u1k.df.list[["all"]]))
# 
# for(name in names(toptags.genes.ensembl)){
#   g4.count.u1k.dummy[,name]<-g4.count.u1k.dummy$gene %in% names(g4.count.u1k.df.list[[name]]) %>% as.integer() %>% as.factor() 
# }
# 
# glm.u1k.g4<-list()
# for(name in names(toptags.genes.ensembl)){
#   glm.u1k.g4[[name]] <- glm(formula = as.formula(paste0("Count~",name)),data=g4.count.u1k.dummy,family = poisson())
# }
# lapply(glm.u1k.g4, summary)
# 
# 
# 
# g4.count.u1k.reg<-data.frame(Category=NULL,gene=NULL,Count=NULL,reg=NULL)
# for(name in names(toptags.genes.ensembl)){
#   tmp<-data.frame(Category=name,gene=g4.count.u1k.dummy$gene,Count=g4.count.u1k.dummy$Count,reg=g4.count.u1k.dummy[,name])
#   g4.count.u1k.reg<-rbind(g4.count.u1k.reg,tmp)
# }
# 
# #write.csv(g4.count.u1k.reg,file="./out/g4.count.u1k.reg.csv",row.names = F,quote = F)
# 
# 
# 
# #-2k~+100bp imitate S.Neidle's mothods-------------------------------------------
# data<-peakAnnoList[["gencode.v29"]]@anno %>% as.data.frame()
# ids<-data$distanceToTSS <=100 & data$distanceToTSS >=-2000
# g4.count.u2kd100<-data[ids,"geneId"] %>% table
# 
# gencode.genes<-genes(gencode.v29.txdb)$gene_id
# notin.genes<-gencode.genes[!(gencode.genes %in% names(g4.count.u2kd100))]
# x<-rep(0,length(notin.genes))
# names(x)<-notin.genes
# g4.count.u2kd100<-c(g4.count.u2kd100,x)
# 
# g4.count.u2kd100.df.list<-list("all"=g4.count.u2kd100)
# g4.count.u2kd100.df<-data.frame(Category="all",Gene=names(g4.count.u2kd100),Count=as.numeric(g4.count.u2kd100))
# for(name in names(toptags.genes.ensembl)){
#   ids<-names(g4.count.u2kd100) %in% toptags.genes.ensembl[[name]]
#   tmp.count<-g4.count.u2kd100[ids]
#   g4.count.u2kd100.df.list[[name]]<-tmp.count
#   tmp.df<-data.frame(Category=name,Gene=names(tmp.count),Count=as.numeric(tmp.count))
#   g4.count.u2kd100.df<-rbind(g4.count.u2kd100.df,tmp.df)
# }
# #g4.cout.df$Category<-as.factor(g4.cout.df$Category)
# boxplot(log10(Count+1)~Category,data=g4.count.u2kd100.df)
# 
# sapply(g4.count.u2kd100.df.list, median)
# sapply(g4.count.u2kd100.df.list, mean)
# 
# library(ggridges)
# 
# p.riges <- ggplot(g4.count.u2kd100.df,aes(x=Count,y=Category,fill=Category)) +
#   geom_density_ridges(quantile_lines = TRUE, quantiles = 2, vline_color = "red",alpha = .7,to=20)+
#   scale_x_continuous(name="PQS Counts") +
#   scale_y_discrete(limits=c("all",names(toptags.genes.ensembl))) +
#   scale_fill_cyclical(values =c("gray",rep(c("#8080ff","#ff8080"),4)))+
#   theme_bw()
# 
# p.riges
# ggsave(file = "./out/geom_density_PQScounts_u2kd100.pdf", plot = p.riges, width = 11, height = 8.5)
# 
# #t.test(g4.cout.list[["all"]],g4.cout.list[["cAQ_mBen.1uM.6h.up"]])
# 
# sapply(g4.count.u2kd100.df.list,length)
# g4.count.u2kd100.dummy<-data.frame(gene=names(g4.count.u2kd100.df.list[["all"]]),Count=as.numeric(g4.count.u2kd100.df.list[["all"]]))
# #name<-"cAQ_mBen.1uM.6h.up"
# for(name in names(toptags.genes.ensembl)){
#   g4.count.u2kd100.dummy[,name]<-g4.count.u2kd100.dummy$gene %in% names(g4.count.u2kd100.df.list[[name]]) %>% as.integer() %>% as.factor() 
# }
# 
# glm.u2kd100.g4<-list()
# #name<-"cAQ_mBen.5uM.6h.down"
# for(name in names(toptags.genes.ensembl)){
#   glm.u2kd100.g4[[name]] <- glm(formula = as.formula(paste0("Count~",name)),data=g4.count.u2kd100.dummy,family = poisson())
# }
# lapply(glm.u2kd100.g4, summary)
# 
# 
# #0~+100bp-------------------------------------------
# data<-peakAnnoList[["gencode.v29"]]@anno %>% as.data.frame()
# ids.d100<-data$distanceToTSS <=100 & data$distanceToTSS >=0 
# g4.count.d100<-data[ids.d100,"geneId"] %>% table
# 
# gencode.genes<-genes(gencode.v29.txdb)$gene_id
# notin.genes<-gencode.genes[!(gencode.genes %in% names(g4.count.d100))]
# x<-rep(0,length(notin.genes))
# names(x)<-notin.genes
# g4.count.d100<-c(g4.count.d100,x)
# 
# g4.count.d100.df.list<-list("all"=g4.count.d100)
# g4.count.d100.df<-data.frame(Category="all",Gene=names(g4.count.d100),Count=as.numeric(g4.count.d100))
# for(name in names(toptags.genes.ensembl)){
#   ids<-names(g4.count.d100) %in% toptags.genes.ensembl[[name]]
#   tmp.count<-g4.count.d100[ids]
#   g4.count.d100.df.list[[name]]<-tmp.count
#   tmp.df<-data.frame(Category=name,Gene=names(tmp.count),Count=as.numeric(tmp.count))
#   g4.count.d100.df<-rbind(g4.count.d100.df,tmp.df)
# }
# 
# boxplot(log10(Count+1)~Category,data=g4.count.d100.df)
# 
# sapply(g4.count.d100.df.list, median)
# sapply(g4.count.d100.df.list, mean)
# 
# 
# p.riges <- ggplot(g4.count.d100.df,aes(x=Count,y=Category,fill=Category)) +
#   geom_density_ridges(quantile_lines = TRUE, quantiles = 2, vline_color = "red",alpha = .7,to=30)+
#   scale_x_continuous(name="PQS Counts") +
#   scale_y_discrete(limits=c("all",names(toptags.genes.ensembl))) +
#   scale_fill_cyclical(values =c("gray",rep(c("#8080ff","#ff8080"),4)))+
#   theme_bw()
# 
# p.riges
# ggsave(file = "./out/geom_density_PQScounts_d100.pdf", plot = p.riges, width = 11, height = 8.5)
# 
# 
# 
# sapply(g4.count.d100.df.list,length)
# g4.count.d100.dummy<-data.frame(gene=names(g4.count.d100.df.list[["all"]]),Count=as.numeric(g4.count.d100.df.list[["all"]]))
# 
# for(name in names(toptags.genes.ensembl)){
#   g4.count.d100.dummy[,name]<-g4.count.d100.dummy$gene %in% names(g4.count.d100.df.list[[name]]) %>% as.integer() %>% as.factor() 
# }
# 
# glm.d100.g4<-list()
# for(name in names(toptags.genes.ensembl)){
#   glm.d100.g4[[name]] <- glm(formula = as.formula(paste0("Count~",name)),data=g4.count.d100.dummy,family = poisson())
# }
# lapply(glm.d100.g4, summary)
# 
# 
# 
# g4.count.d100.reg<-data.frame(Category=NULL,gene=NULL,Count=NULL,reg=NULL)
# for(name in names(toptags.genes.ensembl)){
#   tmp<-data.frame(Category=name,gene=g4.count.d100.dummy$gene,Count=g4.count.d100.dummy$Count,reg=g4.count.d100.dummy[,name])
#   g4.count.d100.reg<-rbind(g4.count.d100.reg,tmp)
# }
# 
# #write.csv(g4.count.d100.reg,file="./out/g4.count.d100.reg.csv",row.names = F,quote = F)
# 
# #-100~0bp-------------------------------------------
# data<-peakAnnoList[["gencode.v29"]]@anno %>% as.data.frame()
# ids.u100<-data$distanceToTSS <=0 & data$distanceToTSS >=-100
# g4.count.u100<-data[ids.u100,"geneId"] %>% table
# 
# gencode.genes<-genes(gencode.v29.txdb)$gene_id
# notin.genes<-gencode.genes[!(gencode.genes %in% names(g4.count.u100))]
# x<-rep(0,length(notin.genes))
# names(x)<-notin.genes
# g4.count.u100<-c(g4.count.u100,x)
# 
# g4.count.u100.df.list<-list("all"=g4.count.u100)
# g4.count.u100.df<-data.frame(Category="all",Gene=names(g4.count.u100),Count=as.numeric(g4.count.u100))
# for(name in names(toptags.genes.ensembl)){
#   ids<-names(g4.count.u100) %in% toptags.genes.ensembl[[name]]
#   tmp.count<-g4.count.u100[ids]
#   g4.count.u100.df.list[[name]]<-tmp.count
#   tmp.df<-data.frame(Category=name,Gene=names(tmp.count),Count=as.numeric(tmp.count))
#   g4.count.u100.df<-rbind(g4.count.u100.df,tmp.df)
# }
# 
# boxplot(log2(Count+1)~Category,data=g4.count.u100.df)
# 
# sapply(g4.count.u100.df.list, median)
# sapply(g4.count.u100.df.list, mean)
# 
# library(ggridges)
# 
# p.riges <- ggplot(g4.count.u100.df,aes(x=Count,y=Category,fill=Category)) +
#   geom_density_ridges(quantile_lines = TRUE, quantiles = 2, vline_color = "red",alpha = .7,to=30)+
#   scale_x_continuous(name="PQS Counts") +
#   scale_y_discrete(limits=c("all",names(toptags.genes.ensembl))) +
#   scale_fill_cyclical(values =c("gray",rep(c("#8080ff","#ff8080"),4)))+
#   theme_bw()
# 
# p.riges
# ggsave(file = "./out/geom_density_PQScounts_u100.pdf", plot = p.riges, width = 11, height = 8.5)
# 
# 
# 
# sapply(g4.count.u100.df.list,length)
# g4.count.u100.dummy<-data.frame(gene=names(g4.count.u100.df.list[["all"]]),Count=as.numeric(g4.count.u100.df.list[["all"]]))
# 
# for(name in names(toptags.genes.ensembl)){
#   g4.count.u100.dummy[,name]<-g4.count.u100.dummy$gene %in% names(g4.count.u100.df.list[[name]]) %>% as.integer() %>% as.factor() 
# }
# 
# glm.u100.g4<-list()
# for(name in names(toptags.genes.ensembl)){
#   glm.u100.g4[[name]] <- glm(formula = as.formula(paste0("Count~",name)),data=g4.count.u100.dummy,family = poisson())
# }
# lapply(glm.u100.g4, summary)
# 
# 
# 
# g4.count.u100.reg<-data.frame(Category=NULL,gene=NULL,Count=NULL,reg=NULL)
# for(name in names(toptags.genes.ensembl)){
#   tmp<-data.frame(Category=name,gene=g4.count.u100.dummy$gene,Count=g4.count.u100.dummy$Count,reg=g4.count.u100.dummy[,name])
#   g4.count.u100.reg<-rbind(g4.count.u100.reg,tmp)
# }
# 
# #write.csv(g4.count.u100.reg,file="./out/g4.count.u100.reg.csv",row.names = F,quote = F)
# 
# 
# 
