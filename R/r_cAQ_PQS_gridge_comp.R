library(ggridges)


p.riges<-g4.count.u1kd1k.reg %>% 
  #filter(str_detect(Category,"1uM"))  %>%
  ggplot(aes(x=Count+1,y=Category,fill=reg,color=reg)) +
  geom_density_ridges(from=0,quantile_lines = TRUE, quantiles = 2,vline_size = 0.3,scale=0.95,bandwidth=0.1)+
  scale_x_continuous(name="PQS Counts",trans='log10') +
  scale_fill_manual(values = c("#0072B250","#D55E0050"), labels = c("DEG","non-DEG")) +
  scale_color_manual(values = c("#0072B2","#D55E00"), guide = "none") +
  scale_discrete_manual("point_color", values = c("#D55E00", "#0072B2"), guide = "none") +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(override.aes = list(fill = c("#D55E00A0", "#0072B2A0"),color = NA, point_color = NA)))+
  labs(title="Promoter(-1kbp ~ +1kbp)")+
  theme_ridges(center=TRUE)

p.riges
ggsave(file = "./out/geom_density_PQScounts2_u1kd1k_bw0.1.pdf", plot = p.riges, width = 11, height = 8.5)


p.riges<-g4.count.d1k.reg %>% 
  #filter(str_detect(Category,"1uM"))  %>%
  ggplot(aes(x=Count+1,y=Category,fill=reg,color=reg)) +
  geom_density_ridges(from=0,quantile_lines = TRUE, quantiles = 2,vline_size = 0.3,scale=0.95,bandwidth=0.1)+
  scale_x_continuous(name="PQS Counts",trans='log10') +
  scale_fill_manual(values = c("#0072B250","#D55E0050"), labels = c("DEG","non-DEG")) +
  scale_color_manual(values = c("#0072B2","#D55E00"), guide = "none") +
  scale_discrete_manual("point_color", values = c("#D55E00", "#0072B2"), guide = "none") +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(override.aes = list(fill = c("#D55E00A0", "#0072B2A0"),color = NA, point_color = NA)))+
  labs(title="Promoter(-1kbp ~ TSS)")+
  theme_ridges(center=TRUE)

p.riges
ggsave(file = "./out/geom_density_PQScounts2_d1k_bw0.1.pdf", plot = p.riges, width = 11, height = 8.5)

p.riges<-g4.count.u1k.reg %>% 
  #filter(str_detect(Category,"1uM"))  %>%
  ggplot(aes(x=Count+1,y=Category,fill=reg,color=reg)) +
  geom_density_ridges(from=0,quantile_lines = TRUE, quantiles = 2,vline_size = 0.3,scale=0.95,bandwidth=0.1)+
  scale_x_continuous(name="PQS Counts",trans='log10') +
  scale_fill_manual(values = c("#0072B250","#D55E0050"), labels = c("DEG","non-DEG")) +
  scale_color_manual(values = c("#0072B2","#D55E00"), guide = "none") +
  scale_discrete_manual("point_color", values = c("#D55E00", "#0072B2"), guide = "none") +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(override.aes = list(fill = c("#D55E00A0", "#0072B2A0"),color = NA, point_color = NA)))+
  labs(title="Promoter(TSS ~ +1kbp)")+
  theme_ridges(center=TRUE)

p.riges
ggsave(file = "./out/geom_density_PQScounts2_u1k_bw0.1.pdf", plot = p.riges, width = 11, height = 8.5)

p.riges<-g4.count.d100.reg %>% 
  #filter(str_detect(Category,"1uM"))  %>%
  ggplot(aes(x=Count+1,y=Category,fill=reg,color=reg)) +
  geom_density_ridges(from=0,quantile_lines = TRUE, quantiles = 2,vline_size = 0.3,scale=0.95,bandwidth=0.1)+
  scale_x_continuous(name="PQS Counts",trans='log10') +
  scale_fill_manual(values = c("#0072B250","#D55E0050"), labels = c("DEG","non-DEG")) +
  scale_color_manual(values = c("#0072B2","#D55E00"), guide = "none") +
  scale_discrete_manual("point_color", values = c("#D55E00", "#0072B2"), guide = "none") +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(override.aes = list(fill = c("#D55E00A0", "#0072B2A0"),color = NA, point_color = NA)))+
  labs(title="Promoter(TSS ~ +100bp)")+
  theme_ridges(center=TRUE)

p.riges
ggsave(file = "./out/geom_density_PQScounts2_d100_bw0.1.pdf", plot = p.riges, width = 11, height = 8.5)

p.riges<-g4.count.u100.reg %>% 
  #filter(str_detect(Category,"1uM"))  %>%
  ggplot(aes(x=Count+1,y=Category,fill=reg,color=reg)) +
  geom_density_ridges(from=0,quantile_lines = TRUE, quantiles = 2,vline_size = 0.3,scale=0.95,bandwidth=0.1)+
  scale_x_continuous(name="PQS Counts",trans='log10') +
  scale_fill_manual(values = c("#0072B250","#D55E0050"), labels = c("DEG","non-DEG")) +
  scale_color_manual(values = c("#0072B2","#D55E00"), guide = "none") +
  scale_discrete_manual("point_color", values = c("#D55E00", "#0072B2"), guide = "none") +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(override.aes = list(fill = c("#D55E00A0", "#0072B2A0"),color = NA, point_color = NA)))+
  labs(title="Promoter(-100bp ~ TSS)")+
  theme_ridges(center=TRUE)

p.riges
ggsave(file = "./out/geom_density_PQScounts2_u100_bw0.1.pdf", plot = p.riges, width = 11, height = 8.5)


#genebody-----------------------------------------------
p.riges<-g4.count.genebody.reg %>% 
  #filter(str_detect(Category,"1uM"))  %>%
  ggplot(aes(x=Count+1,y=Category,fill=reg,color=reg)) +
  geom_density_ridges(from=0,quantile_lines = TRUE, quantiles = 2,vline_size = 0.3,scale=0.95,bandwidth=0.1)+
  scale_x_continuous(name="PQS Counts",trans='log10') +
  scale_fill_manual(values = c("#0072B250","#D55E0050"), labels = c("DEG","non-DEG")) +
  scale_color_manual(values = c("#0072B2","#D55E00"), guide = "none") +
  scale_discrete_manual("point_color", values = c("#D55E00", "#0072B2"), guide = "none") +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(override.aes = list(fill = c("#D55E00A0", "#0072B2A0"),color = NA, point_color = NA)))+
  labs(title="Gene body")+
  theme_ridges(center=TRUE)

p.riges
ggsave(file = "./out/geom_density_PQScounts2_genebody_bw0.1.pdf", plot = p.riges, width = 11, height = 8.5)

p.riges<-g4.count.intron.reg %>% 
  #filter(str_detect(Category,"1uM"))  %>%
  ggplot(aes(x=Count+1,y=Category,fill=reg,color=reg)) +
  geom_density_ridges(from=0,quantile_lines = TRUE, quantiles = 2,vline_size = 0.3,scale=0.95,bandwidth=0.1)+
  scale_x_continuous(name="PQS Counts",trans='log10') +
  scale_fill_manual(values = c("#0072B250","#D55E0050"), labels = c("DEG","non-DEG")) +
  scale_color_manual(values = c("#0072B2","#D55E00"), guide = "none") +
  scale_discrete_manual("point_color", values = c("#D55E00", "#0072B2"), guide = "none") +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(override.aes = list(fill = c("#D55E00A0", "#0072B2A0"),color = NA, point_color = NA)))+
  labs(title="Intron")+
  theme_ridges(center=TRUE)

p.riges
ggsave(file = "./out/geom_density_PQScounts2_intron_bw0.1.pdf", plot = p.riges, width = 11, height = 8.5)

p.riges<-g4.count.1stintron.reg %>% 
  #filter(str_detect(Category,"1uM"))  %>%
  ggplot(aes(x=Count+1,y=Category,fill=reg,color=reg)) +
  geom_density_ridges(from=0,quantile_lines = TRUE, quantiles = 2,vline_size = 0.3,scale=0.95,bandwidth=0.1)+
  scale_x_continuous(name="PQS Counts",trans='log10') +
  scale_fill_manual(values = c("#0072B250","#D55E0050"), labels = c("DEG","non-DEG")) +
  scale_color_manual(values = c("#0072B2","#D55E00"), guide = "none") +
  scale_discrete_manual("point_color", values = c("#D55E00", "#0072B2"), guide = "none") +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(override.aes = list(fill = c("#D55E00A0", "#0072B2A0"),color = NA, point_color = NA)))+
  labs(title="1st Intron")+
  theme_ridges(center=TRUE)

p.riges
ggsave(file = "./out/geom_density_PQScounts2_1stintron_bw0.1.pdf", plot = p.riges, width = 11, height = 8.5)

p.riges<-g4.count.exon.reg %>% 
  #filter(str_detect(Category,"1uM"))  %>%
  ggplot(aes(x=Count+1,y=Category,fill=reg,color=reg)) +
  geom_density_ridges(from=0,quantile_lines = TRUE, quantiles = 2,vline_size = 0.3,scale=0.95,bandwidth=0.1)+
  scale_x_continuous(name="PQS Counts",trans='log10') +
  scale_fill_manual(values = c("#0072B250","#D55E0050"), labels = c("DEG","non-DEG")) +
  scale_color_manual(values = c("#0072B2","#D55E00"), guide = "none") +
  scale_discrete_manual("point_color", values = c("#D55E00", "#0072B2"), guide = "none") +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(override.aes = list(fill = c("#D55E00A0", "#0072B2A0"),color = NA, point_color = NA)))+
  labs(title="Exon")+
  theme_ridges(center=TRUE)

p.riges
ggsave(file = "./out/geom_density_PQScounts2_exon_bw0.1.pdf", plot = p.riges, width = 11, height = 8.5)

p.riges<-g4.count.1stexon.reg %>% 
  #filter(str_detect(Category,"1uM"))  %>%
  ggplot(aes(x=Count+1,y=Category,fill=reg,color=reg)) +
  geom_density_ridges(from=0,quantile_lines = TRUE, quantiles = 2,vline_size = 0.3,scale=0.95,bandwidth=0.1)+
  scale_x_continuous(name="PQS Counts",trans='log10') +
  scale_fill_manual(values = c("#0072B250","#D55E0050"), labels = c("DEG","non-DEG")) +
  scale_color_manual(values = c("#0072B2","#D55E00"), guide = "none") +
  scale_discrete_manual("point_color", values = c("#D55E00", "#0072B2"), guide = "none") +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(override.aes = list(fill = c("#D55E00A0", "#0072B2A0"),color = NA, point_color = NA)))+
  labs(title="1st Exon")+
  theme_ridges(center=TRUE)

p.riges
ggsave(file = "./out/geom_density_PQScounts2_1stexon_bw0.1.pdf", plot = p.riges, width = 11, height = 8.5)

p.riges<-g4.count.5utr.reg %>% 
  #filter(str_detect(Category,"1uM"))  %>%
  ggplot(aes(x=Count+1,y=Category,fill=reg,color=reg)) +
  geom_density_ridges(from=0,quantile_lines = TRUE, quantiles = 2,vline_size = 0.3,scale=0.95,bandwidth=0.1)+
  scale_x_continuous(name="PQS Counts",trans='log10') +
  scale_fill_manual(values = c("#0072B250","#D55E0050"), labels = c("DEG","non-DEG")) +
  scale_color_manual(values = c("#0072B2","#D55E00"), guide = "none") +
  scale_discrete_manual("point_color", values = c("#D55E00", "#0072B2"), guide = "none") +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(override.aes = list(fill = c("#D55E00A0", "#0072B2A0"),color = NA, point_color = NA)))+
  labs(title="5' UTR")+
  theme_ridges(center=TRUE)

p.riges
ggsave(file = "./out/geom_density_PQScounts2_5utr_bw0.1.pdf", plot = p.riges, width = 11, height = 8.5)

p.riges<-g4.count.3utr.reg %>% 
  #filter(str_detect(Category,"1uM"))  %>%
  ggplot(aes(x=Count+1,y=Category,fill=reg,color=reg)) +
  geom_density_ridges(from=0,quantile_lines = TRUE, quantiles = 2,vline_size = 0.3,scale=0.95,bandwidth=0.1)+
  scale_x_continuous(name="PQS Counts",trans='log10') +
  scale_fill_manual(values = c("#0072B250","#D55E0050"), labels = c("DEG","non-DEG")) +
  scale_color_manual(values = c("#0072B2","#D55E00"), guide = "none") +
  scale_discrete_manual("point_color", values = c("#D55E00", "#0072B2"), guide = "none") +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(override.aes = list(fill = c("#D55E00A0", "#0072B2A0"),color = NA, point_color = NA)))+
  labs(title="3' UTR")+
  theme_ridges(center=TRUE)


p.riges
ggsave(file = "./out/geom_density_PQScounts2_3utr_bw0.1.pdf", plot = p.riges, width = 11, height = 8.5)


#div_genelen----------------------------------------
p.riges<-g4.count.u1kd1k.reg %>% 
  #filter(str_detect(Category,"1uM"))  %>%
  ggplot(aes(x=Count/2000,y=Category,fill=reg,color=reg)) +
  geom_density_ridges(from=0,to=0.01,
                      quantile_lines = TRUE, quantiles = 2,vline_size = 0.3,scale=0.95,bandwidth=0.0005)+
  scale_x_continuous(name="PQS Counts/Promoter") +
  scale_fill_manual(values = c("#0072B250","#D55E0050"), labels = c("DEG","non-DEG")) +
  scale_color_manual(values = c("#0072B2","#D55E00"), guide = "none") +
  scale_discrete_manual("point_color", values = c("#D55E00", "#0072B2"), guide = "none") +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(override.aes = list(fill = c("#D55E00A0", "#0072B2A0"),color = NA, point_color = NA)))+
  labs(title="Promoter(-1kbp ~ +1kbp)")+
  theme_ridges(center=TRUE)

p.riges
ggsave(file = "./out/geom_density_PQScounts2_u1kd1k_perbp_bw5e-4.pdf", plot = p.riges, width = 11, height = 8.5)

p.riges<-g4.count.d1k.reg %>% 
  #filter(str_detect(Category,"1uM"))  %>%
  ggplot(aes(x=Count/1000,y=Category,fill=reg,color=reg)) +
  geom_density_ridges(from=0,to=0.01,
                      quantile_lines = TRUE, quantiles = 2,vline_size = 0.3,scale=0.95,bandwidth=0.0005)+
  scale_x_continuous(name="PQS Counts/Promoter") +
  scale_fill_manual(values = c("#0072B250","#D55E0050"), labels = c("DEG","non-DEG")) +
  scale_color_manual(values = c("#0072B2","#D55E00"), guide = "none") +
  scale_discrete_manual("point_color", values = c("#D55E00", "#0072B2"), guide = "none") +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(override.aes = list(fill = c("#D55E00A0", "#0072B2A0"),color = NA, point_color = NA)))+
  labs(title="Promoter(-1kbp ~ TSS)")+
  theme_ridges(center=TRUE)

p.riges
ggsave(file = "./out/geom_density_PQScounts2_d1k_perbp_bw5e-4.pdf", plot = p.riges, width = 11, height = 8.5)

p.riges<-g4.count.u1k.reg %>% 
  #filter(str_detect(Category,"1uM"))  %>%
  ggplot(aes(x=Count/1000,y=Category,fill=reg,color=reg)) +
  geom_density_ridges(from=0,to=0.01,
                      quantile_lines = TRUE, quantiles = 2,vline_size = 0.3,scale=0.95,bandwidth=0.0005)+
  scale_x_continuous(name="PQS Counts/Promoter") +
  scale_fill_manual(values = c("#0072B250","#D55E0050"), labels = c("DEG","non-DEG")) +
  scale_color_manual(values = c("#0072B2","#D55E00"), guide = "none") +
  scale_discrete_manual("point_color", values = c("#D55E00", "#0072B2"), guide = "none") +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(override.aes = list(fill = c("#D55E00A0", "#0072B2A0"),color = NA, point_color = NA)))+
  labs(title="Promoter(TSS ~ +1kbp)")+
  theme_ridges(center=TRUE)

p.riges
ggsave(file = "./out/geom_density_PQScounts2_u1k_perbp_bw5e-4.pdf", plot = p.riges, width = 11, height = 8.5)


p.riges<-g4.count.genebody.reg %>% 
  #filter(str_detect(Category,"1uM"))  %>%
  ggplot(aes(x=Count/Width,y=Category,fill=reg,color=reg)) +
  geom_density_ridges(from=0,to=0.002,
                      quantile_lines = TRUE, quantiles = 2,vline_size = 0.3,scale=0.95,bandwidth=0.0001)+
  scale_x_continuous(name="PQS Counts/Gene length") +
  scale_fill_manual(values = c("#0072B250","#D55E0050"), labels = c("DEG","non-DEG")) +
  scale_color_manual(values = c("#0072B2","#D55E00"), guide = "none") +
  scale_discrete_manual("point_color", values = c("#D55E00", "#0072B2"), guide = "none") +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(override.aes = list(fill = c("#D55E00A0", "#0072B2A0"),color = NA, point_color = NA)))+
  labs(title="Gene body")+
  theme_ridges(center=TRUE)

p.riges
ggsave(file = "./out/geom_density_PQScounts2_genebody_perbp_bw1e-4.pdf", plot = p.riges, width = 11, height = 8.5)

p.riges<-g4.count.intron.reg %>% 
  #filter(str_detect(Category,"1uM"))  %>%
  ggplot(aes(x=Count/Width,y=Category,fill=reg,color=reg)) +
  geom_density_ridges(from=0,to=0.002,
                      quantile_lines = TRUE, quantiles = 2,vline_size = 0.3,scale=0.95,bandwidth=0.0001)+
  scale_x_continuous(name="PQS Counts/Gene length") +
  scale_fill_manual(values = c("#0072B250","#D55E0050"), labels = c("DEG","non-DEG")) +
  scale_color_manual(values = c("#0072B2","#D55E00"), guide = "none") +
  scale_discrete_manual("point_color", values = c("#D55E00", "#0072B2"), guide = "none") +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(override.aes = list(fill = c("#D55E00A0", "#0072B2A0"),color = NA, point_color = NA)))+
  labs(title="Intron")+
  theme_ridges(center=TRUE)

p.riges
ggsave(file = "./out/geom_density_PQScounts2_intron_perbp_perbp_bw1e-4.pdf", plot = p.riges, width = 11, height = 8.5)

p.riges<-g4.count.1stintron.reg %>% 
  #filter(str_detect(Category,"1uM"))  %>%
  ggplot(aes(x=Count/Width,y=Category,fill=reg,color=reg)) +
  geom_density_ridges(from=0,to=0.002,
                      quantile_lines = TRUE, quantiles = 2,vline_size = 0.3,scale=0.95,bandwidth=0.0001)+
  scale_x_continuous(name="PQS Counts/Gene length") +
  scale_fill_manual(values = c("#0072B250","#D55E0050"), labels = c("DEG","non-DEG")) +
  scale_color_manual(values = c("#0072B2","#D55E00"), guide = "none") +
  scale_discrete_manual("point_color", values = c("#D55E00", "#0072B2"), guide = "none") +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(override.aes = list(fill = c("#D55E00A0", "#0072B2A0"),color = NA, point_color = NA)))+
  labs(title="1st Intron")+
  theme_ridges(center=TRUE)

p.riges
ggsave(file = "./out/geom_density_PQScounts2_1stintron_perbp_bw1e-4.pdf", plot = p.riges, width = 11, height = 8.5)

p.riges<-g4.count.exon.reg %>% 
  #filter(str_detect(Category,"1uM"))  %>%
  ggplot(aes(x=Count/Width,y=Category,fill=reg,color=reg)) +
  geom_density_ridges(from=0,to=0.002,
                      quantile_lines = TRUE, quantiles = 2,vline_size = 0.3,scale=0.95,bandwidth=0.0001)+
  scale_x_continuous(name="PQS Counts/Gene length") +
  scale_fill_manual(values = c("#0072B250","#D55E0050"), labels = c("DEG","non-DEG")) +
  scale_color_manual(values = c("#0072B2","#D55E00"), guide = "none") +
  scale_discrete_manual("point_color", values = c("#D55E00", "#0072B2"), guide = "none") +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(override.aes = list(fill = c("#D55E00A0", "#0072B2A0"),color = NA, point_color = NA)))+
  labs(title="Exon")+
  theme_ridges(center=TRUE)

p.riges
ggsave(file = "./out/geom_density_PQScounts2_exon_perbp_bw1e-4.pdf", plot = p.riges, width = 11, height = 8.5)

p.riges<-g4.count.1stexon.reg %>% 
  #filter(str_detect(Category,"1uM"))  %>%
  ggplot(aes(x=Count/Width,y=Category,fill=reg,color=reg)) +
  geom_density_ridges(from=0,to=0.002,
                      quantile_lines = TRUE, quantiles = 2,vline_size = 0.3,scale=0.95,bandwidth=0.0001)+
  scale_x_continuous(name="PQS Counts/Gene length") +
  scale_fill_manual(values = c("#0072B250","#D55E0050"), labels = c("DEG","non-DEG")) +
  scale_color_manual(values = c("#0072B2","#D55E00"), guide = "none") +
  scale_discrete_manual("point_color", values = c("#D55E00", "#0072B2"), guide = "none") +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(override.aes = list(fill = c("#D55E00A0", "#0072B2A0"),color = NA, point_color = NA)))+
  labs(title="1st Exon")+
  theme_ridges(center=TRUE)

p.riges
ggsave(file = "./out/geom_density_PQScounts2_1stexon_perbp_bw1e-4.pdf", plot = p.riges, width = 11, height = 8.5)

p.riges<-g4.count.5utr.reg %>% 
  #filter(str_detect(Category,"1uM"))  %>%
  ggplot(aes(x=Count/(Width+10),y=Category,fill=reg,color=reg)) +
  geom_density_ridges(from=0,to=0.002,
                      quantile_lines = TRUE, quantiles = 2,vline_size = 0.3,scale=0.95,bandwidth=0.0001)+
  scale_x_continuous(name="PQS Counts/Gene length") +
  scale_fill_manual(values = c("#0072B250","#D55E0050"), labels = c("DEG","non-DEG")) +
  scale_color_manual(values = c("#0072B2","#D55E00"), guide = "none") +
  scale_discrete_manual("point_color", values = c("#D55E00", "#0072B2"), guide = "none") +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(override.aes = list(fill = c("#D55E00A0", "#0072B2A0"),color = NA, point_color = NA)))+
  labs(title="5' UTR")+
  theme_ridges(center=TRUE)

p.riges
ggsave(file = "./out/geom_density_PQScounts2_5utr_perbp_bw1e-4.pdf", plot = p.riges, width = 11, height = 8.5)

p.riges<-g4.count.3utr.reg %>% 
  #filter(str_detect(Category,"1uM"))  %>%
  ggplot(aes(x=Count/Width,y=Category,fill=reg,color=reg)) +
  geom_density_ridges(from=0,to=0.002,
                      quantile_lines = TRUE, quantiles = 2,vline_size = 0.3,scale=0.95,bandwidth=0.0001)+
  scale_x_continuous(name="PQS Counts/Gene length") +
  scale_fill_manual(values = c("#0072B250","#D55E0050"), labels = c("DEG","non-DEG")) +
  scale_color_manual(values = c("#0072B2","#D55E00"), guide = "none") +
  scale_discrete_manual("point_color", values = c("#D55E00", "#0072B2"), guide = "none") +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(override.aes = list(fill = c("#D55E00A0", "#0072B2A0"),color = NA, point_color = NA)))+
  labs(title="3' UTR")+
  theme_ridges(center=TRUE)


p.riges
ggsave(file = "./out/geom_density_PQScounts2_3utr_perbp_bw1e-4.pdf", plot = p.riges, width = 11, height = 8.5)
