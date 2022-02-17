library(tidyverse)
library(vroom)

breast_var_tbl<-vroom(file = "~/Documents/multires_bhicect/data/epi_data/VCF/icgc_brca_us.vcf",col_names = F)

breast_var_tbl %>% group_by(X8) %>% 
  summarise(n=n()) %>% 
  ggplot(.,aes(n))+geom_density()


specimen_rank_tbl<-breast_var_tbl %>% group_by(X8) %>% 
  summarise(n=n()) %>% 
  mutate(rank=min_rank(n)) %>% 
  arrange(rank)

gg_rain<-breast_var_tbl %>% 
  left_join(.,specimen_rank_tbl) %>% 
#  filter(X1=="chr17") %>% 
  mutate(X1=fct_relevel(X1,c(paste0("chr",1:22),"chrX","chrY"))) %>%   
  ggplot(.,aes(X2,y=rank))+geom_point(alpha=0.05,size=0.2)+
  xlab("chromosome position")+ylab("specimen")+
  theme_classic()+
  facet_wrap(X1 ~.,scales="free")
gg_rain
ggsave("~/Documents/multires_bhicect/weeklies/Fran_Supek/img/brca_us_mutation_rain_2.png",width = 60,height = 40,units = "cm",dpi = 700)


library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

chr_var<-breast_var_tbl %>% distinct(X1,X2) 

var_Grange<-GRanges(seqnames=chr_var$X1,
                    ranges = IRanges(start=chr_var$X2,
                                     end=chr_var$X2)
)

all_peakAnno <- annotatePeak(var_Grange, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db",verbose = F)

chr_var_fn_dist_tbl<-do.call(bind_rows,lapply(unique(chr_var$X1),function(chromo){
  message(chromo)
  sample_chr_var<-chr_var %>% 
    filter(X1 == chromo)
  
  sample_var_Grange<-GRanges(seqnames=sample_chr_var$X1,
                             ranges = IRanges(start=sample_chr_var$X2,
                                              end=sample_chr_var$X2)
  )
  
  
  peakAnno <- annotatePeak(sample_var_Grange, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db",verbose = F)
  
  return(peakAnno@annoStat %>% mutate(chr=chromo))
}))
chr_var_fn_dist_tbl %>% 
  bind_rows(.,all_peakAnno@annoStat %>% mutate(chr="ALL")) %>% 
  mutate(chr=fct_relevel(chr,c(paste0("chr",1:22),"chrX","chrY","ALL"))) %>%   
  ggplot(.,aes(chr,Frequency,fill=Feature))+
  geom_bar(stat = "identity")+scale_fill_brewer(palette="Paired")
ggsave("~/Documents/multires_bhicect/weeklies/Fran_Supek/img/brca_us_mutation_bar_chr.png",width = 30,height = 20,units = "cm")

chr17_var_tbl<-breast_var_tbl %>% filter(X1=="chr17")
chr17_var<-chr17_var_tbl %>% distinct(X1,X2)
chr17_var_Grange<-GRanges(seqnames=chr_var$X1,
                    ranges = IRanges(start=chr_var$X2,
                                     end=chr_var$X2)
)

chr17_peakAnno <- annotatePeak(chr17_var_Grange, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db",verbose = F)

chr17_var_fn_dist_tbl<-do.call(bind_rows,lapply(unique(chr17_var_tbl$X8),function(specimen){
  message(specimen)
  sample_chr_var<-chr17_var_tbl %>% 
    filter(X8 == specimen)
  
  sample_var_Grange<-GRanges(seqnames=sample_chr_var$X1,
                             ranges = IRanges(start=sample_chr_var$X2,
                                              end=sample_chr_var$X2)
  )
  
  
  peakAnno <- annotatePeak(sample_var_Grange, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db",verbose = F)
  
  return(peakAnno@annoStat %>% mutate(ID=specimen))
}))
chr17_var_fn_dist_tbl %>%   
  mutate(Feature=fct_relevel(Feature,as.character(tibble(chr17_peakAnno@annoStat)$Feature))) %>% 
  bind_rows(.,chr17_peakAnno@annoStat %>% mutate(ID="ALL")) %>% 
  ggplot(.,aes(ID,Frequency,fill=Feature))+
  geom_bar(stat = "identity")+scale_fill_brewer(palette="Paired")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  xlab("specimen")
ggsave("~/Documents/multires_bhicect/weeklies/Fran_Supek/img/brca_us_mutation_bar_chr17_specimen.png",width = 30,height = 20,units = "cm")
#------------------------------
#Mutation kinds
specimen_rank_tbl<-breast_var_tbl %>% group_by(X8) %>% 
  summarise(n=n()) %>% 
  mutate(rank=min_rank(n)) %>% 
  arrange(rank)

breast_var_tbl %>% 
  filter(!(grepl(",",X5))) %>% 
  mutate(mut=paste(X4,X5,sep="->")) %>% 
  group_by(mut) %>% 
  summarise(n=n()) %>% 
  mutate(mut=fct_reorder(mut,n,.desc=T)) %>% 
  ggplot(.,aes(mut,n))+geom_bar(stat="identity")
ggsave("~/Documents/multires_bhicect/weeklies/Fran_Supek/img/brca_us_mutation_kind_bar.png")

spec_mut_n_tbl<-breast_var_tbl %>% 
  filter(!(grepl(",",X5))) %>% 
  mutate(mut=paste(X4,X5,sep="->")) %>% 
  group_by(mut,X8) %>% 
  summarise(n=n()) %>% 
  mutate(mut=fct_reorder(mut,n)) 
spec_order<-spec_mut_n_tbl %>% ungroup %>% 
  group_by(X8) %>% mutate(p=n/sum(n)) %>% 
  filter(mut %in% c("C->T","G->A")) %>% 
  summarise(m=sum(p)) %>% 
  arrange(m) %>% 
  mutate(rank=1:n())
spec_mut_n_tbl %>% 
  left_join(.,spec_order) %>% 
  ggplot(.,aes(rank,n,fill=mut))+
  geom_bar(stat="identity",position="fill")+
  scale_fill_brewer(palette="Paired")+
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank())+
  xlab("specimen")
ggsave("~/Documents/multires_bhicect/weeklies/Fran_Supek/img/brca_us_mutation_kind_bar_by_sample.png")


breast_var_tbl %>% 
  filter(!(grepl(",",X5))) %>% 
  mutate(mut=paste(X4,X5,sep="->")) %>% 
  group_by(mut,X1) %>% 
  summarise(n=n()) %>% 
  mutate(mut=fct_reorder(mut,n)) %>% 
  ggplot(.,aes(X1,n,fill=mut))+geom_bar(stat="identity",position="fill")+scale_fill_brewer(palette="Paired")

breast_var_tbl %>% 
  filter(!(grepl(",",X5))) %>% 
  mutate(mut=paste(X4,X5,sep="->")) %>% 
  left_join(.,specimen_rank_tbl) %>% 
  filter(X1=="chr19") %>% 
  ggplot(.,aes(X2,y=rank))+geom_point(alpha=0.4,size=0.5)+
  xlab("chromosome position")+ylab("specimen")+
  theme_minimal()+facet_wrap(mut~.)

breast_var_kind_tbl<-breast_var_tbl %>% 
  filter(!(grepl(",",X5))) %>% 
  mutate(mut=paste(X4,X5,sep="->"))
  
var_kind_fn_dist_tbl<-do.call(bind_rows,lapply(unique(breast_var_kind_tbl$mut),function(specimen){
  message(specimen)
  sample_chr_var<-breast_var_kind_tbl %>% 
    filter(mut == specimen) %>% distinct(X1,X2)
  
  sample_var_Grange<-GRanges(seqnames=sample_chr_var$X1,
                             ranges = IRanges(start=sample_chr_var$X2,
                                              end=sample_chr_var$X2)
  )
  
  
  peakAnno <- annotatePeak(sample_var_Grange, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db",verbose = F)
  
  return(peakAnno@annoStat %>% mutate(ID=specimen))
}))

var_kind_fn_dist_tbl %>% 
  ggplot(.,aes(ID,Frequency,fill=Feature))+
  geom_bar(stat = "identity")+scale_fill_brewer(palette="Paired")+
  theme(axis.ticks.x=element_blank())+
  xlab("mutation kind")
ggsave("~/Documents/multires_bhicect/weeklies/Fran_Supek/img/brca_us_mutation_kind_bar_fn.png")
