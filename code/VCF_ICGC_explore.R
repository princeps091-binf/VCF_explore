library(tidyverse)
library(vroom)

breast_var_tbl<-vroom(file = "~/Documents/multires_bhicect/data/epi_data/VCF/ICGC_breast_snv_mnv_caller.tsv")

breast_var_tbl %>% distinct(icgc_mutation_id,icgc_donor_id) %>% 
  group_by(icgc_donor_id) %>% 
  summarise(n=n()) %>% 
  ggplot(.,aes(n))+
  scale_x_log10()+
  geom_density()+xlab("number of mutations per specimen")
ggsave("~/Documents/multires_bhicect/weeklies/Fran_Supek/img/breast_variant_burden_density.png")

breast_var_tbl %>% distinct(icgc_mutation_id,icgc_donor_id) %>% 
  group_by(icgc_donor_id) %>% 
  summarise(n=n()) %>% 
  mutate(rank=percent_rank(n)) %>% 
  ggplot(.,aes(rank,n))+geom_line()+scale_y_log10()



breast_var_tbl %>% distinct(chromosome,chromosome_start,chromosome_end) %>% 
  group_by(chromosome) %>% 
  summarise(n=n()) %>% 
  mutate(chromosome=fct_relevel(chromosome,rev(c(1:22,"X")))) %>% 
  ggplot(.,aes(n,chromosome))+geom_point()


#---------------------------------------------
# Examine sample information
breast_sample_info<-vroom(file = "~/Documents/multires_bhicect/data/epi_data/VCF/icgc-dataset-clinical-info-breast/specimen.tsv.gz")
primary_tumor_id<-breast_sample_info %>%
  filter(specimen_type== "Primary tumour - solid tissue") %>% distinct(icgc_specimen_id,icgc_donor_id)
primary_tumor_id<-primary_tumor_id %>% 
  mutate(type="tumour")


breast_cancer_mutation_tbl<-breast_var_tbl %>% 
  filter(icgc_specimen_id %in% primary_tumor_id$icgc_specimen_id) %>% 
  distinct(icgc_mutation_id,chromosome,chromosome_start,chromosome_end,chromosome_strand)


breast_cancer_mutation_tbl<-breast_var_tbl %>% 
  filter(icgc_specimen_id %in% primary_tumor_id$icgc_specimen_id & chromosome == 1) 

specimen_rank_tbl<-breast_cancer_mutation_tbl %>% group_by(icgc_specimen_id) %>% 
  summarise(n=n()) %>% 
  mutate(rank=min_rank(n)) %>% 
  arrange(rank)
breast_cancer_mutation_tbl %>% 
  left_join(.,specimen_rank_tbl) %>% 
#  mutate(categ=ifelse(n<500,"low burden (<500 mutations)","high burden (>500 mutations)")) %>% 
  ggplot(.,aes(chromosome_start,y=rank))+geom_point(alpha=0.02,size=0.5)+
  xlab("chromosome 1 position")+ylab("specimen")+theme_minimal()#+facet_grid(categ ~.,scales="free")
ggsave("~/Documents/multires_bhicect/weeklies/Fran_Supek/img/chr1_breast_variant_waterfall.png")

#---------------------------------------------
# explanation of ICGC specimen metadata: https://docs.cancergenomicscloud.org/docs/icgc-metadata

library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

chr_var<-breast_cancer_mutation_tbl %>% distinct(chromosome,chromosome_start,chromosome_end) %>% mutate(chromosome=paste0("chr",chromosome))

var_Grange<-GRanges(seqnames=chr_var$chromosome,
                    ranges = IRanges(start=chr_var$chromosome_start,
                                     end=chr_var$chromosome_end)
)

sample_var_fn_dist_tbl<-do.call(bind_rows,lapply(unique(breast_cancer_mutation_tbl$icgc_specimen_id),function(ID){
  
  sample_chr_var<-breast_cancer_mutation_tbl %>% 
    filter(icgc_specimen_id == ID) %>% 
    distinct(chromosome,chromosome_start,chromosome_end) %>% mutate(chromosome=paste0("chr",chromosome))
  
  sample_var_Grange<-GRanges(seqnames=sample_chr_var$chromosome,
                             ranges = IRanges(start=sample_chr_var$chromosome_start,
                                              end=sample_chr_var$chromosome_end)
  )
  
  
  peakAnno <- annotatePeak(sample_var_Grange, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db",verbose = F)
  
  return(peakAnno@annoStat %>% mutate(sample=ID))
}))


specimen_rank_tbl<-breast_cancer_mutation_tbl %>% group_by(icgc_specimen_id) %>% 
  summarise(n=n()) %>% 
  mutate(rank=dense_rank(n)) %>% 
  arrange(rank) %>% mutate(pos=1:n())

sample_var_fn_dist_tbl %>%
  left_join(.,specimen_rank_tbl,by=c("sample"="icgc_specimen_id")) %>% 
  left_join(.,breast_sample_info %>% dplyr::select(icgc_specimen_id,specimen_donor_treatment_type ),by=c("sample"="icgc_specimen_id")) %>% 
  mutate(categ=ifelse(n<500,"low burden (<500 mutations)","high burden (>500 mutations)")) %>% 
  ggplot(.,aes(pos,Frequency,fill=Feature))+scale_fill_brewer(palette="Paired")+
  geom_bar(stat="identity")+
  facet_grid(. ~categ,scales="free")+ xlab("specimen")
ggsave("~/Documents/multires_bhicect/weeklies/Fran_Supek/img/chr1_breast_variant_fn_compo.png")


fn_repo<-"~/Documents/multires_bhicect/data/epi_data/fn_BED/"

fn_folder<-paste0(fn_repo,chromo,"/")
fn_file<-grep('BED$',list.files(fn_folder),value = T)
fn_bed_l<-lapply(fn_file,function(f){
  read_bed(paste0(fn_folder,f),n_fields = 3)
})
names(fn_bed_l)<-fn_file
txdb_chr <- txdb
seqlevels(txdb_chr)  <- chromo



rn_annotation<-sample(peakAnno@annoStat$Feature,size = length(chr_feature_Grange),prob = peakAnno@annoStat$Frequency/100,replace = T)
