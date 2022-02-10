library(tidyverse)
library(vroom)

breast_var_tbl<-vroom(file = "~/Documents/multires_bhicect/data/epi_data/VCF/ICGC_breast_snv_mnv_caller.tsv")

breast_var_tbl %>% distinct(icgc_mutation_id,icgc_donor_id) %>% 
  group_by(icgc_donor_id) %>% 
  summarise(n=n()) %>% 
  ggplot(.,aes(n))+
  scale_x_log10()+
  geom_density()

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

breast_cancer_mutation_tbl %>% 
  mutate(chromosome=fct_relevel(chromosome,c(1:22,"X"))) %>% 
  ggplot(.,aes(chromosome_start,y="chr"))+geom_point(alpha=0.01,size=0.1)+facet_grid(chromosome~.,scales="free")
#---------------------------------------------
