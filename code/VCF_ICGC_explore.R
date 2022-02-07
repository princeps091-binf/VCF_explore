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
