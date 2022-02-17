library(tidyverse)
library(GenomicRanges)
library(vroom)
library(rtracklayer)

options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#--------------------------------
hic_dat_in<-function(dat_file,cl_res,chromo){
  chr_dat<-vroom(paste0(dat_file,cl_res,"/",chromo,".txt"),delim = "\t",col_names = F,trim_ws = T,escape_double = F,show_col_types = F)
  return(chr_dat%>%mutate(X3=as.numeric(X3))%>%filter(!(is.nan(X3)))%>%filter(X1!=X2)%>%mutate(d=abs(X1-X2))%>%mutate(lw=log10(X3),ld=log10(d)))
}


breast_var_tbl<-vroom(file = "~/Documents/multires_bhicect/data/epi_data/VCF/icgc_brca_us.vcf",col_names = F)

var_Grange<-GRanges(seqnames=breast_var_tbl$X1,
                    ranges = IRanges(start=breast_var_tbl$X2,
                                     end=breast_var_tbl$X2)
)
mcols(var_Grange)<-breast_var_tbl %>% dplyr::select(X8)

chr_var<-breast_var_tbl %>% distinct(X1,X2) 

unique_var_Grange<-GRanges(seqnames=chr_var$X1,
                    ranges = IRanges(start=chr_var$X2,
                                     end=chr_var$X2)
)

#Convert hg38 to hg19
ch = import.chain("~/Documents/multires_bhicect/data/epi_data/hg38ToHg19.over.chain")

var_Grange_hg19<-unlist(liftOver(var_Grange, ch))
unique_var_Grange_hg19<-unlist(liftOver(unique_var_Grange, ch))
#--------------------------------
dat_file<-"~/Documents/multires_bhicect/data/HMEC/"
chromo<-"chr1"
tmp_res<-"5kb"
chr_set<-unlist(lapply(strsplit(list.files(paste(dat_file,tmp_res,sep='/')),split='\\.'),'[',1))
chr_bins<-do.call(bind_rows,lapply(chr_set,function(chromo){
  message(chromo)
  return(hic_dat_in(dat_file,tmp_res,chromo)%>% 
    summarise(bins=list(unique(c(X1,X2)))) %>% 
    unnest(cols=c(bins)) %>% arrange(bins) %>% 
    mutate(end=bins+res_num[tmp_res]-1,chrom=chromo)
  )
}))

chr_GRange<-GRanges(seqnames=chr_bins$chrom,
                    ranges = IRanges(start=chr_bins$bins,
                                     end=chr_bins$end)
)

#hist(log10(countOverlaps(chr_GRange,unique_var_Grange_hg19)))
spec_bin_count<-do.call(bind_rows,lapply(unique(mcols(var_Grange_hg19)$X8),function(ID){
  tibble(bin=1:length(chr_GRange),count=countOverlaps(chr_GRange,var_Grange_hg19[mcols(var_Grange_hg19)$X8==ID]),specimen=ID)
  
}))
bin_stat_tbl<-spec_bin_count %>% 
  group_by(bin) %>% 
  summarise(med=sum(count>0)/n(),mean=mean(count),var=var(count))

bin_stat_tbl %>% 
  #  filter(med>0.3)
  ggplot(.,aes(med))+geom_density()

bin_stat_tbl %>% 
  #  filter(med>0.3)
  ggplot(.,aes(mean))+geom_density()

bin_stat_tbl %>%
  ggplot(.,aes(x=med,y=var/mean))+geom_point(alpha=0.1)

bin_stat_tbl %>% 
  ggplot(.,aes(x=var/mean))+geom_histogram()+geom_vline(xintercept = 1)+scale_x_log10()
#-------------------------------------------------
bin_stat_tbl %>% 
  mutate(union.count=countOverlaps(chr_GRange,unique_var_Grange_hg19)) %>% 
#  filter(med>0) %>% 
  ggplot(.,aes(union.count,med))+
  geom_point(alpha=0.1)

bin_stat_tbl %>% 
  mutate(union.count=countOverlaps(chr_GRange,unique_var_Grange_hg19)) %>% 
#  filter(med>0) %>% 
  ggplot(.,aes(union.count,mean))+
  geom_point(alpha=0.1)


