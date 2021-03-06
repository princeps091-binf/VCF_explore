library(renv)

renv::init()

renv::install("tidyverse")
renv::install("YuLab-SMU/ggtree")

renv::install("bioc::ChIPseeker")
renv::install("bioc::GenomicRanges")
renv::install("bioc::org.Hs.eg.db")
renv::snapshot()
