library(magrittr)
library(UpSetR)
library(readr)
library(dplyr)
library(grid)

setwd("./compare_svim_output/")

tbl <- read_tsv("key_combination_sv_cnt.tsv",skip=1,col_names = TRUE)
tbl <- as.data.frame(tbl)
colnames(tbl)[1] <- "sv_type"
tbl$sv_type %>% unique()

upset(tbl %>% filter(`sv_type` == "DEL"),sets = c("105","3A","6A"),keep.order = TRUE, text.scale=2)
grid.text("DEL",x = 0.10, y=0.95, gp=gpar(fontsize=15))

upset(tbl %>% filter(`sv_type` == "INS"),sets = c("105","3A","6A"),keep.order = TRUE, text.scale=2)
grid.text("INS",x = 0.10, y=0.95, gp=gpar(fontsize=15))

upset(tbl %>% filter(`sv_type` == "DUP_TAN"),sets = c("105","3A","6A"),keep.order = TRUE, text.scale=2)
grid.text("DUP_TAN",x = 0.10, y=0.95, gp=gpar(fontsize=15))

upset(tbl %>% filter(`sv_type` == "DUP_INT"),sets = c("105","3A","6A"),keep.order = TRUE, text.scale=2)
grid.text("DUP_INT",x = 0.10, y=0.95, gp=gpar(fontsize=15))
          