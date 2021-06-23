###########################<Common processing>##############################################
#load library
library(tidyverse)
library(ggsci)
library(ggtree)
library(ggrepel)
library(cowplot)

#load dataset
meta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_sample_210322.txt') %>% filter(current_final_for_lineage == 'Y')
raw_meta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_sample_210322.txt')
lng_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_lineage_210322.txt')
DB_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_DB_210322.txt')
rate_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/Mutation_rate_simulation/df_forest_plot.txt')
pmeta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_TargetSeq_210325.txt')



# path assign
nwk_tbl <- tribble(
  ~deadbody, ~nwk_path,
  "DB2", '/home/users/sypark/00_Project/06_LineageTracing/db2/10_point_mutations/DB2_Lineage_count_table.txt.nwk', 
  "DB3", '/home/users/sypark/00_Project/06_LineageTracing/db3/14_point_mutations/DB3_Lineage_count_table.txt.nwk',
  "DB5",'/home/users/sypark/00_Project/06_LineageTracing/db5/11_point_mutations/DB5_Lineage_count_table.txt.nwk',
  "DB6",'/home/users/sypark/00_Project/06_LineageTracing/db6/12_point_mutations/DB6_Lineage_count_table.txt.nwk',
  "DB8",'/home/users/sypark/00_Project/06_LineageTracing/db8/09_point_mutations/DB8_Lineage_count_table.txt.nwk',
  "DB9",'/home/users/sypark/00_Project/06_LineageTracing/db9/09_point_mutations/DB9_Lineage_count_table.txt.nwk',
  "DB10",'/home/users/sypark/00_Project/06_LineageTracing/db10/09_point_mutations/DB10_Lineage_count_table.txt.nwk'
)
pointmt_path_tbl <- tribble(
  ~deadbody, ~pointmt_path,
  "DB2","/home/users/sypark/00_Project/06_LineageTracing/db2/10_point_mutations/",
  "DB3","/home/users/sypark/00_Project/06_LineageTracing/db3/14_point_mutations/",
  "DB5","/home/users/sypark/00_Project/06_LineageTracing/db5/11_point_mutations/",
  "DB6","/home/users/sypark/00_Project/06_LineageTracing/db6/12_point_mutations/",
  "DB8","/home/users/sypark/00_Project/06_LineageTracing/db8/09_point_mutations/",
  "DB9","/home/users/sypark/00_Project/06_LineageTracing/db9/09_point_mutations/",
  "DB10","/home/users/sypark/00_Project/06_LineageTracing/db10/09_point_mutations/"
)

#edit meta_dt
meta_dt <- meta_dt %>% mutate(Source3 = ifelse(is.na(Source_side)==T, Source2, paste0(Source_side,'_',Source2)))
meta_dt$deadbody <- factor(meta_dt$deadbody, levels =c('DB2','DB3','DB5','DB6','DB8','DB9','DB10'))

#edit lng_dt
lng_dt$deadbody <- factor(lng_dt$deadbody, levels =c('DB2','DB3','DB5','DB6','DB8','DB9','DB10'))

#edit DB_dt
DB_dt$deadbody <- factor(DB_dt$deadbody, levels =c('DB2','DB3','DB5','DB6','DB8','DB9','DB10'))



#Fx source
source('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/merged_analysis_functions.R')

# --------

library(colorRamps)

this_db="DB9"

total_break=100
col_breaks = seq(0, 100, length.out=total_break)
a <- length(col_breaks[col_breaks < 1] )  #change color at 1
b <- length(col_breaks[col_breaks < 10])-a  #change color at 10
#rampCol1 <- colorRampPalette(c("gray", "yellow"))(n=a)
rampCol1 <- colorRampPalette(c("gray", rgb(250/255, 218/255, 107/255)))(n=a)
#rampCol2 <- colorRampPalette(c("yellow","red"))(n=b)
rampCol2 <- colorRampPalette(c(rgb(250/255, 218/255, 107/255), rgb(198/255,40/255,89/255)))(n=b)
#rampCol3 <- colorRampPalette(c("red", "blue"))(n=total_break-(a+b))
rampCol3 <- colorRampPalette(c(rgb(198/255,40/255,89/255), rgb(45/255,65/255,167/255)))(n=total_break-(a+b))
my_palette<- c(rampCol1, rampCol2, rampCol3)

file_list <- list.files(paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == this_db], 'perLineage/'), pattern='^L.*.pointmts.txt$', full.names = T)
dt <- do.call(rbind,lapply(file_list, function(x) read_tsv(x, col_types = cols(`#CHROM`='c')) %>% mutate(lineage_id2 = paste0(this_db, '_', gsub('.pointmts.txt','',unlist(strsplit(x,'//'))[2])))))
dt <- dt %>% mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep="_")) %>% separate(`blood_dp;blood_var;blood_vafpct`, c('blood_dp','blood_var','blood_vafpct'), sep=';',convert=T)
dt <- left_join(dt, lng_dt %>% dplyr::select(lineage_id2, start_mtn, time_of_lineage))
dt <- dt %>% mutate(adjblVAF = ifelse(time_of_lineage =='early' | blood_var > 1, blood_vafpct,0))   # ignore when the lineage is not early and blood_var == 1
for (t_lid in unique(dt$lineage_id2)){   #ignore when the lowest 3 blVAFs of the upper stream are 0
  if (dt %>% filter(lineage_id2 == t_lid) %>% arrange(adjblVAF) %>% head(3) %>% pull(adjblVAF) %>% max() == 0){
    dt$adjblVAF[grepl(paste0(t_lid, '-'),dt$lineage_id2)] <- 0
  }
}
max_VAF_len=20
res_tbl <- as.tibble(data.frame(matrix(vector(), 0, max_VAF_len + 1,
                                       dimnames=list(c(), c('lineage_id2',paste0('VAF',1:max_VAF_len)))),
                                stringsAsFactors=F))

n=0
for (lid in unique(dt$lineage_id2)){
  n=n+1
  VAF_v <- dt %>% filter(lineage_id2 == lid) %>% arrange(desc(adjblVAF)) %>% pull(adjblVAF)%>% .[1:max_VAF_len]
  #names(VAF_v) <- paste0('VAF', 1:length(VAF_v))
  NA_v <- rep(NA, max_VAF_len - length(VAF_v))
  #names(NA_v) <- pate0('VAF', (length(VAF_v)+1):max_VAF_len)
  temp <- c(lid,VAF_v, NA_v)
  names(temp) <- c('lineage_id2',paste0('VAF',1:max_VAF_len))
  #print(enframe(temp) %>% spread(name, value) %>% dplyr::select(lineage_id2, paste0('VAF', 1:max_VAF_len)))
  res_tbl <- bind_rows(res_tbl, enframe(temp) %>% spread(name, value) %>% dplyr::select(lineage_id2, paste0('VAF', 1:max_VAF_len)))
}


res_tbl <- res_tbl %>% mutate_at(vars(starts_with('VAF')), list(as.numeric))
tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
m_lng_dt <- lng_dt %>% filter(deadbody == this_db) %>% rename(taxa = lineage_id) %>% dplyr::select(-deadbody)
m_lng_dt <- left_join(m_lng_dt, res_tbl)

ggtree(tree, size=0.3) %<+% m_lng_dt + theme_tree2()+
  geom_nodelab() +
  coord_cartesian(xlim = c(-1,30))+
  geom_segment2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt, xend = x -n_pointmt +1, y = y, yend = y, color = VAF1), size=1.5)+
  geom_segment2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt+1, xend = x -n_pointmt +2, y = y, yend = y, color = VAF2), size=1.5)+
  geom_segment2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt+2, xend = x -n_pointmt +3, y = y, yend = y, color = VAF3), size=1.5)+
  geom_segment2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt+3, xend = x -n_pointmt +4, y = y, yend = y, color = VAF4), size=1.5)+
  geom_segment2(aes(subset = is.na(VAF5) == F, x=x-n_pointmt+4, xend = x -n_pointmt +5, y = y, yend = y, color = VAF5), size=1.5)+
  geom_segment2(aes(subset = is.na(VAF6) == F, x=x-n_pointmt+5, xend = x -n_pointmt +6, y = y, yend = y, color = VAF6), size=1.5)+
  geom_segment2(aes(subset = is.na(VAF7) == F, x=x-n_pointmt+6, xend = x -n_pointmt +7, y = y, yend = y, color = VAF7), size=1.5)+
  geom_segment2(aes(subset = is.na(VAF8) == F, x=x-n_pointmt+7, xend = x -n_pointmt +8, y = y, yend = y, color = VAF8), size=1.5)+
  geom_segment2(aes(subset = is.na(VAF9) == F, x=x-n_pointmt+8, xend = x -n_pointmt +9, y = y, yend = y, color = VAF9), size=1.5)+
  geom_segment2(aes(subset = is.na(VAF10) == F, x=x-n_pointmt+9, xend = x -n_pointmt +10, y = y, yend = y, color = VAF10), size=1.5)+
  geom_segment2(aes(subset = is.na(VAF11) == F, x=x-n_pointmt+10, xend = x -n_pointmt +11, y = y, yend = y, color = VAF11), size=1.5)+
  geom_segment2(aes(subset = is.na(VAF12) == F, x=x-n_pointmt+11, xend = x -n_pointmt +12, y = y, yend = y, color = VAF12), size=1.5)+
  geom_segment2(aes(subset = is.na(VAF13) == F, x=x-n_pointmt+12, xend = x -n_pointmt +13, y = y, yend = y, color = VAF13), size=1.5)+
  geom_segment2(aes(subset = is.na(VAF14) == F, x=x-n_pointmt+13, xend = x -n_pointmt +14, y = y, yend = y, color = VAF14), size=1.5)+
  geom_segment2(aes(subset = is.na(VAF15) == F, x=x-n_pointmt+14, xend = x -n_pointmt +15, y = y, yend = y, color = VAF15), size=1.5)+
  geom_segment2(aes(subset = is.na(VAF16) == F, x=x-n_pointmt+15, xend = x -n_pointmt +16, y = y, yend = y, color = VAF16), size=1.5)+
  geom_segment2(aes(subset = is.na(VAF17) == F, x=x-n_pointmt+16, xend = x -n_pointmt +17, y = y, yend = y, color = VAF17), size=1.5)+
  geom_segment2(aes(subset = is.na(VAF18) == F, x=x-n_pointmt+17, xend = x -n_pointmt +18, y = y, yend = y, color = VAF18), size=1.5)+
  geom_segment2(aes(subset = is.na(VAF19) == F, x=x-n_pointmt+18, xend = x -n_pointmt +19, y = y, yend = y, color = VAF19), size=1.5)+
  geom_segment2(aes(subset = is.na(VAF20) == F, x=x-n_pointmt+19, xend = x -n_pointmt +20, y = y, yend = y, color = VAF20), size=1.5)+
  scale_color_gradientn(colors = my_palette, name="VAF", breaks=c(0,25,50))+
  theme(
    axis.text = element_text(size=10, angle=90), 
    plot.title = element_text(size=20), 
    plot.margin= unit(c(0.0,0,0.0,0.0),"cm"),
    axis.line.x=element_blank(), 
    legend.text = element_text(size=10, angle=90),
    legend.title = element_blank(),
    legend.position=c(0.10, 0.8))




m_lng_dt %>% filter(end_mtn>1) %>% 
  dplyr::select(lineage_id2, starts_with("VAF")) %>% 
  mutate(max=apply(.[,2:ncol(.)], 1, max)*2) %>% 
  dplyr::select(lineage_id2, max) %>% arrange(desc(max)) %>% 
  mutate(lineage_id2=factor(lineage_id2, levels=lineage_id2)) %>% 
  ggplot(aes(x=lineage_id2, y=max)) + 
  #geom_bar(stat="identity", width=1.0, fill=rgb(203/255,104/255,41/255), color=rgb(203/255,104/255,41/255))+ #DB2
  #geom_bar(stat="identity", width=1.0, fill=rgb(36/255,122/255,191/255), color=rgb(36/255,122/255,191/255)) + #DB3
  #geom_bar(stat="identity", width=1.0, fill=rgb(141/255,34/255,40/255), color=rgb(141/255,34/255,40/255)) + #DB5
  #geom_bar(stat="identity", width=1.0, fill=rgb(87/255,145/255,114/255), color=rgb(87/255,145/255,114/255)) + #DB6
  #geom_bar(stat="identity", width=1.0, fill=rgb(69/255,60/255,116/255), color=rgb(69/255,60/255,116/255)) + #DB8
  geom_bar(stat="identity", width=1.0, fill=rgb(89/255,83/255,70/255), color=rgb(89/255,83/255,70/255)) + #DB9
  #geom_bar(stat="identity", width=1.0, fill=rgb(32/255,47/255,52/255), color=rgb(32/255,47/255,52/255)) + #DB10
  geom_hline(yintercept=100/nrow(m_lng_dt %>% filter(end_mtn>30)), color="red", linetype="dashed") + 
  theme_minimal() +
  #ylab("Maximum VAF") + 
  scale_y_continuous(expand=c(0,0))+
  theme(text=element_text(size=20), 
        axis.text.x=element_blank(), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.border = element_rect(color="black", fill=NA)) 
