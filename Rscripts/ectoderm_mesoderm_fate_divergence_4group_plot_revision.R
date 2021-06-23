library(tidyverse)
library(ggsci)
library(scales)

#load data
merged_m_vaf_tbl <- readRDS('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/merged_m_vaf_tbl.rds')
lng_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_lineage_210322.txt')
f_vaf_tbl <- readRDS('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/f_vaf_tbl.rds')



this_db = 'DB3'
fm_vaf_tbl <- f_vaf_tbl %>% filter(deadbody == this_db & time_of_lineage == 'early' & tissue_class != 'cancer_tissue')
fm_vaf_tbl <- fm_vaf_tbl %>% mutate(group2 = dominant_layer2)

ct_tbl <- fm_vaf_tbl %>% group_by(deadbody,group2, var_id) %>% count() %>% group_by(group2) %>% summarise(mean_n = mean(n))
print(ct_tbl)
incl_ids <- ct_tbl %>% filter(mean_n >=2) %>% pull(group2)
print(incl_ids)
fm_vaf_tbl <- fm_vaf_tbl %>% filter(group2 %in% incl_ids)

med_tbl <- fm_vaf_tbl %>% group_by(deadbody,group2, var_id) %>% summarise(medVAF = median(VAF), q0 = quantile(VAF)[2], q1=quantile(VAF)[4])
med_tbl <- left_join(med_tbl, fm_vaf_tbl %>% select(deadbody, var_id, lineage_id2, step_n) %>% unique())
med_tbl <- med_tbl %>% ungroup() %>%  mutate(first_group = apply(.["lineage_id2"], 1,function(x) unlist(strsplit(x, '-'))[1]))
var_order <- fm_vaf_tbl %>% group_by(deadbody,var_id) %>% summarise(medVAF = median(VAF)) %>% arrange(desc(medVAF)) %>% pull(var_id)

med_tbl %>% filter(group2%in%c("ectoderm", "endoderm", "mesoderm")) %>% 
ggplot(aes(x=var_id, y=log10(medVAF+0.0001), color = group2))+
  geom_point(size=3, alpha=0.1)+
  stat_smooth(aes(group = group2), alpha=0.3, size=0.5)+
  coord_cartesian(ylim=log10(c(0,0.5)+0.001))+
  #scale_color_manual(values=group_pal)+
  scale_color_manual(values=c("endoderm"=rgb(192/255,61/255,60/255), "mesoderm"=rgb(32/255,93/255,181/255),
                              "ectoderm"=rgb(234/255,181/255,10/255)))+
  scale_x_discrete(limits = var_order)+
  scale_y_continuous(expand = c(0,0.05), labels = c(1, 5,10, 0.1,0.01)*0.01*100, breaks= log10(c(1, 5,10,0.1,0.01)*0.01+0.001))+
  xlab("Early mutations") + 
  ylab("Log10(median VAF)") + 
  theme_minimal() + 
  theme(text=element_text(size=15), axis.text.x = element_blank(), panel.grid.major = element_blank(),
        #axis.line=element_line(color="black"), 
        panel.border = element_rect(color="black", fill=NA),
        axis.ticks.y=element_line(color="black"),
        legend.position="none"
  )

# 78 x 59 
# DB3:64
# DB6:522
# DB8:144
# DB9:150
# DB10:102
