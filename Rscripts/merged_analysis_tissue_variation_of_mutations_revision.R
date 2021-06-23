###########################<Common processing>##############################################
#load library
library(tidyverse)
library(ggsci)
library(ggtree)
library(ggrepel)
library(cowplot)
library(binom)

#file paths
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

targetseq_path_tbl <- tribble(
  ~deadbody, ~folder_path, ~Call_name, 
  "DB3",'/home/users/sypark/00_Project/06_LineageTracing/db3/04-5_Target_seq/05_re_processing/','DB3_TargetSeq_Target_list.GSNPadded.60sCall',
  "DB6",'/home/users/sypark/00_Project/06_LineageTracing/db6/14_TargetSeq/','DB6_TargetSeq_Target_list.merged.GSNPadded.143sCall',
  "DB8",'/home/users/sypark/00_Project/06_LineageTracing/db8/11_TargetSeq/','DB8_TargetSeq_Target_list.merged.GSNPadded.44sCall',
  "DB9",'/home/users/sypark/00_Project/06_LineageTracing/db9/11_TargetSeq/','DB9_TargetSeq_Target_list.merged.GSNPadded.43sCall',
  "DB10",'/home/users/sypark/00_Project/06_LineageTracing/db10/11_TargetSeq/','DB10_TargetSeq_Target_list.merged.GSNPadded.89sCall'
)




#load data
pmeta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_TargetSeq_210228.txt')
lng_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_lineage_210322.txt')
rate_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/Mutation_rate_simulation/df_forest_plot.txt')
meta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_sample_210322.txt') %>% filter(current_final_for_lineage == 'Y')
meta_dt <- meta_dt %>% mutate(Source3 = ifelse(is.na(Source_side)==T, Source2, paste0(Source_side,'_',Source2)))
meta_dt$deadbody <- factor(meta_dt$deadbody, levels =c('DB2','DB3','DB5','DB6','DB8','DB9','DB10'))

#Sourcing functions
source('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/TargetSeq_functions.R')

#load data
merged_m_vaf_tbl <- readRDS('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/merged_m_vaf_tbl.rds')



# --------


m_meta_dt <- meta_dt %>% mutate_at(vars(starts_with("SBS")), list(~ .*n_snv*0.01)) %>% mutate(SBS7total = SBS7a + SBS7b + SBS7c + SBS7d) %>%
  mutate(SBS5cor = SBS5 - (SBS7total*0.1087), SBS18cor = SBS18 + (SBS7total*0.0106)) %>%
  mutate(n_snv_endo = SBS1 + SBS5cor + SBS18cor)
f_meta_dt <- m_meta_dt %>% filter(Cell_type == 'skin_fb')
multi_ids <- f_meta_dt %>% group_by(tissue_id) %>% count() %>% filter(n>1) %>% pull(tissue_id)
f_meta_dt <- f_meta_dt %>% filter(tissue_id %in% multi_ids)
fm_meta_dt <- f_meta_dt %>% dplyr::select(tissue_id, sample_id, Source_class, Source_side, SBS7total, n_snv_endo) %>% gather(key='group', value='n_muts', 5:6)
x_order <- f_meta_dt %>% group_by(deadbody, tissue_id) %>% summarise(med7 = median(SBS7total)) %>% arrange(deadbody, desc(med7)) %>% pull(tissue_id)


g1 <- fm_meta_dt %>% mutate(tissue_id=factor(tissue_id, levels=x_order)) %>% 
  ggplot(aes(x=tissue_id, y=n_muts, fill=group, color=group))+
  geom_boxplot()+
  #scale_x_discrete(limits= x_order)+
  scale_fill_nejm()+
  scale_color_nejm()+
  theme_minimal()+
  ylab("Number of mutations") + 
  theme(
    legend.position="none",
    text=element_text(size=15),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    #axis.text.x=element_text(angle=90),
    #axis.line=element_line(color="black"),
    panel.border = element_rect(color="black", fill=NA))
  
g2 <- fm_meta_dt %>% dplyr::select(tissue_id:Source_side) %>% separate(tissue_id, into=c("DB", "Anatomic_site"), extra="drop", remove=F) %>% 
  dplyr::select(tissue_id, DB) %>% distinct() %>% 
  mutate(tissue_id=factor(tissue_id, x_order)) %>% 
  ggplot(aes(x=tissue_id, y=factor(1), fill=DB, color=DB)) +
  geom_bar(stat="identity", width=1.0) + 
  scale_y_discrete(expand=c(0,0))+
  scale_x_discrete(limits= x_order)+
  scale_fill_jama()+
  scale_color_jama()+
  theme_void() + 
  theme(
    legend.position="none",
    panel.border = element_rect(color="black", fill=NA))
  

g3 <- fm_meta_dt %>% dplyr::select(tissue_id,Source_class) %>% distinct() %>% 
  mutate(tissue_id=factor(tissue_id, x_order)) %>% 
ggplot(aes(x=tissue_id, y=factor(1), fill=Source_class, color=Source_class)) +
  #scale_fill_manual(values=c("HN"=rgb(41/255,61/255,67/255),
                             #"LE"=rgb(40/255,142/255,203/255),
                             #"trunk"=rgb(160/255,51/255,52/255),
                             #"UE"=rgb(104/255,161/255,133/255))) +
  #scale_color_manual(values=c("HN"=rgb(41/255,61/255,67/255),
                             #"LE"=rgb(40/255,142/255,203/255),
                             #"trunk"=rgb(160/255,51/255,52/255),
                             #"UE"=rgb(104/255,161/255,133/255))) +  
  scale_fill_jco()+
  scale_color_jco()+
  geom_bar(stat="identity", width=1.0) + 
  scale_y_discrete(expand=c(0,0))+
  theme_void() + 
  theme(
    legend.position="none",
    panel.border = element_rect(color="black", fill=NA))

g4 <- fm_meta_dt %>% dplyr::select(tissue_id:Source_side) %>% distinct() %>% 
  mutate(tissue_id=factor(tissue_id, x_order)) %>% 
  ggplot(aes(x=tissue_id, y=factor(1), fill=Source_side, color=Source_side)) +  
  geom_bar(stat="identity", width=1.0) + 
  scale_y_discrete(expand=c(0,0))+
  theme_void() + 
  theme(
    legend.position="none",
    panel.border = element_rect(color="black", fill=NA))

g6 <- fm_meta_dt %>% mutate(tissue_id=factor(tissue_id, levels=x_order)) %>% group_by(tissue_id) %>% mutate(n=1:n()) %>% 
  ggplot(aes(x=tissue_id, y=n)) + geom_point(size=0.3) + theme_minimal() + theme(text=element_text(size=10), axis.title.x=element_blank(), axis.text.x=element_blank()) 

plot_grid(g6, g2, g1, g3, ncol=1, rel_heights = c(0.20, 0.05, 0.70, 0.05), align="v")




fm_meta_dt %>% filter(grepl("6_", tissue_id)) %>% #View()
  #filter(tissue_id=="6_LLArm_SF") %>%
  filter(grepl("LLLA", sample_id)) %>% 
  mutate(group=factor(group, levels=c("SBS7total", "n_snv_endo"))) %>% 
           ggplot(aes(x=sample_id, y=n_muts, fill=group)) + geom_bar(stat='identity') + theme_minimal() + 
  ylab("Count") + 
  scale_fill_manual(values=c("n_snv_endo"=rgb(172/255,41/255,31/255), "SBS7total"=rgb(30/255,92/255,165/255))) +
  theme(text=element_text(size=15), legend.position="none", axis.title.x=element_blank(), axis.text.x=element_blank())

         