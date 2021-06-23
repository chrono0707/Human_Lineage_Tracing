#load libraries
library(tidyverse)
library(ggtree)
library(ggsci)
library(scales)
library(cowplot)
library(ggrepel)

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
merged_m_vaf_tbl <- readRDS('/home/users/sypark/00_Project/06_LineageTracing/merged_data/Target_seq_related/merged_m_vaf_tbl.rds')

#first_occuring_variants
first_variants <- tribble(~deadbody, ~lineage_id2, ~var_id,
                          'DB3','DB3_L1','1_74024235_G_T',
                          'DB3','DB3_L1','2_4923793_T_C',
                          'DB3','DB3_L1','7_32925447_G_C',
                          'DB3','DB3_L1','7_129376347_C_T',
                          'DB3','DB3_L1','4_34076338_C_A',
                          'DB3','DB3_L1','16_85372995_G_A',
                          'DB3','DB3_L2','19_19279396_T_G',
                          'DB3','DB3_L2','18_6114473_T_C',
                          'DB6','DB6_L1','13_61660346_G_A',
                          'DB6','DB6_L2','2_82379165_T_A',
                          'DB6','DB6_L2','16_16266502_G_A',
                          'DB8','DB8_L1','5_85367961_G_T',
                          'DB8','DB8_L1','8_123054768_C_T',
                          'DB8','DB8_L1','13_102309939_G_A',
                          'DB8','DB8_L1','22_17072858_C_G', 
                          'DB8','DB8_L2','1_162486606_C_A',
                          'DB9','DB9_L1','11_50359839_A_C',
                          'DB9','DB9_L2','8_58127053_G_A',
                          'DB10','DB10_L1','1_166645826_C_T',
                          'DB10','DB10_L1','4_1293246_C_T', 
                          'DB10','DB10_L1','8_93863618_C_A', 
                          'DB10','DB10_L1','15_94394073_C_A',
                          'DB10','DB10_L2','1_56148328_C_A',
                          'DB10','DB10_L2','6_44692877_C_A', 
                          'DB10','DB10_L2','8_139121031_C_T',
                          'DB10','DB10_L2','9_94617243_T_C', 
                          'DB10','DB10_L2','11_96493593_G_T',
                          'DB10','DB10_L2','15_62973855_A_C'
)

#---------------------------------------------------------------script-------------------------------------------------------------------
GSNP_mVAF <- merged_m_vaf_tbl %>% filter(grepl('GSNP', lineage_id2)) %>% pull(VAF) %>% mean(., na.rm=T)

p_list <- list()
n=0
for (this_db in c('DB3','DB6','DB8','DB9','DB10')){
  n=n+1
  m_vaf_tbl <- merged_m_vaf_tbl %>% filter(deadbody == this_db & tissue_class != 'cancer_tissue')
  fm_vaf_tbl <- m_vaf_tbl %>% filter(var_id %in% first_variants$var_id[first_variants$deadbody == this_db]) %>% filter(DP > 100)
  tmp_tbl <- fm_vaf_tbl %>% group_by(sample_id, lineage_id2) %>% summarise(meanVAF = mean(VAF)) 
  m_vaf_tbl %>% filter(lineage_id2 == paste0(this_db,'_GSNP')) %>% filter(DP > 100) %>% dplyr::select(sample_id, lineage_id2, VAF)
  tmp_tbl <- left_join(tmp_tbl,pmeta_dt %>% dplyr::select(sample_id, dominant_layer, anatomy_class1)) %>% separate(lineage_id2, c('deadbody','lineage_id'), sep="_") %>% 
    spread(key=lineage_id, value=meanVAF) %>% mutate(L1ratio=L1/(L1+L2), L2ratio=L2/(L1+L2)) %>% dplyr::select(-L1, -L2) %>% filter(is.na(L1ratio)==F & is.na(L2ratio) == F) %>%
    gather(-c(sample_id,deadbody,dominant_layer, anatomy_class1), key=lineage_id, value=ratio)
  
  if(this_db=="DB8") x_order <- tmp_tbl %>% filter(lineage_id == 'L1ratio') %>% arrange(ratio) %>% pull(sample_id)
  else               x_order <- tmp_tbl %>% filter(lineage_id == 'L1ratio') %>% arrange(desc(ratio)) %>% pull(sample_id)
  
  exp_v1 <- m_vaf_tbl %>% filter(lineage_id2 == paste0(this_db, '_L2')) %>% pull(expected_VAF) %>% unique()
  exp_v2 <- m_vaf_tbl %>% filter(lineage_id2 == paste0(this_db, '_L1')) %>% pull(expected_VAF) %>% unique()
  tmp_tbl %<>% mutate(ratio=ratio*100) 
  
  if(this_db=="DB8") {
    tmp_tbl %<>% mutate(lineage_id=factor(lineage_id, levels=c("L2ratio", "L1ratio"))) 
    #tmp_tbl %<>% ungroup() %>% mutate(sample_id=factor(sample_id, levels=rev(unique(sample_id))))
    yinter <- exp_v2/(exp_v1+exp_v2)*100
  } else {
    yinter <- exp_v1/(exp_v1+exp_v2)*100
  }

  p_list[[n]] <-ggplot(tmp_tbl, aes(x=sample_id, y=ratio, fill=lineage_id, color=lineage_id))+
    geom_bar(stat="identity", width=1.0)+
    geom_hline(yintercept=yinter, linetype="dashed")+
    scale_x_discrete(limits=x_order, expand=c(0.01,0.01))+
    #scale_x_discrete(limits=x_order, expand=c(0.01,0.01))+
    scale_y_continuous(expand=c(0.00,0.00))+
    #scale_y_continuous(expand=c(0.01,0.01))+ 
    scale_fill_simpsons() + 
    scale_color_simpsons() + 
    theme_void() + 
    #xlab("Proportion") +
    ylab("Proportion") + 
    coord_flip() + 
    theme(#axis.text.y = element_text(size=15), 
          axis.text.x = element_text(size=20), 
          axis.title.x = element_text(size=20),  
          #axis.title.y = element_text(size=20, angle=90),
          legend.position = 'none',
          plot.margin = unit(c(0,0.0,0,0.0), "cm"),
          panel.spacing = unit(c(0,0,0,0), "cm"),  
          #axis.ticks.y = element_line(color="black"),
          axis.ticks.x = element_line(color="black"),
          axis.ticks.length=unit(.05, "cm"),
          panel.border=element_rect(color='black', fill=NA, size=2) 
          #panel.background = element_rect(color="black", size=1)
          )
  
  p_list[[n]] <- plot_grid(p_list[[n]], 
                           tmp_tbl %>% 
                             ggplot(aes(x=sample_id, y=factor(1), fill=anatomy_class1)) + geom_tile() + theme_void() + 
                             scale_fill_manual(values=c("HN"=rgb(42/255,61/255,68/255),
                                                        "internal_organ"=rgb(214/255,124/255,52/255),
                                                        "LE"=rgb(40/255,142/255,203/255),
                                                        "trunk"=rgb(160/255,51/255,52/255),
                                                        "UE"=rgb(104/255,161/255,133/255))) + 
                             scale_x_discrete(expand=c(0.01,0.01)) + 
                             scale_y_discrete(expand=c(0,0)) + 
                             coord_flip() +
                             theme(legend.position="none", plot.margin=unit(c(0,0,0,0.1), unit="cm")),
                           tmp_tbl %>% 
                             ggplot(aes(x=sample_id, y=factor(1), fill=dominant_layer)) + geom_tile() + theme_void() + 
                             #scale_fill_igv() + 
                             scale_fill_manual(values=c("ectoderm"=rgb(62/255,48/255,255/255),
                                                        "mesoderm"=rgb(193/255,40/255,38/255),
                                                        "mixed"=rgb(98/255,140/255,70/255))) +
                             scale_x_discrete(expand=c(0.01,0.01)) + 
                             scale_y_discrete(expand=c(0,0)) + 
                             coord_flip() +
                             theme(legend.position="none", plot.margin=unit(c(0,0,0,0.1), unit="cm")),                           
                           nrow=1,
                           align="hv", rel_widths = c(0.90,0.05, 0.05)
                           )
  
  print(exp_v1/(exp_v1+exp_v2))
}

p_list[[1]]

