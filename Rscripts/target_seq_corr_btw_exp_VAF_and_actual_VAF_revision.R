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
merged_m_vaf_tbl <- readRDS('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/merged_m_vaf_tbl.rds')


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

#Draw correlation between expected VAF and actual VAF
x_order <- fm_vaf_tbl %>% arrange(expected_VAF) %>% pull(expected_VAF) %>% unique()
  fm_vaf_tbl <- merged_m_vaf_tbl %>% filter(DP >= 100 & time_of_lineage == 'early')
  fm_vaf_tbl %>% mutate(deadbody=factor(deadbody, levels=c("DB3", "DB6", "DB8", "DB9", "DB10"))) %>% 
    #mutate(expected_VAF=factor(expected_VAF, levels=x_order)) %>% 
    ggplot(aes(x=expected_VAF, y=VAF))+
    geom_jitter(aes(color=deadbody), alpha=0.3, size=1, width=0.005)+
    #geom_violin(aes(fill=lineage_id2), position=position_dodge2(padding=0)) + 
    geom_boxplot(aes(fill=lineage_id2), position="dodge", width=0.01, alpha=0, color="white", size=0.3) + 
    annotate("segment", x=0, xend=0.6, y=0, yend=0.6, color="red") +
    scale_color_manual(values=c("DB3"=rgb(40/255,142/255,203/255),
                               "DB6"=rgb(104/255,161/255,133/255),
                               "DB8"=rgb(87/255,80/255,135/255),
                               "DB9"=rgb(108/255,102/255,88/255),
                               "DB10"=rgb(42/255,61/255,68/255))) + 
    theme_minimal() +
    #scale_x_continuous(expand=c(0,0), limits=c(0,0.67)) + 
    scale_y_continuous(expand=c(0,0), limits=c(0,0.67))+
    ylab("VAF in bulk tissues") + xlab("Expected VAF") +
    theme(legend.position="none", 
          legend.title=element_blank(), 
          text=element_text(size=15), 
          axis.line.x=element_line(color="black"),
          axis.line.y=element_line(color="black")) 

