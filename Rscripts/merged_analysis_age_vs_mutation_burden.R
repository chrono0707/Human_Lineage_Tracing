###########################<Common processing>##############################################
#load library
library(tidyverse)
library(ggsci)
library(ggtree)
library(ggrepel)
library(cowplot)
library(binom)

#load dataset
meta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_sample_200908.txt') %>% filter(current_final_for_lineage == 'Y')
lng_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_lineage_200908.txt')
DB_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_DB_200522.txt')

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

#color setting
library(scales)
#show_col(pal_npg('nrc')(10))
db_pal <- pal_npg('nrc')(10)[c(1:5,9,7)]
names(db_pal) <- unique(meta_dt$deadbody)
cell_pal <- pal_npg('nrc')(6)
names(cell_pal) <- unique(meta_dt$Cell_type)
source2_pal <- pal_d3('category10')(6)
names(source2_pal) <- c('HN','lt_LE','lt_UE','rt_LE','trunk','rt_UE')
timing_pal <- pal_npg('nrc')(10)[c(1,3,4)]
names(timing_pal) <- c('early','early_to_late','late')

#Fx source
source('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/merged_analysis_functions.R')


# --------

m_meta_dt <- meta_dt %>% mutate_at(vars(starts_with("SBS")), list(~ .*n_snv*0.01)) %>% mutate(SBS7total = SBS7a + SBS7b + SBS7c + SBS7d)
m_meta_dt$Cell_type %>% table()
f_meta_dt <- m_meta_dt %>% filter(Cell_type == 'skin_fb')


#total SNV
ggplot(f_meta_dt, aes(x=age, y=n_snv))+
  geom_boxplot(aes(fill = deadbody))+
  geom_smooth(method="lm")+
  scale_x_continuous(limits=c(0,100))+
  scale_fill_manual(values=db_pal)+
  theme_syp+
  ggtitle('total SNV')
lm(f_meta_dt$n_snv ~ f_meta_dt$age) %>% summary


#v3_1,5,18
ggplot(f_meta_dt, aes(x=age, y=SBS1 + SBS5 + SBS18))+
  geom_boxplot(aes(fill = deadbody))+
  geom_smooth(method="lm")+
  scale_x_continuous(limits=c(0,100))+
  scale_y_continuous(limits=c(0,4500))+
  scale_fill_manual(values=db_pal)+
  theme_syp+
  ggtitle('sig1+5+18')

f_meta_dt <- f_meta_dt %>% mutate(v3_1518 = v3_1 + v3_5 + v3_18)
lm(f_meta_dt$v3_1518 ~ f_meta_dt$age) %>% summary()

#signature count correction
f_meta_dt <- f_meta_dt %>% mutate(SBS5cor = SBS5 - (SBS7total*0.1087)) %>% 
  mutate(SBS18cor = SBS18 + (SBS7total*0.0106))


ggplot(f_meta_dt, aes(x=age, y=SBS1 + SBS5cor + SBS18cor))+
  geom_jitter(width=1, color='gray80') + 
  geom_boxplot(aes(fill = deadbody), outlier.shape = NA)+
  geom_smooth(method="lm", color="red", size=0.3, fill="red", alpha=0.2)+
  scale_fill_manual(values=c("DB2"=rgb(203/255,104/255,41/255),
                             "DB3"=rgb(36/255,122/255,191/255),
                             "DB5"=rgb(141/255,34/255,40/255),
                             "DB6"=rgb(87/255,145/255,114/255),
                             "DB8"=rgb(69/255,60/255,116/255),
                             "DB9"=rgb(89/255,83/255,70/255),
                             "DB10"=rgb(32/255,47/255,52/255))) +
  theme_minimal() + 
  xlab("Age (years old)") + 
  ylab("Corrected endogenous\nmutation burden per sample") +
  theme(
    legend.position="bottom",
    text=element_text(size=15), 
    legend.title=element_blank(),
    panel.border = element_rect(color="black", fill=NA))
  
f_meta_dt <- f_meta_dt %>% mutate(v3_15c18c = v3_1 + v3_5cor + v3_18cor)
lm(f_meta_dt$v3_15c18c ~ f_meta_dt$age) %>% summary()
