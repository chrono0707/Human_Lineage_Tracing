###########################<Common processing>##############################################
#load library
library(tidyverse)
library(ggsci)
library(ggtree)
library(ggrepel)
library(cowplot)

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

library(pheatmap)
library(grid)
data_path=tribble(~deadbody, ~path,
                  'DB2','/home/users/sypark/00_Project/06_LineageTracing/db2/10_point_mutations/perLineage/DB2_early_pointmts.txt.27sCall',
                  'DB3','/home/users/sypark/00_Project/06_LineageTracing/db3/14_point_mutations/perLineage/DB3_early_pointmts.txt.43sCall',
                  'DB5','/home/users/sypark/00_Project/06_LineageTracing/db5/11_point_mutations/perLineage/DB5_early_pointmts.txt.28sCall',
                  'DB6','/home/users/sypark/00_Project/06_LineageTracing/db6/12_point_mutations/perLineage/DB6_early_pointmts.txt.115sCall',
                  'DB8','/home/users/sypark/00_Project/06_LineageTracing/db8/09_point_mutations/perLineage/DB8_early_pointmts.txt.47sCall',
                  'DB9','/home/users/sypark/00_Project/06_LineageTracing/db9/09_point_mutations/perLineage/DB9_early_pointmts.txt.35sCall',
                  'DB10','/home/users/sypark/00_Project/06_LineageTracing/db10/09_point_mutations/perLineage/DB10_early_pointmts.txt.39sCall'
)

plot_list=list()

scoreCol = function(x) {
  score = 0
  for (i in 1:length(x)) {
    if (x[i] != 0) {
      score = score + 2^(length(x) - i * 1/x[i])
    }
  }
  return(score)
}

for (this_db in c('DB3','DB6','DB8','DB9','DB10')){
  print(this_db)
  dt <- read_tsv(data_path$path[data_path$deadbody == this_db])
  dt <- dt %>% mutate(var_id = paste0(`#CHROM`,':',POS, ' ',REF,'>',ALT)) %>% select(var_id, lineage_id2, ends_with('vafpct'), -starts_with('blood'))
  dt <- left_join(dt,lng_dt %>% select(lineage_id2,n_samples))
  
  f_dt <- dt #%>% filter(n_samples >=3)
  colnames(f_dt) <- gsub('_vafpct','',colnames(f_dt))
  s_order <- meta_dt %>% filter(deadbody == this_db) %>% arrange(lineage_id_cindel) %>% pull(sample_id)
  f_dt <- f_dt %>% arrange(lineage_id2) %>% select(var_id, s_order)
  
  mx <- f_dt %>% as.data.frame() %>%column_to_rownames('var_id') %>% as.matrix()
  #plot_list[[this_db]] <- pheatmap(mx, cluster_cols = F, cluster_rows = F, show_rownames = F, show_colnames = F,main = this_db)
  #plot_list[[this_db]] <-
  
  mx <- mx[order(rowSums(mx!=0), decreasing=T),]
  mx[is.na(mx)] <- 0
  scores = apply(mx!=0, 2, scoreCol)
  
  n_row <- 50
  if(this_db=="DB6") n_row <- 100
  plot_list[[this_db]] <- Heatmap(mx[1:n_row,], col=rev(brewer.pal(10,"Spectral")), 
          column_order = order(scores, decreasing = TRUE),
          cluster_rows=F, cluster_columns = F,
          show_column_names = F, show_row_names = F,
          cell_fun = function(j,i,x,y,width,height,fill){
            grid.rect(x=x, y=y, width=width, height=height,gp = gpar(col=fill,fill=fill))
          }, name="VAF", border="black")
 
 
}




grid.newpage()
pushViewport(viewport(x=0, y=0.5, width = 0.3, height = 0.5, just = c('left', 'bottom')))
print(plot_list[['DB3']], newpage = F)
popViewport(1)
pushViewport(viewport(x=0.3, y=0.5, width = 0.3, height = 0.5, just = c('left', 'bottom')))
print(plot_list[['DB8']], newpage = F)
popViewport(1)
pushViewport(viewport(x=0.6, y=0, width = 0.4, height = 1, just = c('left', 'bottom')))
print(plot_list[['DB6']], newpage = F)
popViewport(1)
pushViewport(viewport(x=0, y=0, width = 0.3, height = 0.5, just = c('left', 'bottom')))
print(plot_list[['DB9']], newpage = F)
popViewport(1)
pushViewport(viewport(x=0.3, y=0, width = 0.3, height = 0.5, just = c('left', 'bottom')))
print(plot_list[['DB10']], newpage = F)
popViewport(1)

