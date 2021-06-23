library(tidyverse)
library(ggsci)
library(scales)
library(ggtree)
library(ggpubr)
library(cowplot)

#load data
merged_m_vaf_tbl <- readRDS('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/merged_m_vaf_tbl.rds')
lng_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_lineage_210322.txt')
meta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_sample_210322.txt') %>% filter(current_final_for_lineage == 'Y')
f_vaf_tbl <- readRDS('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/f_vaf_tbl_rank.rds')

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


#################################################germ layer###################################################


#wilcox test per variant group
pval_co=0.05
this_db='DB6'

Wilcox_Test_on_Phylo_Ecto_MesoEndo <- function(f_vaf_tbl, this_db, lng_dt, nwk_tbl, pval_co){
  fm_vaf_tbl <- f_vaf_tbl %>% filter(deadbody == this_db)
  f_vaf_tbl <- fm_vaf_tbl %>% filter(dominant_layer2 %in% c('mesoderm','endoderm','ectoderm') & is.na(rank)== F)
  print(f_vaf_tbl %>% select(lineage_id2, rank))
  f_vaf_tbl <- f_vaf_tbl %>% mutate(var_g_id = paste0(lineage_id2,'_',rank))
  res_tbl = tibble(var_g_id = as.character(), wilcox_p_Eclow =as.numeric(), wilcox_p_Echi = as.numeric(), wilcox_p_both = as.numeric(),
                   ecto_mean = as.numeric(), mesoendo_mean = as.numeric(), ecto_n = as.integer(), mesoendo_n = as.integer())
  n=0
  
  for(t_varid in unique(f_vaf_tbl$var_g_id)){  #group id
    n=n+1
    a <- f_vaf_tbl %>% filter(var_g_id == t_varid & dominant_layer2 == 'ectoderm') %>% pull(VAF)
    b <- f_vaf_tbl %>% filter(var_g_id == t_varid & dominant_layer2 != 'ectoderm') %>% pull(VAF)
    if (length(a) > 0 & length(b) > 0 ){
      res1 <- wilcox.test(a,b, alternative = 'less')
      pval1 <- res1$p.value
      res2 <- wilcox.test(a,b, alternative = 'greater')
      pval2 <- res2$p.value
      res3 <- wilcox.test(a,b, alternative = 'two.sided')
      pval3 <- res3$p.value
      
    } else {
      pval1 <- NA
      pval2 <- NA
      pval3 <- NA
    }
    res_tbl[n,] <- list(t_varid, pval1, pval2, pval3, mean(a), mean(b), length(a), length(b))
  }
  res_tbl <- res_tbl %>% mutate(wilcox_p = ifelse(wilcox_p_Eclow < wilcox_p_Echi, wilcox_p_both, -1*wilcox_p_both), discrete = floor(wilcox_p*100))
  res_tbl2 <- res_tbl %>% select(var_g_id, wilcox_p)
  res_tbl <- res_tbl %>% mutate(discrete = as.factor(discrete)) %>% select(var_g_id, discrete)
  
  #wilcox pvalue on phylogenetic tree
  if(T){ #order in rank (random rank in ties)
    m_res_tbl <- left_join(res_tbl,
                           f_vaf_tbl %>% filter(deadbody == this_db & time_of_lineage == 'early' & is.na(rank) == F ) %>% select(lineage_id2, var_id, rank, rank_ctie_first) %>% unique() %>% 
                             mutate(var_g_id = paste0(lineage_id2, '_', rank))
    )
    m_res_tbl <- m_res_tbl %>% mutate(rank2 = paste0('VAF', rank_ctie_first)) %>% select(-var_id, -var_g_id, -rank_ctie_first, -rank) %>% spread(key=rank2, value=discrete)
    max_VAF_len =20
    for(id in setdiff(paste0('VAF',1:max_VAF_len), colnames(m_res_tbl))){
      m_res_tbl <- m_res_tbl %>% mutate(this = NA)
      colnames(m_res_tbl) <- gsub('this',id,colnames(m_res_tbl))
    }
  }
  
  if(T){ #continuous variable
    library(colorRamps)
    total_break=50
    col_breaks = seq(0, 0.5, length.out=total_break)
    
    a <- length(col_breaks[col_breaks < 0.05] )  #change color at 1
    b <- length(col_breaks[col_breaks < 0.1])-a  #change color at 5
    #rampCol1 <- colorRampPalette(c("red", "orange"))(n=a) #229 167 13
    #rampCol2 <- colorRampPalette(c("orange","gray"))(n=b)
    #rampCol3 <- colorRampPalette(c("gray", "gray"))(n=total_break-(a+b)-1)
    #my_palette1<- c(rampCol1, rampCol2, rampCol3)
    
    #rampCol1 <- colorRampPalette(c(rgb(28/255,70/255,166/255), "gray"))(n=a+b)
    rampCol1 <- colorRampPalette(c(rgb(126/255,71/255,142/255), "gray"))(n=a+b) 
    rampCol2 <- colorRampPalette(c("gray", "gray"))(n=total_break-(a+b)-1)
    my_palette1 <- c(rampCol1, rampCol2)
    
    #legend..................................
    if(T){
      plot.new()
      plot(c(0,2),c(0,20),type = 'n', axes = F,xlab = '', ylab = '', main = 'Pvalues from Meso-Endoderm > Ectoderm')
      colfunc1 <- rev(rampCol1)
      legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
      rasterImage(legend_image1, 0, 0, 0.3,5)
      #colfunc2 <- rev(rampCol2)
      #legend_image2 <- as.raster(matrix(colfunc2, ncol=1))
      #rasterImage(legend_image2, 0, 5, 0.3,10)
      text(x=0.5, y = c(0,5,10,20), labels =c(0,0.05, 0.1, 0.2), cex=2)
    }
    
    total_break=50
    col_breaks = seq(0, 0.5, length.out=total_break)
    
    a <- length(col_breaks[col_breaks < 0.05] )  #change color at 1
    b <- length(col_breaks[col_breaks < 0.1])-a  #change color at 5
    #rampCol1 <- colorRampPalette(c(pal_combi[4], pal_combi[2]))(n=a)
    #rampCol2 <- colorRampPalette(c(pal_combi[2],"gray"))(n=b)
    #rampCol3 <- colorRampPalette(c("gray", "gray"))(n=total_break-(a+b)-1)
    #my_palette2<- c(rampCol1, rampCol2, rampCol3)
    #my_palette2<- colorRampPalette(c(rgb(28/255,70/255,166/255), "gray"))(n=total_break-1)
    
    rampCol1 <- colorRampPalette(c(rgb(229/255,167/255,13/255), "gray"))(n=a+b)   
    rampCol2 <- colorRampPalette(c("gray", "gray"))(n=total_break-(a+b)-1)
    my_palette2 <- c(rampCol1, rampCol2)
    
    #legend..................................
    if(T){
      plot.new()
      plot(c(0,2),c(0,20),type = 'n', axes = F,xlab = '', ylab = '', main = 'Pvalues from Meso-Endoderm < Ectoderm')
      colfunc1 <- rev(rampCol1)
      legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
      rasterImage(legend_image1, 0, 0, 0.3,5)
      colfunc2 <- rev(rampCol2)
      legend_image2 <- as.raster(matrix(colfunc2, ncol=1))
      rasterImage(legend_image2, 0, 5, 0.3,10)
      colfunc3 <- rev(rampCol3)
      legend_image3 <- as.raster(matrix(colfunc3, ncol=1))
      rasterImage(legend_image3, 0, 10, 0.3,20)
      text(x=0.5, y = c(0,5,10,50), labels =c(0,0.05, 0.1, 0.2), cex=2)
    }
    
    my_palette <- c(rev(my_palette2), my_palette1)
    names(my_palette) <- seq(-49,48,1)
  }
  
  #assign point changing to <0.05 and >= 0.05
  if(T){
    m_res_tbl2 <- left_join(res_tbl2,
                            f_vaf_tbl %>% filter(deadbody == this_db & time_of_lineage == 'early' & is.na(rank) == F ) %>% select(lineage_id2, var_id, rank, n_vars_same_rank, start_mtn) %>% unique() %>% 
                              mutate(var_g_id = paste0(lineage_id2, '_', rank)))
    
    el_ids <- lng_dt %>% filter(deadbody == this_db & time_of_lineage == 'early_to_late') %>% pull(lineage_id2)
    changing_lines = tibble()
    for(t_id in el_ids){
      find_v = 'N'
      tmp_v <- unlist(strsplit(t_id, '-'))
      t_id_list <- c()
      for (i in 1:length(tmp_v)){
        t_id_list <- c(t_id_list, paste(tmp_v[1:i], collapse="-"))
      }
      for(n in 1:(length(t_id_list)-1)){
        rank_list <- m_res_tbl2 %>% filter(lineage_id2 == t_id_list[n]) %>% pull(rank) %>% unique() %>% sort()
        if(length(rank_list) == 0){
          next
        }
        for (m in 1:length(rank_list)){
          t_pv <- m_res_tbl2 %>% filter(lineage_id2 == t_id_list[n] & rank == rank_list[m]) %>% pull(wilcox_p) %>% unique() %>% abs()
          if(is.na(t_pv)){
            next
          }
          if (t_pv < pval_co){
            changing_lines <- bind_rows(changing_lines, m_res_tbl2 %>% filter(lineage_id2 == t_id_list[n] & rank == rank_list[m]) %>% select(-var_id) %>% unique())
            find_v = 'Y'
            break
          }
        }
        if (find_v == 'Y'){
          break
        }
      }
    }
    changing_lines <- changing_lines %>% unique()
    changing_lines <- changing_lines %>% mutate(point0.05 = (start_mtn + rank - 1 + start_mtn + rank - 1 + n_vars_same_rank)/2)
  }
  my_func <- function(num){
    if (is.na(num)){
      return (NA)
    } else if(num < -49){
      return (-49)
    } else if ( num > 48){
      return (48) 
    } else {
      return (num)
    }
  }
  m_res_tbl <- m_res_tbl %>% mutate_at(vars(starts_with('VAF')), as.character) %>% mutate_at(vars(starts_with('VAF')), as.numeric) %>%
    rowwise() %>% mutate_at(vars(starts_with('VAF')),my_func) %>% mutate_at(vars(starts_with('VAF')), as.factor)
  m_lng_dt <- left_join(lng_dt %>% filter(deadbody == this_db) ,m_res_tbl)
  m_lng_dt <- left_join(m_lng_dt, changing_lines %>% select(lineage_id2, point0.05))
  m_lng_dt <- m_lng_dt %>% dplyr::rename(taxa = lineage_id) %>% select(-deadbody)
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
  g <- ggtree(tree) %<+% m_lng_dt + theme_tree2()+
    coord_cartesian(xlim = c(-1,30))+
    geom_segment2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt, xend = x -n_pointmt +1, y = y, yend = y,color = VAF1) ,size=2)+
    geom_segment2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt+1, xend = x -n_pointmt +2, y = y, yend = y, color = VAF2), size=2)+
    geom_segment2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt+2, xend = x -n_pointmt +3, y = y, yend = y, color = VAF3), size=2)+
    geom_segment2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt+3, xend = x -n_pointmt +4, y = y, yend = y, color = VAF4), size=2)+
    geom_segment2(aes(subset = is.na(VAF5) == F, x=x-n_pointmt+4, xend = x -n_pointmt +5, y = y, yend = y, color = VAF5), size=2)+
    geom_segment2(aes(subset = is.na(VAF6) == F, x=x-n_pointmt+5, xend = x -n_pointmt +6, y = y, yend = y, color = VAF6), size=2)+
    geom_segment2(aes(subset = is.na(VAF7) == F, x=x-n_pointmt+6, xend = x -n_pointmt +7, y = y, yend = y, color = VAF7), size=2)+
    geom_segment2(aes(subset = is.na(VAF8) == F, x=x-n_pointmt+7, xend = x -n_pointmt +8, y = y, yend = y, color = VAF8), size=2)+
    geom_segment2(aes(subset = is.na(VAF9) == F, x=x-n_pointmt+8, xend = x -n_pointmt +9, y = y, yend = y, color = VAF9), size=2)+
    geom_segment2(aes(subset = is.na(VAF10) == F, x=x-n_pointmt+9, xend = x -n_pointmt +10, y = y, yend = y, color = VAF10), size=2)+
    geom_segment2(aes(subset = is.na(VAF11) == F, x=x-n_pointmt+10, xend = x -n_pointmt +11, y = y, yend = y, color = VAF11), size=2)+
    geom_segment2(aes(subset = is.na(VAF12) == F, x=x-n_pointmt+11, xend = x -n_pointmt +12, y = y, yend = y, color = VAF12), size=2)+
    geom_segment2(aes(subset = is.na(VAF13) == F, x=x-n_pointmt+12, xend = x -n_pointmt +13, y = y, yend = y, color = VAF13), size=2)+
    geom_segment2(aes(subset = is.na(VAF14) == F, x=x-n_pointmt+13, xend = x -n_pointmt +14, y = y, yend = y, color = VAF14), size=2)+
    geom_segment2(aes(subset = is.na(VAF15) == F, x=x-n_pointmt+14, xend = x -n_pointmt +15, y = y, yend = y, color = VAF15), size=2)+
    geom_segment2(aes(subset = is.na(VAF16) == F, x=x-n_pointmt+15, xend = x -n_pointmt +16, y = y, yend = y, color = VAF16), size=2)+
    geom_segment2(aes(subset = is.na(VAF17) == F, x=x-n_pointmt+16, xend = x -n_pointmt +17, y = y, yend = y, color = VAF17), size=2)+
    geom_segment2(aes(subset = is.na(VAF18) == F, x=x-n_pointmt+17, xend = x -n_pointmt +18, y = y, yend = y, color = VAF18), size=2)+
    geom_segment2(aes(subset = is.na(VAF19) == F, x=x-n_pointmt+18, xend = x -n_pointmt +19, y = y, yend = y, color = VAF19), size=2)+
    geom_segment2(aes(subset = is.na(VAF20) == F, x=x-n_pointmt+19, xend = x -n_pointmt +20, y = y, yend = y, color = VAF20), size=2)+
    #geom_point(aes(x=point0.05), shape=1)+
    scale_color_manual(values = my_palette)+
    theme(axis.text = element_text(size=15, angle=90), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"),
          legend.position="none", axis.line.x=element_blank())
  
  print(g)
}





ff_vaf_tbl <- f_vaf_tbl %>% filter(time_of_lineage == 'early' & tissue_class !='cancer_tissue')
pval_co=0.05
this_db='DB6'
Wilcox_Test_on_Phylo_Ecto_MesoEndo(ff_vaf_tbl, this_db, lng_dt, nwk_tbl, pval_co) # Main Figure 3b





#Ectoderm-mesoendoderm bias timing curve!!!
m_ch_dt <- tibble()
for (this_db in c('DB6','DB8','DB9','DB10')){
  changing_lines <- Wilcox_Test_on_Phylo_Ecto(merged_m_vaf_tbl, this_db, lng_dt, nwk_tbl, pval_co)
  m_ch_dt <- bind_rows(m_ch_dt, changing_lines)
}

tmp_tbl <- left_join(
  m_ch_dt,
  lng_dt %>% select(deadbody, lineage_id2, early_desc_n)
)
tmp_tbl <- tmp_tbl %>% mutate(new_mtn = start_mtn + rank) %>% group_by(deadbody, new_mtn) %>% summarise(added_lineage_n = sum(early_desc_n))
tmp_tbl <- tmp_tbl %>% mutate(cum_sum = cumsum(added_lineage_n))
tmp_tbl <- bind_rows(tmp_tbl, tibble(deadbody = tmp_tbl$deadbody %>% unique(), new_mtn = 0,added_lineage_n = 0, cum_sum=0))
tmp_tbl <- left_join(tmp_tbl, lng_dt %>% filter(time_of_lineage== 'early_to_late') %>% group_by(deadbody) %>% summarise(total_lineage_n=n()))
tmp_tbl <- tmp_tbl %>% ungroup() %>%  mutate(deadbody=factor(deadbody, levels=c("DB6", "DB8", "DB9", "DB10")))

ggplot(tmp_tbl, aes(x=new_mtn, y=cum_sum/total_lineage_n, fill=deadbody, color=deadbody))+
  geom_line(aes(group = deadbody), size=0.5)+
  geom_point(size=3, pch=21, stroke=0.5, color="black") +
  xlab("Number of mutations")+ylab("Proportion of biased lineages")+
  scale_y_continuous(limits=c(0,1))+
  scale_x_continuous(limits=c(0,20))+
  scale_fill_manual(values=c("DB2"=rgb(203/255,104/255,41/255),
                              "DB3"=rgb(36/255,122/255,191/255),
                              "DB5"=rgb(141/255,34/255,40/255),
                              "DB6"=rgb(87/255,145/255,114/255),
                              "DB8"=rgb(69/255,60/255,116/255),
                              "DB9"=rgb(89/255,83/255,70/255),
                              "DB10"=rgb(32/255,47/255,52/255)))+ 
  scale_color_manual(values=c("DB2"=rgb(203/255,104/255,41/255),
                             "DB3"=rgb(36/255,122/255,191/255),
                             "DB5"=rgb(141/255,34/255,40/255),
                             "DB6"=rgb(87/255,145/255,114/255),
                             "DB8"=rgb(69/255,60/255,116/255),
                             "DB9"=rgb(89/255,83/255,70/255),
                             "DB10"=rgb(32/255,47/255,52/255))) +
  theme_minimal() +
  theme(text=element_text(size=15), legend.title=element_blank(), legend.key=element_blank()) 
  


