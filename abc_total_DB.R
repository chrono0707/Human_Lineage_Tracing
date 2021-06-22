df_sim_DB3  <- getSimResults("/home/users/chrono0707/projects/02_DB/01_DB3/14_Simulations/01_initial_screening")
df_sim_DB10 <- getSimResults("/home/users/chrono0707/projects/02_DB/07_DB10/14_Simulation/01_initial_screening/")
df_sim_DB6  <- getSimResults("/home/users/chrono0707/projects/02_DB/02_DB6/14_Simulation/03_initial_screening_narrowed_at/")
df_sim_DB8  <- getSimResults("/home/users/chrono0707/projects/02_DB/05_DB8/14_simulation/02_narrowed_range/")
df_sim_DB9  <- getSimResults("/home/users/chrono0707/projects/02_DB/06_DB9/14_Simulation/01_initial_screening/")


temp_DB3 <- df_sim_DB3 %>% 
  filter(number_of_1st_lineage==2) %>% 
  separate(contribution_of_1st_lineage, c("lineage1", "lineage2")) %>% 
  separate(branch_length_of_1st_lineage, c("shared1", "shared2")) %>% 
  mutate(lin1=ifelse(lineage1>=lineage2, lineage1, lineage2), lin2=ifelse(lineage1>=lineage2, lineage2, lineage1)) %>% 
  mutate(sha1=ifelse(lineage1>=lineage2, shared1,  shared2),  sha2=ifelse(lineage1>=lineage2, shared2,  shared1)) %>%
  mutate(lin1=as.numeric(lin1), lin2=as.numeric(lin2), sha1=as.numeric(sha1), sha2=as.numeric(sha2)) %>% 
  mutate(distance=((number_of_shared_mutations - target_DB3["number_of_shared_mutations"]))^2 + 
           ((multifurcation_score       - target_DB3["multifurcation_score"]))^2 + 
           sqrt(((lin1-target_DB3["lin1"]))^2 + ((lin2-target_DB3["lin2"]))^2 + ((sha1-target_DB3["sha1"]))^2 + ((sha2-target_DB3["sha2"]))^2)^2) %>% 
  arrange(distance) %>% dplyr::slice(1:50)


temp_DB10 <- df_sim_DB10 %>% 
  filter(number_of_1st_lineage==2) %>% #filter(at==3 & n==3) %>% 
  separate(contribution_of_1st_lineage, c("lineage1", "lineage2")) %>% 
  separate(branch_length_of_1st_lineage, c("shared1", "shared2")) %>% 
  mutate(lin1=ifelse(lineage1>=lineage2, lineage1, lineage2), lin2=ifelse(lineage1>=lineage2, lineage2, lineage1)) %>% 
  mutate(sha1=ifelse(lineage1>=lineage2, shared1,  shared2),  sha2=ifelse(lineage1>=lineage2, shared2,  shared1)) %>%
  mutate(lin1=as.numeric(lin1), lin2=as.numeric(lin2), sha1=as.numeric(sha1), sha2=as.numeric(sha2)) %>% 
  mutate(distance=((number_of_shared_mutations - target_DB10["number_of_shared_mutations"]))^2 + 
           ((multifurcation_score       - target_DB10["multifurcation_score"]))^2 + 
           sqrt(((lin1-target_DB10["lin1"]))^2 + ((lin2-target_DB10["lin2"]))^2 + ((sha1-target_DB10["sha1"]))^2 + ((sha2-target_DB10["sha2"]))^2)^2
  ) %>% 
  arrange(distance) %>% 
  dplyr::slice(1:50)

temp_DB6 <- df_sim_DB6 %>% 
  filter(number_of_1st_lineage==2) %>%
  separate(contribution_of_1st_lineage, c("lineage1", "lineage2")) %>% 
  separate(branch_length_of_1st_lineage, c("shared1", "shared2")) %>% 
  mutate(lin1=ifelse(lineage1>=lineage2, lineage1, lineage2), lin2=ifelse(lineage1>=lineage2, lineage2, lineage1)) %>% 
  mutate(sha1=ifelse(lineage1>=lineage2, shared1,  shared2),  sha2=ifelse(lineage1>=lineage2, shared2,  shared1)) %>%
  mutate(lin1=as.numeric(lin1), lin2=as.numeric(lin2), sha1=as.numeric(sha1), sha2=as.numeric(sha2)) %>% 
  mutate(distance=((number_of_shared_mutations - target_DB6["number_of_shared_mutations"]))^2 +
           #((number_of_1st_lineage      - target_DB6["number_of_1st_lineage"]))^2 + 
           ((multifurcation_score       - target_DB6["multifurcation_score"]))^2 + 
           sqrt(((lin1-target_DB6["lin1"]))^2 + ((lin2-target_DB6["lin2"]))^2 + ((sha1-target_DB6["sha1"]))^2 + ((sha2-target_DB6["sha2"]))^2)^2
  ) %>% 
  arrange(distance) %>%
  dplyr::slice(1:50)

temp_DB9 <- df_sim_DB9 %>% 
  filter(number_of_1st_lineage==2) %>% #filter(at==3 & n==3) %>% 
  separate(contribution_of_1st_lineage, c("lineage1", "lineage2")) %>% 
  separate(branch_length_of_1st_lineage, c("shared1", "shared2")) %>% 
  mutate(lin1=ifelse(lineage1>=lineage2, lineage1, lineage2), lin2=ifelse(lineage1>=lineage2, lineage2, lineage1)) %>% 
  mutate(sha1=ifelse(lineage1>=lineage2, shared1,  shared2),  sha2=ifelse(lineage1>=lineage2, shared2,  shared1)) %>%
  mutate(lin1=as.numeric(lin1), lin2=as.numeric(lin2), sha1=as.numeric(sha1), sha2=as.numeric(sha2)) %>% 
  mutate(distance=((number_of_shared_mutations - target_DB9["number_of_shared_mutations"]))^2 + 
           ((multifurcation_score       - target_DB9["multifurcation_score"]))^2 + 
           sqrt(((lin1-target_DB9["lin1"]))^2 + ((lin2-target_DB9["lin2"]))^2 + ((sha1-target_DB9["sha1"]))^2 + ((sha2-target_DB9["sha2"]))^2)^2
  ) %>% 
  arrange(distance) %>% 
  dplyr::slice(1:40) 

temp_DB8 <- df_sim_DB8 %>% #filter(at==3 & n==3) %>%
  filter(number_of_1st_lineage==2) %>%
  separate(contribution_of_1st_lineage, c("lineage1", "lineage2")) %>% 
  separate(branch_length_of_1st_lineage, c("shared1", "shared2")) %>% 
  mutate(lin1=ifelse(lineage1>=lineage2, lineage1, lineage2), lin2=ifelse(lineage1>=lineage2, lineage2, lineage1)) %>% 
  mutate(sha1=ifelse(lineage1>=lineage2, shared1,  shared2),  sha2=ifelse(lineage1>=lineage2, shared2,  shared1)) %>%
  mutate(lin1=as.numeric(lin1), lin2=as.numeric(lin2), sha1=as.numeric(sha1), sha2=as.numeric(sha2)) %>% 
  mutate(distance=((number_of_shared_mutations - target_DB8["number_of_shared_mutations"]))^2 + 
           #((number_of_multifurcation   - target_DB8["number_of_multifurcation"]))^2 + 
           #(number_of_1st_lineage      - target_DB8["number_of_1st_lineage"])^2 + 
           #(distance_1st_lineage       - target_DB8["distance_1st_lineage"])^2 +
           ((multifurcation_score       - target_DB8["multifurcation_score"]))^2 + 
           #(number_of_shared_mutations_in_1st_lineage - target_DB8["number_of_shared_mutations_in_1st_lineage"])^2 +
           #(minmax_ratio_1st_shared_mutations         - target_DB8["minmax_ratio_1st_shared_mutations"])^2
           #((number_of_1st_lineage -2 ))^2 + 
           sqrt(((lin1-target_DB8["lin1"]))^2 + ((lin2-target_DB8["lin2"]))^2 + ((sha1-target_DB8["sha1"]))^2 + ((sha2-target_DB8["sha2"]))^2)^2
  ) %>% 
  arrange(distance) %>% 
  dplyr::slice(1:40)

df_temp <- bind_rows(temp_DB3 %>% mutate(ID="DB3"),
                     temp_DB6 %>% mutate(ID="DB6"),
                     temp_DB8 %>% mutate(ID="DB8"),
                     temp_DB9 %>% mutate(ID="DB9"),
                     temp_DB10 %>% mutate(ID="DB10")
                     )

df_temp %>% 
  mutate(ID=factor(ID, levels=c("DB3", "DB6", "DB8", "DB9", "DB10"))) %>%
  ggplot() + 
  #geom_density_2d_filled(mapping=aes(x=at,y=n, color=ID), contour=T, alpha=0.2, size=0.5) +
  stat_density_2d_filled(geom="polygon", mapping=aes(x=at,y=n, color=ID, fill=after_stat(level)), contour=T, alpha=0.2, size=0.5) +
  scale_x_continuous(expand=c(0,0), limits = c(1,6), labels = c("2 cells", "4 cells", "8 cells", "16 cells", "32 cells", "64 cells")) + 
  scale_y_continuous(expand=c(0,0), limits = c(1,10), breaks=seq(1,10)) +
  scale_color_jama() + 
  #scale_fill_distiller("Blues") + 
  theme_minimal() + 
  xlab("Stage (s)") + 
  ylab("Number of ICM precursor cells (n)") + 
  guides(fill=FALSE) + 
  scale_color_manual(values=c("DB3"=rgb(40/255,142/255,203/255),
                              "DB6"=rgb(104/255,161/255,133/255),
                              "DB8"=rgb(87/255,80/255,135/255),
                              "DB9"=rgb(108/255,102/255,88/255),
                              "DB10"=rgb(42/255,61/255,68/255))) +
  theme(text=element_text(size=15), legend.title=element_blank(), panel.border = element_rect(color="black", fill=NA),
        legend.position="none",
        axis.ticks=element_line(color="black"),
        axis.text.x=element_text(hjust=-.1),
        plot.margin = unit(c(.5,2,.5,.5), "cm")
        )



