# abc for epiblast selection

df_sim_DB3 <- getSimResults("/home/users/chrono0707/projects/02_DB/01_DB3/14_Simulations/06_epi_hypoblast")
df_sim_DB8 <- getSimResults("/home/users/chrono0707/projects/02_DB/05_DB8/14_simulation/08_epi_hypo/")
df_sim_DB9 <- getSimResults("/home/users/chrono0707/projects/02_DB/06_DB9/14_Simulation/03_epihypo/")

target_DB3 <- c("number_of_shared_mutations"=49, 
                "multifurcation_score"=50/18,
                "lin1" = 22,
                "lin2" = 11,
                "sha1" = 7,
                "sha2" = 7,
                "average_min_max_ratio"=2.792157)

target_DB8 <- c("number_of_shared_mutations"=80, 
                "multifurcation_score"=74/31,
                "lin1" = 26,
                "lin2" = 18,
                "sha1" = 1,
                "sha2" = 4,
                "average_min_max_ratio"=2.862455)

target_DB9 <- c("number_of_shared_mutations"=85, 
                "multifurcation_score"=49/21,
                "lin1" = 21,
                "lin2" = 8,
                "sha1" = 3,
                "sha2" = 5,
                "average_min_max_ratio"=2.80754)

# determine at & n 
df_sim_DB3 %>% 
  filter(number_of_1st_lineage==2) %>% filter(at==3 & n==3) %>% 
  separate(contribution_of_1st_lineage, c("lineage1", "lineage2")) %>% 
  separate(branch_length_of_1st_lineage, c("shared1", "shared2")) %>% 
  mutate(lin1=ifelse(lineage1>=lineage2, lineage1, lineage2), lin2=ifelse(lineage1>=lineage2, lineage2, lineage1)) %>% 
  mutate(sha1=ifelse(lineage1>=lineage2, shared1,  shared2),  sha2=ifelse(lineage1>=lineage2, shared2,  shared1)) %>%
  mutate(lin1=as.numeric(lin1), lin2=as.numeric(lin2), sha1=as.numeric(sha1), sha2=as.numeric(sha2)) %>% 
  mutate(distance=
           (number_of_shared_mutations - target_DB3["number_of_shared_mutations"])^2 + 
           (multifurcation_score       - target_DB3["multifurcation_score"])^2 + 
           (average_min_max_ratio      - target_DB3["average_min_max_ratio"])^2 + 
           #sqrt(
           (lin1-target_DB3["lin1"])^2 + 
           (lin2-target_DB3["lin2"])^2 + 
           (sha1-target_DB3["sha1"])^2 + 
           (sha2-target_DB3["sha2"])^2#)
  ) %>% 
  arrange(distance) %>% 
  dplyr::slice(1:30) %>% 
  ggplot(aes(x=at_epihypo, y=fraction_epiblast)) + 
  stat_density2d(aes(fill=..density..), geom="raster", contour = FALSE) + scale_fill_distiller(palette="Spectral", direction=-1) +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  #theme_test() + 
  xlab("Cell stage") + 
  ylab("Fraciton of epiblast") + 
  theme(legend.pos="none", text=element_text(size=15), panel.border=element_rect(color="black", fill=NA)) 


df_sim_DB8 %>% 
  filter(number_of_1st_lineage==2) %>% #filter(at==3 & n==3) %>% 
  separate(contribution_of_1st_lineage, c("lineage1", "lineage2")) %>% 
  separate(branch_length_of_1st_lineage, c("shared1", "shared2")) %>% 
  mutate(lin1=ifelse(lineage1>=lineage2, lineage1, lineage2), lin2=ifelse(lineage1>=lineage2, lineage2, lineage1)) %>% 
  mutate(sha1=ifelse(lineage1>=lineage2, shared1,  shared2),  sha2=ifelse(lineage1>=lineage2, shared2,  shared1)) %>%
  mutate(lin1=as.numeric(lin1), lin2=as.numeric(lin2), sha1=as.numeric(sha1), sha2=as.numeric(sha2)) %>% 
  mutate(distance=
           (number_of_shared_mutations - target_DB8["number_of_shared_mutations"])^2 + 
           (multifurcation_score       - target_DB8["multifurcation_score"])^2 + 
           (average_min_max_ratio      - target_DB8["average_min_max_ratio"])^2 + 
           #sqrt(
           (lin1-target_DB8["lin1"])^2 + 
           (lin2-target_DB8["lin2"])^2 + 
           (sha1-target_DB8["sha1"])^2 + 
           (sha2-target_DB8["sha2"])^2#)
  ) %>% 
  arrange(distance) %>% 
  dplyr::slice(1:30) %>% 
  ggplot(aes(x=at_epihypo, y=fraction_epiblast)) + 
  stat_density2d(aes(fill=..density..), geom="raster", contour = FALSE) + scale_fill_distiller(palette="Spectral", direction=-1) +
  #scale_x_continuous(expand=c(0,0), limits = c(6,12)) + 
  #scale_y_continuous(expand=c(0,0), limits = c(1,10)) +
  theme_test() + 
  xlab("cleavage") + 
  ylab("Number of ICM precursor cells") + 
  theme(legend.pos="none", text=element_text(size=15))


df_sim_DB9 %>% 
  filter(number_of_1st_lineage==2) %>% filter(at==3 & n==3) %>% 
  separate(contribution_of_1st_lineage, c("lineage1", "lineage2")) %>% 
  separate(branch_length_of_1st_lineage, c("shared1", "shared2")) %>% 
  mutate(lin1=ifelse(lineage1>=lineage2, lineage1, lineage2), lin2=ifelse(lineage1>=lineage2, lineage2, lineage1)) %>% 
  mutate(sha1=ifelse(lineage1>=lineage2, shared1,  shared2),  sha2=ifelse(lineage1>=lineage2, shared2,  shared1)) %>%
  mutate(lin1=as.numeric(lin1), lin2=as.numeric(lin2), sha1=as.numeric(sha1), sha2=as.numeric(sha2)) %>% 
  mutate(distance=
           (number_of_shared_mutations - target_DB9["number_of_shared_mutations"])^2 + 
           (multifurcation_score       - target_DB9["multifurcation_score"])^2 + 
           (average_min_max_ratio      - target_DB9["average_min_max_ratio"])^2 + 
           #sqrt(
           (lin1-target_DB9["lin1"])^2 + 
           (lin2-target_DB9["lin2"])^2 + 
           (sha1-target_DB9["sha1"])^2 + 
           (sha2-target_DB9["sha2"])^2#)
  ) %>% 
  arrange(distance) %>% 
  dplyr::slice(1:30) %>% 
  ggplot(aes(x=at_epihypo, y=fraction_epiblast)) + 
  stat_density2d(aes(fill=..density..), geom="raster", contour = FALSE) + scale_fill_distiller(palette="Spectral", direction=-1) +
  #scale_x_continuous(expand=c(0,0), limits = c(6,12)) + 
  #scale_y_continuous(expand=c(0,0), limits = c(1,10)) +
  theme_test() + 
  xlab("cleavage") + 
  ylab("Number of ICM precursor cells") + 
  theme(legend.pos="none", text=element_text(size=15)) 


# goodness of fit

res <- gfit(target=target_DB3,
            sumstat=df_sim_DB3 %>% 
              filter(number_of_1st_lineage==2) %>% 
              separate(contribution_of_1st_lineage, c("lineage1", "lineage2")) %>% 
              separate(branch_length_of_1st_lineage, c("shared1", "shared2")) %>% 
              mutate(lin1=ifelse(lineage1>=lineage2, lineage1, lineage2), lin2=ifelse(lineage1>=lineage2, lineage2, lineage1)) %>% 
              mutate(sha1=ifelse(lineage1>=lineage2, shared1,  shared2),  sha2=ifelse(lineage1>=lineage2, shared2,  shared1)) %>%
              mutate(lin1=as.numeric(lin1), lin2=as.numeric(lin2), sha1=as.numeric(sha1), sha2=as.numeric(sha2)) %>%
              dplyr::select(names(target_DB3)), nb.replicate=100)
summary(res)
plot(res)



res <- gfit(target=target_DB8,
            sumstat=df_sim_DB8 %>% 
              filter(number_of_1st_lineage==2) %>% 
              separate(contribution_of_1st_lineage, c("lineage1", "lineage2")) %>% 
              separate(branch_length_of_1st_lineage, c("shared1", "shared2")) %>% 
              mutate(lin1=ifelse(lineage1>=lineage2, lineage1, lineage2), lin2=ifelse(lineage1>=lineage2, lineage2, lineage1)) %>% 
              mutate(sha1=ifelse(lineage1>=lineage2, shared1,  shared2),  sha2=ifelse(lineage1>=lineage2, shared2,  shared1)) %>%
              mutate(lin1=as.numeric(lin1), lin2=as.numeric(lin2), sha1=as.numeric(sha1), sha2=as.numeric(sha2)) %>%
              dplyr::select(names(target_DB8)), nb.replicate=100)
summary(res)
plot(res)




# cross validation

cv.res.rej <- cv4abc(data.frame(at_epihypo=df_sim_DB3 %>% filter(number_of_1st_lineage==2) %>% filter(at==3 & n ==3) %>%
                                  separate(contribution_of_1st_lineage, c("lineage1", "lineage2")) %>% 
                                  separate(branch_length_of_1st_lineage, c("shared1", "shared2")) %>% 
                                  mutate(lin1=ifelse(lineage1>=lineage2, lineage1, lineage2), lin2=ifelse(lineage1>=lineage2, lineage2, lineage1)) %>% 
                                  mutate(sha1=ifelse(lineage1>=lineage2, shared1,  shared2),  sha2=ifelse(lineage1>=lineage2, shared2,  shared1)) %>%
                                  mutate(lin1=as.numeric(lin1), lin2=as.numeric(lin2), sha1=as.numeric(sha1), sha2=as.numeric(sha2)) %>% .$at_epihypo,
                                fraction_epiblast=df_sim_DB3 %>% filter(number_of_1st_lineage==2) %>% filter(at==3 & n ==3) %>%
                                  separate(contribution_of_1st_lineage, c("lineage1", "lineage2")) %>% 
                                  separate(branch_length_of_1st_lineage, c("shared1", "shared2")) %>% 
                                  mutate(lin1=ifelse(lineage1>=lineage2, lineage1, lineage2), lin2=ifelse(lineage1>=lineage2, lineage2, lineage1)) %>% 
                                  mutate(sha1=ifelse(lineage1>=lineage2, shared1,  shared2),  sha2=ifelse(lineage1>=lineage2, shared2,  shared1)) %>%
                                  mutate(lin1=as.numeric(lin1), lin2=as.numeric(lin2), sha1=as.numeric(sha1), sha2=as.numeric(sha2)) %>% .$fraction_epiblast),
                     df_sim_DB3 %>% filter(number_of_1st_lineage==2) %>% filter(at==3 & n ==3) %>%
                       separate(contribution_of_1st_lineage, c("lineage1", "lineage2")) %>% 
                       separate(branch_length_of_1st_lineage, c("shared1", "shared2")) %>% 
                       mutate(lin1=ifelse(lineage1>=lineage2, lineage1, lineage2), lin2=ifelse(lineage1>=lineage2, lineage2, lineage1)) %>% 
                       mutate(sha1=ifelse(lineage1>=lineage2, shared1,  shared2),  sha2=ifelse(lineage1>=lineage2, shared2,  shared1)) %>%
                       mutate(lin1=as.numeric(lin1), lin2=as.numeric(lin2), sha1=as.numeric(sha1), sha2=as.numeric(sha2)) %>% 
                       dplyr::select(names(target_DB9)) %>% as.data.frame(),
                     nval=50, tols=c(0.05), statistic="mean", method="rejection")
summary(cv.res.rej)
plot(cv.res.rej)


cv.res.rej <- cv4abc(data.frame(at_epihypo=df_sim_DB8 %>% filter(number_of_1st_lineage==2) %>% 
                                  separate(contribution_of_1st_lineage, c("lineage1", "lineage2")) %>% 
                                  separate(branch_length_of_1st_lineage, c("shared1", "shared2")) %>% 
                                  mutate(lin1=ifelse(lineage1>=lineage2, lineage1, lineage2), lin2=ifelse(lineage1>=lineage2, lineage2, lineage1)) %>% 
                                  mutate(sha1=ifelse(lineage1>=lineage2, shared1,  shared2),  sha2=ifelse(lineage1>=lineage2, shared2,  shared1)) %>%
                                  mutate(lin1=as.numeric(lin1), lin2=as.numeric(lin2), sha1=as.numeric(sha1), sha2=as.numeric(sha2)) %>% .$at_epihypo,
                                fraction_epiblast=df_sim_DB8 %>% filter(number_of_1st_lineage==2) %>% 
                                  separate(contribution_of_1st_lineage, c("lineage1", "lineage2")) %>% 
                                  separate(branch_length_of_1st_lineage, c("shared1", "shared2")) %>% 
                                  mutate(lin1=ifelse(lineage1>=lineage2, lineage1, lineage2), lin2=ifelse(lineage1>=lineage2, lineage2, lineage1)) %>% 
                                  mutate(sha1=ifelse(lineage1>=lineage2, shared1,  shared2),  sha2=ifelse(lineage1>=lineage2, shared2,  shared1)) %>%
                                  mutate(lin1=as.numeric(lin1), lin2=as.numeric(lin2), sha1=as.numeric(sha1), sha2=as.numeric(sha2)) %>% .$fraction_epiblast,
                                at=df_sim_DB8 %>% filter(number_of_1st_lineage==2) %>% 
                                  separate(contribution_of_1st_lineage, c("lineage1", "lineage2")) %>% 
                                  separate(branch_length_of_1st_lineage, c("shared1", "shared2")) %>% 
                                  mutate(lin1=ifelse(lineage1>=lineage2, lineage1, lineage2), lin2=ifelse(lineage1>=lineage2, lineage2, lineage1)) %>% 
                                  mutate(sha1=ifelse(lineage1>=lineage2, shared1,  shared2),  sha2=ifelse(lineage1>=lineage2, shared2,  shared1)) %>%
                                  mutate(lin1=as.numeric(lin1), lin2=as.numeric(lin2), sha1=as.numeric(sha1), sha2=as.numeric(sha2)) %>% .$at,
                                n=df_sim_DB8 %>% filter(number_of_1st_lineage==2) %>% 
                                  separate(contribution_of_1st_lineage, c("lineage1", "lineage2")) %>% 
                                  separate(branch_length_of_1st_lineage, c("shared1", "shared2")) %>% 
                                  mutate(lin1=ifelse(lineage1>=lineage2, lineage1, lineage2), lin2=ifelse(lineage1>=lineage2, lineage2, lineage1)) %>% 
                                  mutate(sha1=ifelse(lineage1>=lineage2, shared1,  shared2),  sha2=ifelse(lineage1>=lineage2, shared2,  shared1)) %>%
                                  mutate(lin1=as.numeric(lin1), lin2=as.numeric(lin2), sha1=as.numeric(sha1), sha2=as.numeric(sha2)) %>% .$n,
                                rate=df_sim_DB8 %>% filter(number_of_1st_lineage==2) %>% 
                                  separate(contribution_of_1st_lineage, c("lineage1", "lineage2")) %>% 
                                  separate(branch_length_of_1st_lineage, c("shared1", "shared2")) %>% 
                                  mutate(lin1=ifelse(lineage1>=lineage2, lineage1, lineage2), lin2=ifelse(lineage1>=lineage2, lineage2, lineage1)) %>% 
                                  mutate(sha1=ifelse(lineage1>=lineage2, shared1,  shared2),  sha2=ifelse(lineage1>=lineage2, shared2,  shared1)) %>%
                                  mutate(lin1=as.numeric(lin1), lin2=as.numeric(lin2), sha1=as.numeric(sha1), sha2=as.numeric(sha2)) %>% .$rate
                                ),
                     
                     df_sim_DB8 %>% filter(number_of_1st_lineage==2) %>% 
                       separate(contribution_of_1st_lineage, c("lineage1", "lineage2")) %>% 
                       separate(branch_length_of_1st_lineage, c("shared1", "shared2")) %>% 
                       mutate(lin1=ifelse(lineage1>=lineage2, lineage1, lineage2), lin2=ifelse(lineage1>=lineage2, lineage2, lineage1)) %>% 
                       mutate(sha1=ifelse(lineage1>=lineage2, shared1,  shared2),  sha2=ifelse(lineage1>=lineage2, shared2,  shared1)) %>%
                       mutate(lin1=as.numeric(lin1), lin2=as.numeric(lin2), sha1=as.numeric(sha1), sha2=as.numeric(sha2)) %>% 
                       dplyr::select(names(target_DB8)) %>% as.data.frame(),
                     nval=50, tols=c(0.05), statistic="mean", method="rejection")
summary(cv.res.rej)
plot(cv.res.rej)




# parameter estimation
res <- abc(target  = target_DB3,
           param   = data.frame(at_epihypo=df_sim_DB3 %>% filter(at==3 & n==3) %>% 
                                  filter(number_of_1st_lineage==2) %>% 
                                  .$at_epihypo,
                                fraction_epiblast=df_sim_DB3 %>% filter(at==3 & n==3) %>% 
                                  filter(number_of_1st_lineage==2) %>% 
                                  .$fraction_epiblast
           ),
           sumstat = df_sim_DB3 %>% filter(number_of_1st_lineage==2) %>% filter(at==3 & n==3) %>% 
             separate(contribution_of_1st_lineage, c("lineage1", "lineage2")) %>% 
             separate(branch_length_of_1st_lineage, c("shared1", "shared2")) %>% 
             mutate(lin1=ifelse(lineage1>=lineage2, lineage1, lineage2), lin2=ifelse(lineage1>=lineage2, lineage2, lineage1)) %>% 
             mutate(sha1=ifelse(lineage1>=lineage2, shared1,  shared2),  sha2=ifelse(lineage1>=lineage2, shared2,  shared1)) %>%
             mutate(lin1=as.numeric(lin1), lin2=as.numeric(lin2), sha1=as.numeric(sha1), sha2=as.numeric(sha2)) %>% 
             dplyr::select(names(target_DB3)) %>% as.data.frame(),
           tol = 0.05, 
           transf=c("none"), 
           numnet=50,
           sizenet=10,
           method="neuralnet"
)
summary(res)
hist(res)




