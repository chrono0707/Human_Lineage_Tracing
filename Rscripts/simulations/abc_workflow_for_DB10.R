# ABC simulation work flow

# read sim results
df_sim_DB10 <- getSimResults("/home/users/chrono0707/projects/02_DB/07_DB10/14_Simulation/01_initial_screening/")

# define target values
target_DB10 <- c("number_of_shared_mutations"=60, 
                "number_of_1st_lineage"=2,
                "number_of_multifurcation"=25, 
                "multifurcation_score"=65/25,
                "lin1" = 26,
                "lin2" = 12,
                "sha1" = 4,
                "sha2" = 6
                )

# determine at & n 
df_sim_DB10 %>% 
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
  dplyr::slice(1:40) %>% 
  ggplot(aes(x=at, y=n)) + stat_density2d(aes(fill=..density..), geom="raster", contour = FALSE) + scale_fill_distiller(palette="Spectral", direction=-1) +
  scale_x_continuous(expand=c(0,0), limits = c(1,6)) + 
  scale_y_continuous(expand=c(0,0), limits = c(1,10)) + 
  theme_test() + 
  xlab("cleavage") + 
  ylab("Number of ICM precursor cells") + 
  theme(legend.pos="none", text=element_text(size=15)) 


# read sim results with specific at & n
df_sim_DB10 <- getSimResults("/home/users/chrono0707/projects/02_DB/07_DB10/14_Simulation/02_at_3_n_3/")

# check goodness of fit
res <- gfit(target=target_DB10, #[c("number_of_shared_mutations", "number_of_1st_lineage", "distance_1st_lineage", "multifurcation_score")], 
            sumstat=df_sim_DB10 %>% 
              filter(number_of_1st_lineage==2) %>% 
              separate(contribution_of_1st_lineage, c("lineage1", "lineage2")) %>% 
              separate(branch_length_of_1st_lineage, c("shared1", "shared2")) %>% 
              mutate(lin1=ifelse(lineage1>=lineage2, lineage1, lineage2), lin2=ifelse(lineage1>=lineage2, lineage2, lineage1)) %>% 
              mutate(sha1=ifelse(lineage1>=lineage2, shared1,  shared2),  sha2=ifelse(lineage1>=lineage2, shared2,  shared1)) %>%
              mutate(lin1=as.numeric(lin1), lin2=as.numeric(lin2), sha1=as.numeric(sha1), sha2=as.numeric(sha2)) %>%
              dplyr::select(names(target_DB10)), nb.replicate=100)
summary(res)
plot(res)


# scatter plot 
df_sim_DB10 %>%
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
  dplyr::slice(1:50) %>% 
  ggplot(aes(x=initial_rate, y=rate)) + stat_density2d(aes(fill=..density..), geom="raster", contour = FALSE) + scale_fill_distiller(palette="Spectral", direction=-1) +
  #scale_x_continuous(expand=c(0,0), limits = c(1,5)) + 
  #scale_y_continuous(expand=c(0,0), limits = c(1,5)) + 
  theme_test() + 
  xlab("cleavage") + 
  ylab("Number of ICM precursor cells") + 
  theme(legend.pos="none", text=element_text(size=15)) 

# cross validation
cv.res.rej <- cv4abc(data.frame(rate=df_sim_DB9 %>% filter(number_of_1st_lineage==2) %>% filter(at==3 & n ==3) %>%
                       separate(contribution_of_1st_lineage, c("lineage1", "lineage2")) %>% 
                       separate(branch_length_of_1st_lineage, c("shared1", "shared2")) %>% 
                       mutate(lin1=ifelse(lineage1>=lineage2, lineage1, lineage2), lin2=ifelse(lineage1>=lineage2, lineage2, lineage1)) %>% 
                       mutate(sha1=ifelse(lineage1>=lineage2, shared1,  shared2),  sha2=ifelse(lineage1>=lineage2, shared2,  shared1)) %>%
                       mutate(lin1=as.numeric(lin1), lin2=as.numeric(lin2), sha1=as.numeric(sha1), sha2=as.numeric(sha2)) %>% .$rate),
                     df_sim_DB9 %>% filter(number_of_1st_lineage==2) %>% filter(at==3 & n ==3) %>%
                       separate(contribution_of_1st_lineage, c("lineage1", "lineage2")) %>% 
                       separate(branch_length_of_1st_lineage, c("shared1", "shared2")) %>% 
                       mutate(lin1=ifelse(lineage1>=lineage2, lineage1, lineage2), lin2=ifelse(lineage1>=lineage2, lineage2, lineage1)) %>% 
                       mutate(sha1=ifelse(lineage1>=lineage2, shared1,  shared2),  sha2=ifelse(lineage1>=lineage2, shared2,  shared1)) %>%
                       mutate(lin1=as.numeric(lin1), lin2=as.numeric(lin2), sha1=as.numeric(sha1), sha2=as.numeric(sha2)) %>% 
                       dplyr::select(names(target_DB9)) %>% as.data.frame(),
                     nval=50, tols=c(0.05), statistic="mean", method="rejection")
summary(cv.res.rej)
plot(cv.res.rej)

# parameter inference
res <- abc(target  = target_DB10,
           param   = data.frame(rate=df_sim_DB10 %>% 
                                  filter(number_of_1st_lineage==2) %>% 
                                  .$rate,
                                initia_rate=df_sim_DB10 %>% 
                                  filter(number_of_1st_lineage==2) %>% 
                                  .$initial_rate
                                ),
           sumstat = df_sim_DB10 %>% filter(number_of_1st_lineage==2) %>%
             separate(contribution_of_1st_lineage, c("lineage1", "lineage2")) %>% 
             separate(branch_length_of_1st_lineage, c("shared1", "shared2")) %>% 
             mutate(lin1=ifelse(lineage1>=lineage2, lineage1, lineage2), lin2=ifelse(lineage1>=lineage2, lineage2, lineage1)) %>% 
             mutate(sha1=ifelse(lineage1>=lineage2, shared1,  shared2),  sha2=ifelse(lineage1>=lineage2, shared2,  shared1)) %>%
             mutate(lin1=as.numeric(lin1), lin2=as.numeric(lin2), sha1=as.numeric(sha1), sha2=as.numeric(sha2)) %>% 
             dplyr::select(names(target_DB10)) %>% as.data.frame(),
           tol = 0.05, 
           #transf=c("log"), 
           method="neuralnet"
           )
summary(res)
hist(res)
