# ABC simulation work flow

# read sim results
df_sim_DB6 <- getSimResults("/home/users/chrono0707/projects/02_DB/02_DB6/14_Simulation/03_initial_screening_narrowed_at")


# define target values
target_DB6 <- c("number_of_shared_mutations"=186, 
                #"number_of_1st_lineage"=2,
                "multifurcation_score"=174/64,
                "lin1" = 99,
                "lin2" = 15,
                "sha1" = 1,
                "sha2" = 2
                )

# determine at & n 
df_sim_DB6 %>% 
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
  dplyr::slice(1:50) %>% 
  ggplot(aes(x=at, y=n)) + stat_density2d(aes(fill=..density..), geom="raster", contour = FALSE) + scale_fill_distiller(palette="Spectral", direction=-1) +
  scale_x_continuous(expand=c(0,0), limits=c(1,6)) + 
  scale_y_continuous(expand=c(0,0), limits=c(1,10)) + 
  theme_test() + 
  xlab("cleavage") + 
  ylab("Number of ICM precursor cells") + 
  theme(legend.pos="none", text=element_text(size=15)) 


# read sim results with specific at & n
df_sim_DB6 <- getSimResults("/home/users/chrono0707/projects/02_DB/02_DB6/14_Simulation/04_at_4_n_6/")

# check goodness of fit
res <- gfit(target=target_DB6, #[c("number_of_shared_mutations", "number_of_1st_lineage", "distance_1st_lineage", "multifurcation_score")], 
            sumstat=df_sim_DB6 %>% 
              filter(number_of_1st_lineage==2) %>% 
              separate(contribution_of_1st_lineage, c("lineage1", "lineage2")) %>% 
              separate(branch_length_of_1st_lineage, c("shared1", "shared2")) %>% 
              mutate(lin1=ifelse(lineage1>=lineage2, lineage1, lineage2), lin2=ifelse(lineage1>=lineage2, lineage2, lineage1)) %>% 
              mutate(sha1=ifelse(lineage1>=lineage2, shared1,  shared2),  sha2=ifelse(lineage1>=lineage2, shared2,  shared1)) %>%
              mutate(lin1=as.numeric(lin1), lin2=as.numeric(lin2), sha1=as.numeric(sha1), sha2=as.numeric(sha2)) %>%
              dplyr::select(names(target_DB6)), nb.replicate=100)
summary(res)
plot(res)


# scatter plot 
df_sim_DB8 %>%
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
           sqrt(((sha1-1))^2 + ((sha2-4))^2) # + (number_of_1st_lineage - 2)^2
  ) %>% 
  arrange(distance) %>% 
  filter(number_of_1st_lineage==2) %>% 
  #dplyr::slice(1:100) %>% 
  filter(distance < 0.01) %>% 
  ggplot(aes(x=initial_rate, y=rate, color=distance)) + geom_point(alpha=0.3, size=3)

# cross validation
cv.res.rej <- cv4abc(data.frame(rate=df_sim_DB6 %>% filter(number_of_1st_lineage==2) %>% 
                       separate(contribution_of_1st_lineage, c("lineage1", "lineage2")) %>% 
                       separate(branch_length_of_1st_lineage, c("shared1", "shared2")) %>% 
                       mutate(lin1=ifelse(lineage1>=lineage2, lineage1, lineage2), lin2=ifelse(lineage1>=lineage2, lineage2, lineage1)) %>% 
                       mutate(sha1=ifelse(lineage1>=lineage2, shared1,  shared2),  sha2=ifelse(lineage1>=lineage2, shared2,  shared1)) %>%
                       mutate(lin1=as.numeric(lin1), lin2=as.numeric(lin2), sha1=as.numeric(sha1), sha2=as.numeric(sha2)) %>% .$rate),
                     df_sim_DB6 %>% filter(number_of_1st_lineage==2) %>%
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
res <- abc(target  = target_DB6,
           param   = data.frame(rate=df_sim_DB6 %>% 
                                  filter(number_of_1st_lineage==2) %>% 
                                  .$rate,
                                initia_rate=df_sim_DB6 %>% 
                                  filter(number_of_1st_lineage==2) %>% 
                                  .$initial_rate
           ),
           sumstat = df_sim_DB6 %>% filter(number_of_1st_lineage==2) %>%
             separate(contribution_of_1st_lineage, c("lineage1", "lineage2")) %>% 
             separate(branch_length_of_1st_lineage, c("shared1", "shared2")) %>% 
             mutate(lin1=ifelse(lineage1>=lineage2, lineage1, lineage2), lin2=ifelse(lineage1>=lineage2, lineage2, lineage1)) %>% 
             mutate(sha1=ifelse(lineage1>=lineage2, shared1,  shared2),  sha2=ifelse(lineage1>=lineage2, shared2,  shared1)) %>%
             mutate(lin1=as.numeric(lin1), lin2=as.numeric(lin2), sha1=as.numeric(sha1), sha2=as.numeric(sha2)) %>% 
             dplyr::select(names(target_DB6)) %>% as.data.frame(),
           tol = 0.05, 
           transf=c("log"), 
           method="neuralnet"
)

summary(res)
hist(res)
