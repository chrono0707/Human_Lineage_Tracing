library(tidyverse)
library(magittr)
library(abc)

# heteroplasy in oocyte simulation

# load data
df_temp <- getSimResults("/home/users/chrono0707/projects/02_DB/07_DB10/14_Simulation/10_heteroplasmy_in_oocyte/01_1000copies")
df_temp <- getSimResults("/home/users/chrono0707/projects/02_DB/01_DB3/14_Simulations/10_heteroplasmy_in_oocyte/")
df_temp <- getSimResults("/home/users/chrono0707/projects/02_DB/02_DB6/14_Simulation/10_heteroplasmy_in_oocyte/")
df_temp <- getSimResults("/home/users/chrono0707/projects/02_DB/05_DB8/14_simulation/10_heteroplasmy_in_oocyte/")
df_temp <- getSimResults("/home/users/chrono0707/projects/02_DB/06_DB9/14_Simulation/10_heteroplamsy_in_oocyte/")

# define target
target <- c("n_pos_cells"=14, "mean_VAF"=0.8742857) # DB10 16256 C>T
target <- c("n_pos_cells"=1,  "mean_VAF"=0.99)      # homoplasmic singleton

target <- c("n_pos_cells"=2,  "mean_VAF"=1.0)       # DB3 15649 A>G
target <- c("n_pos_cells"=2,  "mean_VAF"=0.025)     # DB3 2876  G>A
target <- c("n_pos_cells"=2,  "mean_VAF"=0.020)     # DB3 9554  G>A

target <- c("n_pos_cells"=2,  "mean_VAF"=0.52)      # DB6 11642 G>A
target <- c("n_pos_cells"=2,  "mean_VAF"=0.63)      # DB6 10161 A>G
target <- c("n_pos_cells"=2,  "mean_VAF"=0.205)     # DB6 7075  G>A
target <- c("n_pos_cells"=2,  "mean_VAF"=0.565)     # DB6 1990  G>A
target <- c("n_pos_cells"=2,  "mean_VAF"=0.055)     # DB6 12889 G>A

target <- c("n_pos_cells"=4,  "mean_VAF"=0.3924654) # DB8 5894 A>G
target <- c("n_pos_cells"=2,  "mean_VAF"=0.3818774) # DB8 16172 T>C
target <- c("n_pos_cells"=2,  "mean_VAF"=0.3358135) # DB8 13681 A>G

target <- c("n_pos_cells"=2,  "mean_VAF"=0.8915) # DB9 12359 C>T
target <- c("n_pos_cells"=2,  "mean_VAF"=0.13598) # DB9 9478 T>G
target <- c("n_pos_cells"=2,  "mean_VAF"=0.141629) # DB9 2573 G>A



# check goodness of fit
res <- gfit(target=target,
            sumstat=df_temp %>% dplyr::select(names(target)), nb.replicate=500, tol=0.10)
summary(res)
plot(res)

res$dist.sim %>% enframe() %>% ggplot(aes(x=value)) + 
  geom_histogram(na.rm=T, fill=rgb(32/255,93/255,181/255)) + 
  geom_vline(xintercept = res$dist.obs, color="red") + 
  theme_minimal() + 
  xlab("Distance") + ylab("Count") + 
  theme(text=element_text(size=15))

# cross validation
cv.res.rej <- cv4abc(param=data.frame(initial_fraction=df_temp$initial_fraction),
                     sumstat=df_temp %>% dplyr::select(names(target)) %>% as.data.frame(),
                     nval=50, tols=c(0.05), statistic="mean", method="neuralnet")
summary(cv.res.rej)
plot(cv.res.rej)

as.tibble(cv.res.rej$true) %>% mutate(estim=cv.res.rej$estim$tol0.05) %>% 
  ggplot(aes(x=initial_fraction, y=estim)) +
  geom_point(size=3) + 
  annotate("segment", x=0, y=0, xend=1, yend=1, color="red", size=0.3) + 
  theme_minimal() + 
  xlab("True") + ylab("Estimated") + 
  theme(text=element_text(size=15))

# parameter inference
res <- abc(target  = target,
           param   = data.frame(initial_fraction=df_temp$initial_fraction),
           sumstat = df_temp %>% dplyr::select(names(target)) %>% as.data.frame(),
           tol = 0.05, 
           transf=c("log"), 
           method="neuralnet"
)

summary(res)


Rtibble(initial_fraction=res$adj.values) %>% ggplot(aes(x=initial_fraction)) + 
  geom_density(fill=rgb(234/255,181/255,10/255)) + 
  xlab("Heteroplamic level of 16256 C>T") + 
  ylab("Density") + 
  theme_minimal() + 
  scale_y_continuous(expand=c(0,0)) + #, limits=c(0,25))+
  theme(text=element_text(size=13), axis.line=element_line(color="black"))




m <- matrix(res$adj.values, ncol=res$numparam)
mean <- apply(m, 2, weighted.mean, w=res$weights)
quants <- apply(m,2, function(x) rq(x~1, tau=c((1-0.95)/2, 0.5, 1-(1-0.95)/2), weights=res$weights)$coef)






# heteroplasy contour

df_temp <- getSimResults("/home/users/chrono0707/projects/02_DB/07_DB10/14_Simulation/10_heteroplasmy_in_oocyte/01_1000copies/01_contour")
