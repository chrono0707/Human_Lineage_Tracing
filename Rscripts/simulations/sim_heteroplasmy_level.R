library(tidyverse)
library(parallel)
library(magrittr)
# Embryogenesis simulation 
# ====================================================================================================================================

# structure containing descriptions of reactions
Cell <- setClass("Cell",
                 slots = c(
                   barcode       = "character",  # to mark the cellular origin
                   mtDNA         = "data.frame" # mitochondrial DNA 
                 ))

# simulation per generation (input: mother cell, output: two daughter cells)
Cell_division <- function(mother) { 
  
  divide_at <- sample(1:nrow(mother@mtDNA), 1) # define the cleavage furrow
  
  # clustering two mtDNA groups that will be transferred to each daughter cells 
  if(divide_at <= nrow(mother@mtDNA)/2) {
    daughter_mtDNA_1 <- bind_rows(mother@mtDNA %>% head(divide_at),
                                  mother@mtDNA %>% tail(nrow(mother@mtDNA)/2 - divide_at))
    daughter_mtDNA_2 <- mother@mtDNA %>% dplyr::slice((divide_at+1):(nrow(mother@mtDNA) - (nrow(mother@mtDNA)/2 - divide_at - 1) - 1))
  } else {
    daughter_mtDNA_1 <- bind_rows(mother@mtDNA %>% tail(nrow(mother@mtDNA) - divide_at),
                                  mother@mtDNA %>% head(nrow(mother@mtDNA)/2 - (nrow(mother@mtDNA) - divide_at)))
    daughter_mtDNA_2 <- mother@mtDNA %>% dplyr::slice((nrow(mother@mtDNA)/2 - (nrow(mother@mtDNA) - divide_at) + 1):(divide_at))
  }
  
  daughters <- list(   
    Cell(barcode = paste0(mother@barcode, "1"), 
         mtDNA   = bind_rows(daughter_mtDNA_1, daughter_mtDNA_1) %>% arrange(clone)),
    Cell(barcode = paste0(mother@barcode, "2"),
         mtDNA   = bind_rows(daughter_mtDNA_2, daughter_mtDNA_2) %>% arrange(clone))
  )
  
  invisible(daughters)
}



# simulate embryogenesis 
Embryogenesis <- function(mother,               # egg
                          generations = 10,      # how many cell generataions will be simulated?
                          mc.cores = 2) {       # how many cores will be employed?
  
  # create buffer and register MRCA
  cells         <- vector(mode="list", length=generations + 1) 
  names(cells)  <- str_c("G", c(0, seq_len(generations)))
  
  # register the fertilized egg
  cells[["G0"]] <- list(mother) 
  
  # cell division
  for(g in 2:(generations + 1)) {
    print(sprintf("generation: %d", g))
    cells[[g]] <- mclapply(cells[[g-1]], Cell_division, mc.cores = mc.cores) %>% unlist() # parallelize the generation of daughter cells
    #cells[[g]] <- lapply(cells[[g-1]], Cell_division) %>%  unlist() # parallelize the generation of daughter cells
  }
  
  invisible(cells)
}

Simulate_sequencing <- function(cells,                # output of Embryogenesis function 
                                sampling,             # how many samples will be taken from the cells?
                                mc.cores = 2          # How many cores will be employed?
                                ) { 
  
  terminal_cells <- cells[[length(cells)]] # cells from terminal stage 
  terminal_cells <- terminal_cells[sort(sample(seq_len(length(terminal_cells)), sampling, replace = FALSE))] # sampling 
  
  # convert list to dataframe 
  df <- do.call(rbind, mclapply(terminal_cells, function(x) {
    return(tibble(cell=x@barcode, VAF=sum(x@mtDNA$clone)/nrow(x@mtDNA))) 
  }, mc.cores = mc.cores))
  
  df_pos <- df %>% filter(VAF!=0)
  
  return(list(n_pos_cells          = nrow(df_pos),
              proportion_pos_cells = nrow(df_pos)/nrow(df),
              median_VAF           = median(df_pos$VAF),
              mean_VAF             = mean(df_pos$VAF),
              sd_VAF               = sd(df_pos$VAF),
              max_VAF              = max(df_pos$VAF),
              min_VAF              = min(df_pos$VAF)))
}


args             <- commandArgs(trailingOnly=TRUE)
iter             <- as.integer(args[1])
total_mtDNA      <- as.integer(args[2])
min_initial_fraction <- as.numeric(args[3])
max_initial_fraction <- as.numeric(args[4])
generations      <- as.integer(args[5])
sampling         <- as.integer(args[6])
out              <- args[7]

cat("initial_fraction\tsampling\ttotal_mtDNA\tn_pos_cells\tproportion_pos_cells\tmedian_VAF\tmean_VAF\tsd_VAF\tmax_VAF\tmin_VAF\n", file=out, append=TRUE)

for(i in 1:iter) {

  initial_fraction <- runif(1, min=min_initial_fraction, max=max_initial_fraction)

  print(sprintf("initial fraction: %.3f, iter: %d", initial_fraction, i))

  cells <- Embryogenesis(Cell(barcode  = "0",
                              mtDNA    = tibble(clone=rep(c(1,0), times=c(as.integer(total_mtDNA * initial_fraction), 
                                                                          total_mtDNA - as.integer(total_mtDNA * initial_fraction))))),
                         generations = generations)

  #result <- Simulate_sequencing(cells, length(cells[[length(cells)]]))
  result <- Simulate_sequencing(cells, sampling)

  cat(sprintf("%f\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n",
    initial_fraction,
    sampling,
    total_mtDNA,
    result$n_pos_cells,
    result$proportion_pos_cells,
    result$median_VAF,
    result$median_VAF,
    result$sd_VAF,
    result$max_VAF,
    result$median_VAF), file=out, append=TRUE)

}