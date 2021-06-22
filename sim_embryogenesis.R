library(tidyverse)
library(magrittr)
library(parallel)

# Embryogenesis simulation 
# ====================================================================================================================================

# structure containing descriptions of reactions
Cell <- setClass("Cell",
                 slots = c(
                   barcode       = "character",  # to mark the cellular origin
                   nuclear       = "data.frame", # nuclear genome mutation 
                   mrate_nuclear = "numeric"    # nuclear genome mutation rate (단위: mutations/division)
                   #mrate_1st_div = "numeric"     # mutation rate of 1st division (단위: mutations/division)
                 ))

# simulation per generation (input: mother cell, output: two daughter cells)
Cell_division <- function(mother) { 
  
  # nuclear mutation follows poisson distribution 
  n_nuclear_muts <- rpois(1, mother@mrate_nuclear)
  if(n_nuclear_muts) {
    mother@nuclear %<>% bind_rows(tibble(id=str_c(mother@barcode, ":", 1:n_nuclear_muts))) # register novel mutations 
  }
  
  daughters <- list(   
    Cell(barcode = paste0(mother@barcode, "1"), 
         mrate_nuclear = mother@mrate_nuclear,
         nuclear       = mother@nuclear),
    Cell(barcode = paste0(mother@barcode, "2"),
         mrate_nuclear = mother@mrate_nuclear, 
         nuclear       = mother@nuclear)
  )
  
  invisible(daughters)
}


# simulate embryogenesis 
Embryogenesis <- function(mother,               # egg
                          mrate_1st_div,        # mutation rate of 1st division 
                          generations = 7,      # how many cell generataions will be simulated?
                          transition = 2,       # From which stage ZGA will occur?
                          mc.cores = 2) {       # how many cores will be employed?

  # create buffer and register MRCA
  cells         <- vector(mode="list", length=generations + 1) 
  names(cells)  <- str_c("G", c(0, seq_len(generations)))
  
  # divisions before ZGA
  mrate <- mother@mrate_nuclear
  mother@mrate_nuclear <- mrate_1st_div   # initial division rate 
  cells[["G0"]] <- list(mother) 
  
  for(g in 2:(transition + 1))                       cells[[g]] <- mclapply(cells[[g - 1]], Cell_division, mc.cores = mc.cores) %>% unlist() # division with initial mutation rate 
  for(e in seq_len(length(cells[[transition + 1]]))) cells[[transition + 1]][[e]]@mrate_nuclear <- mrate # restore zygotic mutation rate
  for(g in (transition + 2):(generations + 1))       cells[[g]] <- mclapply(cells[[g-1]], Cell_division, mc.cores = mc.cores) %>% unlist() # parallelize the generation of daughter cells
  
  invisible(cells)
}

Simulate_sequencing <- function(cells,                # output of Embryogenesis function 
                                sampling,             # how many samples will be taken from the cells?
                                at,                   # At which generation (with 2^cells) ICM will be selected? 
                                n,                    # How many embryonic cell will be ICM precursor? 
                                mc.cores = 2,         # How many cores will be employed?
                                draw.tree   = TRUE) { # draw phylogenetic tree  
  
  precursors <- sample(lapply(cells[[at + 1]], function(x) x@barcode), n) # ICM precursor selection
  
  terminal_cells <- cells[[length(cells)]] # cells from terminal stage 
  terminal_cells <- terminal_cells[which(unlist(mclapply(terminal_cells, function(x) str_sub(x@barcode, end=(at + 1)), mc.cores = mc.cores)) %in% precursors)] # select cells originated from ICM precursors
  #terminal_cells <- sample(terminal_cells, sampling, replace = FALSE) # sampling
  terminal_cells <- terminal_cells[sort(sample(seq_len(length(terminal_cells)), sampling, replace = FALSE))] # sampling 
  
  # convert list to dataframe 
  df <- do.call(rbind, mclapply(terminal_cells, function(x) {
    x@nuclear %<>% mutate(barcode = x@barcode)
    return(x@nuclear) 
    }, mc.cores = mc.cores))
  
  df %<>% mutate(bin=1) %>% spread(key=barcode, value=bin, fill=0) # convert to wide-form dataset 
  df %<>% arrange(desc(rowSums(.[, 2:ncol(.)])))                   # sorting
  
  # save number of shared mutations before further processing 
  n_shared_mutations <- df %>% select(-id) %>% filter(rowSums(.) > 1 & rowSums(.) != sampling) %>% nrow()
  
  
  df %<>% unite("barcode", 2:ncol(.), remove = FALSE) %>%               # sharing pattern of each mutation
    group_by(barcode) %>% mutate(count=n(), node=NA) %>% ungroup() %>%  # count: branch length
    filter(!duplicated(barcode)) %>% 
    separate(id, c("id", "id_num")) %>% 
    select(-id_num, -barcode) %>% 
    select(id, count, node, everything())
  
  # to annotate tip 
  identity           <- diag(ncol(df)-3) # create identity matrix 
  colnames(identity) <- colnames(df)[4:ncol(df)]
  identity           <- as_tibble(identity) %>% mutate(id=colnames(df)[4:ncol(df)], count=0, node=NA) %>% select(id, count, node, everything()) 
  df %<>% bind_rows(identity)
  
  # find nodes of each branch
  for(i in 2:nrow(df)) {
    check <- 0
    for(j in (i-1):1) {
      if(sum(df[i, 4:ncol(df)] * df[j, 4:ncol(df)]) == sum(df[i, 4:ncol(df)])) { check <- 1; break }
    }
    if(check) df$node[i] <- df$id[j]
  }
  df %<>% mutate(node=ifelse(is.na(node), "G0", node))
  
  
  # count contribution of lineages emerged during the first 2 two generations
  n_samples_from_1st_branch <- integer()
  
  if(sum(df$node == "G0") > 1) {
    first_lineages <- df$id[df$node=="G0"]  
    shared_mutations_in_1st_lineage <- df$count[df$node=="G0"]
  } else {
    first_lineages <- df$id[df$node %in% df$id[df$node=="G0"]]
    shared_mutations_in_1st_lineage <- df$count[df$node %in% df$id[df$node=="G0"]]
  }
  
  for(l in first_lineages) {
    n_samples_from_1st_branch <- c(n_samples_from_1st_branch, sum(startsWith(colnames(df), l)))
  }
  
  # generate newick format
  df_newick <- df %>% mutate(line=str_c(id, ":", count)) %>% 
    group_by(node) %>% 
    mutate(lines=paste0(line, collapse=","), branches=n(), count=sum(count)) %>%  # branches: branches per node, count: mutation number assigned to each node
    ungroup() %>% 
    select(node, branches, count, lines) %>% 
    filter(!duplicated(node)) 
  
  for(i in 1:nrow(df_newick)) {
    df_newick %<>% mutate(lines=str_replace(lines, eval(parse(text=sprintf("'%s:'", df_newick$node[i]))), str_c("(", df_newick$lines[i], ")", df_newick$node[i], ":")))
  }
  
  # save results in list   
  result <- list(tree                                      = str_c("(", df_newick$lines[df_newick$node == "G0"], ")G0:0;"), # tree structure in newick format 
                 multifurcation_score                      = sum(df_newick$branches[df_newick$branches > 1]) / nrow(df_newick[df_newick$branches > 1,]), # branch가 여러개인 node에 대해서 multifurcation score 계산 (G0 node 여도 branch가 하나이면 계산에서 제외됨)
                 number_of_multifurcation                  = nrow(df_newick[df_newick$branches  > 1, ]), # multifurcation (include bifurcation)
                 number_of_bifurcation                     = nrow(df_newick[df_newick$branches == 2, ]), # bifurcation 
                 from_first_cell                           = sum(startsWith(colnames(df), "01")), # cells from one of the first two cells 
                 from_second_cell                          = sum(startsWith(colnames(df), "02")), # cells from one of the first two cells 
                 number_of_1st_lineage                     = length(n_samples_from_1st_branch), # degree of multifurcation at the very first lineage formation 
                 contribution_of_1st_lineage               = paste0(n_samples_from_1st_branch, collapse=";"), # number of samples contributed from each lineages
                 distance_1st_lineage                      = sum((n_samples_from_1st_branch/sampling - 1/length(n_samples_from_1st_branch))^2)/length(n_samples_from_1st_branch), # distance from uniform distribution
                 number_of_shared_mutations                = n_shared_mutations, # shared mutations (exclude mutations in G0 branch if G0 branch has only one lineage)
                 branch_length_of_1st_lineage              = paste0(shared_mutations_in_1st_lineage, collapse=";"), 
                 number_of_shared_mutations_in_1st_lineage = sum(shared_mutations_in_1st_lineage), # shared mutations in 1st lineage
                 minmax_ratio_1st_shared_mutations         = min(shared_mutations_in_1st_lineage)/max(shared_mutations_in_1st_lineage), # to quantify asymmetricity
                 minmax_ratio_number_of_1st_lineage        = min(n_samples_from_1st_branch)/max(n_samples_from_1st_branch) # to quantify asymmetricity of 1st lineage
                 ) 
  
  # draw tree
  if(draw.tree) {
    tr             <- read.newick(textConnection(result$tree))
    p              <- ggtree(tr, ladderize=T) + geom_tiplab()
    edge           <- data.frame(tr$edge, edge_length=tr$edge.length)
    colnames(edge) <- c("parent", "node", "edge_length")
    p              <- p %<+% edge + geom_label(aes(x=branch, label=edge_length), na.rm = TRUE)
    print(p)
  }
  
  return(result)
}


args     <- commandArgs(trailingOnly=TRUE)
iter     <- as.integer(args[1])
sampling <- as.integer(args[2])
min_initial_rate <- as.numeric(args[3])
max_initial_rate <- as.numeric(args[4])
min_rate <- as.numeric(args[5])
max_rate <- as.numeric(args[6])
out      <- args[7]

cat("initial_rate\trate\tat\tn\tmultifurcation_score\tnumber_of_multifurcation\tnumber_of_bifurcation\tfrom_first_cell\tfrom_second_cell\tnumber_of_1st_lineage\tcontribution_of_1st_lineage\tdistance_1st_lineage\tnumber_of_shared_mutations\tbranch_length_of_1st_lineage\tnumber_of_shared_mutations_in_1st_lineage\tminmax_ratio_1st_shared_mutations\tminmax_ratio_number_of_1st_lineage\ttree\n", file=out, append=TRUE)

      for(i in 1:iter) {

        mutation_rate <- runif(1, min=min_rate, max=max_rate)
        initial_rate  <- runif(1, min=min_initial_rate, max=max_initial_rate)
        at            <- sample(2:7, 1)
        n             <- sample(1:(2^at/2), 1)

        print(sprintf("initial rate: %.1f, rate: %.1f, at: %d, n: %d, iter: %d", initial_rate, mutation_rate, at, n, i))
        
        cells <- Embryogenesis(Cell(barcode  = "0",
          mrate_nuclear = mutation_rate, 
          nuclear  = tibble(id=character())), generation = at + 10, mrate_1st_div = initial_rate)

        result <- Simulate_sequencing(cells, at=at, n=n, sampling=sampling, draw.tree=FALSE)

        if(result$from_first_cell == 0 | result$from_second_cell == 0) ratio <- 0
        else {
          ratio <- result$from_first_cell/result$from_second_cell
          if(ratio < 1) ratio <- 1/ratio
        }
        
        cat(sprintf("%.3f\t%.3f\t%d\t%d\t%.3f\t%d\t%d\t%d\t%d\t%d\t%s\t%.3f\t%d\t%s\t%d\t%.3f\t%.3f\t%s\n", 
          initial_rate, 
          mutation_rate,
          at, 
          n, 
          result$multifurcation_score,
          result$number_of_multifurcation,
          result$number_of_bifurcation,
          result$from_first_cell,
          result$from_second_cell,
          result$number_of_1st_lineage,
          result$contribution_of_1st_lineage,
          result$distance_1st_lineage,
          result$number_of_shared_mutations,
          result$branch_length_of_1st_lineage,
          result$number_of_shared_mutations_in_1st_lineage,
          result$minmax_ratio_1st_shared_mutations,
          result$minmax_ratio_number_of_1st_lineage,
          result$tree),
        file=out, append=TRUE)
      }

