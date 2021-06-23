# VAF를 mpileup 값으로 대체

convertVAF_from_mpileup <- function(df_varmat, path_to_mp) {
  
  samples <- df_varmat %>% select(-(Variant:ExonicFunc_ensGene)) %>% colnames()
  
  for(s in samples) {
    
    df_mp <- read_tsv(str_c(path_to_mp, s, ".mp")) 
    
    for(i in 1:nrow(df_varmat)) {
      
      if(str_length(df_varmat$REF[i]) > str_length(df_varmat$ALT[i])) { # deletion
        ALT <- str_c("-", str_sub(df_varmat$REF[i], start = str_length(df_varmat$ALT[i]) + 1))
        
        indel <- gsub("[\\{\\}]", "", df_mp$Indel[which(df_mp$POS==df_varmat$POS[i])])
        indel <- gsub(":", "=", indel)
        indel <- eval(parse(text=paste("list(", indel, ")")))
          
        upper <- if(!is.null(indel[[toupper(ALT)]])) indel[[toupper(ALT)]] else 0
        lower <- if(!is.null(indel[[tolower(ALT)]])) indel[[tolower(ALT)]] else 0
        
        VAF <- (upper + lower)/df_mp$depth[which(df_mp$POS==df_varmat$POS[i])] * 100
        df_varmat[i, s] <- VAF
      } else if(str_length(df_varmat$REF[i]) < str_length(df_varmat$ALT[i])) { # insertion
        ALT <- str_c("+", str_sub(df_varmat$ALT[i], start = str_length(df_varmat$REF[i]) + 1))
        
        indel <- gsub("[\\{\\}]", "", df_mp$Indel[which(df_mp$POS==df_varmat$POS[i])])
        indel <- gsub(":", "=", indel)
        indel <- eval(parse(text=paste("list(", indel, ")")))
        
        upper <- if(!is.null(indel[[toupper(ALT)]])) indel[[toupper(ALT)]] else 0
        lower <- if(!is.null(indel[[tolower(ALT)]])) indel[[tolower(ALT)]] else 0
        
        VAF <- (upper + lower)/df_mp$depth[which(df_mp$POS==df_varmat$POS[i])] * 100
        df_varmat[i, s] <- VAF      
        
      } else if(str_length(df_varmat$REF[i]) == str_length(df_varmat$ALT[i])) {  # substitution

        VAF <- as.numeric(df_mp[which(df_mp$POS==df_varmat$POS[i]), df_varmat$ALT[i]]) / as.numeric(df_mp[which(df_mp$POS==df_varmat$POS[i]), "depth"]) * 100
        df_varmat[i, s] <- VAF
        #print(sprintf("%s %d %s %s %s %s %f", s, df_varmat$POS[i], df_varmat$REF[i], df_varmat$ALT[i], df_mp[which(df_mp$POS==df_varmat$POS[i]), df_varmat$ALT[i]], df_mp[which(df_mp$POS==df_varmat$POS[i]), "depth"],VAF))
        #print(df_mp[which(df_mp$POS==df_varmat$POS[i]),])
      }
    }
  }
  
  invisible(df_varmat)
}

getVarMat("/home/users/chrono0707/projects/02_DB/02_DB6/04_MT_merged/02_Filtered_bam_BQ10_MQ10/01_annotated/01_baminfo") %>% convertVAF_from_mpileup("/home/users/chrono0707/projects/02_DB/02_DB6/08_Mpileup/01_Q10q10/") %>% write_tsv("./output/temp.tsv")
#df_temp2 <- convertVAF_from_mpileup(df_temp, "/home/users/chrono0707/projects/02_DB/97_merged_mpileups/")
