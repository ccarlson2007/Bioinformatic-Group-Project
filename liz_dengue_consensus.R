setwd("~/Desktop/EECBData/BioinfoRData/")
library(seqinr)

dengue_basic <- read.fasta("DengueVirus1.fasta_pruned.mu.trim05")
number_of_seqs <- length(dengue_basic)
dengue_align <- read.alignment("denguesmall.txt", format = "fasta", forceToLower = T)
dengue_consensus <- seqinr :: consensus(dengue_align, method = "majority")
dengue_consensus_matrix <- seqinr :: consensus(dengue_align, method = "profile")

consensus_length <- length(dengue_consensus)
number_column <- seq(1, consensus_length)

Dengue_DF <- data.frame("num" = number_column, "MeanFreq" = 0, "WTnt" = dengue_consensus)

base_count <- ncol(dengue_consensus_matrix)

for(x in 1:base_count){
  current_base <- dengue_consensus[x]
  current_matrix_base_count <- dengue_consensus_matrix[,x]
  ts_count <- 0
  
  if(current_base == "a"){
    ts_count <- current_matrix_base_count[["g"]]
  }
  if(current_base == "g"){
    ts_count <- current_matrix_base_count[["a"]]
  }
  if(current_base == "c"){
    ts_count <- current_matrix_base_count[["t"]]
  }
  if(current_base == "t"){
    ts_count <- current_matrix_base_count[["c"]]
  }
  
  Dengue_DF[x, 2] <- ts_count/number_of_seqs
}




CpG_finder <- function(new_virus_data){
  virus_data <- new_virus_data
  WTnt <- virus_data$WTnt
  data_length <- nrow(virus_data)
  new_CpG <- vector(mode = "numeric", length = data_length)
  
  for(x in 1:data_length){
    current_nuc <- toupper(WTnt[x])
    current_neighbor <- toupper(WTnt[x + 1])
    
    if(current_nuc == "T"){
      if(current_neighbor == "G"){
        new_CpG[x] = 1
      }
      else{}
    }
    
    if(current_nuc == "C"){
      if(current_neighbor == "A"){
        new_CpG[x+1] = 1
      }
      else{}
    }
    else{}
  }
  
  virus_data$new_CpG <- new_CpG
  
  return(virus_data)
}

CpG_finder(Dengue_DF)
