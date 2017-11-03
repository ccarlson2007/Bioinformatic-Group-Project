setwd("~/Desktop/EECB Classes Fall 2017/Bioinformatics Fall 2017/bioinfoRliz/")

virus_dataframe <- function(fasta_file){
  library(seqinr)

  dengue_basic <- read.fasta(fasta_file)
  number_of_seqs <- length(dengue_basic)
  
  dengue_align <- read.alignment(fasta_file, format = "fasta", forceToLower = T)
  dengue_consensus <- seqinr :: consensus(dengue_align, method = "majority")
  dengue_consensus_matrix <- seqinr :: consensus(dengue_align, method = "profile")
  
  consensus_length <- length(dengue_consensus)
  
  number_column <- seq(1, consensus_length)
  
  Dengue_DF <- data.frame("num" = number_column, "MeanFreq" = 0, "WTnt" = dengue_consensus)
  
  for(x in 1:consensus_length){
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
  
  return(Dengue_DF)
}


CpG_finder <- function(new_virus_data){
  
  virus_data <- new_virus_data
  
  WTnt <- virus_data$WTnt
  
  data_length <- nrow(virus_data)
  
  makesCpG <- vector(mode = "numeric", length = data_length)
  
  for(x in 1:(data_length-1)){
    
    current_nuc <- toupper(WTnt[x])
    
    current_neighbor <- toupper(WTnt[x + 1])
    
    if(current_nuc == "T"){
      if(current_neighbor == "G"){
        makesCpG[x] = 1
      }
      else{}
    }
    if(current_nuc == "C"){
      if(current_neighbor == "A"){
        makesCpG[x+1] = 1
      }
      else{}
    }
    else{}
  }
  
  virus_data$makesCpG <- makesCpG
  
  return(virus_data)
}

test_virus <- virus_dataframe("DengueVirus1.fasta_pruned.mu.trim05")
test_virus <- CpG_finder(test_virus)

my_consensus <- as.vector(test_virus$WTnt, mode = "character")

MUTAA <- list()
WTAA <- list()

for(x in seq(1, length(my_consensus), 3)){
  first <- my_consensus[x]
  second <- my_consensus[x+1]
  third <- my_consensus[x+2]
  codon <- c(first, second, third)
  
  mutated_codon <- codon
  
  if(codon[1] == "a"){
    mutated_codon <- replace(x=mutated_codon, values=c("g", codon[2], codon[3]))
  }
  if(codon[1] == "g"){
    mutated_codon <- replace(x=mutated_codon, values=c("a", codon[2], codon[3]))
  }
  if(codon[1] == "c"){
    mutated_codon <- replace(x=mutated_codon, values=c("t", codon[2], codon[3]))
  }
  if(codon[1] == "t"){
    mutated_codon <- replace(x=mutated_codon, values=c("c", codon[2], codon[3]))
  }
  MUTAA[x] <- translate(mutated_codon)
  
  if(codon[2] == "a"){
    mutated_codon <- replace(x=mutated_codon, values=c(codon[1], "g", codon[3]))
  }
  if(codon[2] == "g"){
    mutated_codon <- replace(x=mutated_codon, values=c(codon[1], "a", codon[3]))
  }
  if(codon[2] == "c"){
    mutated_codon <- replace(x=mutated_codon, values=c(codon[1], "t", codon[3]))
  }
  if(codon[2] == "t"){
    mutated_codon <- replace(x=mutated_codon, values=c(codon[1], "c", codon[3]))
  }
  MUTAA[x+1] <- translate(mutated_codon)
  
  if(codon[3] == "a"){
    mutated_codon <- replace(x=mutated_codon, values=c(codon[1], codon[2], "g"))
  }
  if(codon[3] == "g"){
    mutated_codon <- replace(x=mutated_codon, values=c(codon[1], codon[2], "a"))
  }
  if(codon[3] == "c"){
    mutated_codon <- replace(x=mutated_codon, values=c(codon[1], codon[2], "t"))
  }
  if(codon[3] == "t"){
    mutated_codon <- replace(x=mutated_codon, values=c(codon[1], codon[2], "c"))
  }
  MUTAA[x+2] <- translate(mutated_codon)
}

for(x in seq(1, length(my_consensus) - 2, 3)){
  codon <- c(my_consensus[x], my_consensus[x+1], my_consensus[x+2])
  new_AA <- translate(codon)
  WTAA[x] <- new_AA
  WTAA[x+1] <- new_AA
  WTAA[x+2] <- new_AA
}

DF <- cbind(WTAA, MUTAA)

DF <- data.frame(cbind(WTAA, MUTAA))
DF$TypeOfSite <- "NA" 

for(x in seq(1, nrow(DF), 1)){
  if(DF[[x, "WTAA"]] == DF[[x, "MUTAA"]]){
    DF[[x, "TypeOfSite"]] <- "syn"
  }
  if(DF[[x, "WTAA"]] != DF[[x, "MUTAA"]]){
    DF[[x, "TypeOfSite"]] <- "nonsyn"
  }
  if(DF[[x, "WTAA"]] != DF[[x, "MUTAA"]] &&  DF[[x, "MUTAA"]] == "*"){
    DF[[x, "TypeOfSite"]] <- "stop"
  }
}

DF$bigAAChange <- 0

pos <- c("R", "H", "K")
neg <- c("D", "E")
unc <- c("S", "T", "N", "Q")
spe <- c("C", "U", "G", "P")
hyd <- c("A", "I", "L", "F", "M", "W", "Y", "V")

for(x in 1:nrow(DF)){
  if(DF[[x, "WTAA"]] %in% pos){
    wtaa_group <- "pos"
  }
  if(DF[[x, "WTAA"]] %in% neg){
    wtaa_group <- "neg"
  }
  if(DF[[x, "WTAA"]] %in% unc){
    wtaa_group <- "unc"
  }
  if(DF[[x, "WTAA"]] %in% spe){
    wtaa_group <- "spe"
  }
  if(DF[[x, "WTAA"]] %in% hyd){
    wtaa_group <- "hyd"
  }
  
  if(DF[[x, "MUTAA"]] %in% pos){
    mutaa_group <- "pos"
  }
  if(DF[[x, "MUTAA"]] %in% neg){
    mutaa_group <- "neg"
  }
  if(DF[[x, "MUTAA"]] %in% unc){
    mutaa_group <- "unc"
  }
  if(DF[[x, "MUTAA"]] %in% spe){
    mutaa_group <- "spe"
  }
  if(DF[[x, "MUTAA"]] %in% hyd){
    mutaa_group <- "hyd"
  }
  
  if(wtaa_group != mutaa_group){
    DF[[x, "bigAAChange"]] <- 1
  }
}

big_DF <- cbind(test_virus, DF)