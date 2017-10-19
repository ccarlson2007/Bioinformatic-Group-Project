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
  #reads data into function as CSV file
  virus_data <- new_virus_data
  
  #singles out column of nucleotides -- this assumes that the new datafile has the same column headers such as WTnt as the HIV file does
  WTnt <- virus_data$WTnt
  
  #gets the length of the data file for use in a loop
  data_length <- nrow(virus_data)
  
  #creates an empty vector (like a list) of the same length as the data file to be used to record the results of the loop below
  new_CpG <- vector(mode = "numeric", length = data_length)
  
  #loop that determines if a CpG could occur due to mutation at each spot in the list of nucleotides (WTnt)
  #loops from row 1 to the the last row of the column WTnt
  for(x in 1:data_length){
    
    #assigns a name (current_nuc) to the nucleotide at row x in WTnt and makes the nucleotide capitalized, in case the data uses lower case letters
    current_nuc <- toupper(WTnt[x])
    
    #assigns a name (current_neighbor) to the nucleotide in the next row down in WTnt and makes the nucleotide capitalized, in case the data uses lower case letters
    current_neighbor <- toupper(WTnt[x + 1])
    
    #begins a loop that compares the two nucleotides that were just isolated
    #this if statement is only evaluated if the current nucleotide is a T
    if(current_nuc == "T"){
      #this if statement is only evaluated if the neighbor is a G
      if(current_neighbor == "G"){
        #if the current nucleotide is a T and the neighbor is a G (thus a mutation can cause a CpG), then the of the current nucleotide in the list of 0s we made is changed to a 1
        new_CpG[x] = 1
      }
      #otherwise, nothing is changed
      else{}
    }
    #repeat the same loop as before, but modified slightly to look for CA pairs
    if(current_nuc == "C"){
      if(current_neighbor == "A"){
        #here, instead of changing the position of the current nuc to a 1, we change the neighbor's position because that is where the mutation would be, hence x+1
        new_CpG[x+1] = 1
      }
      else{}
    }
    else{}
    #print(c("Data length is ", nrow(virus_data), "List length is ", length(new_CpG)))
  }
  
  #append the list of 0s and 1s we created to the end of the data file we imported
  virus_data$new_CpG <- new_CpG
  
  #when the entire function is done, return the amended data file
  return(virus_data)
}

CpG_finder(Dengue_DF)
