CpG_finder <- function(virus_data){
  #reads data into function from WTnt column
  WTnt <- virus_data$WTnt
  
  #gets the length of the data file for use in a loop
  data_length <- nrow(virus_data)
  
  #creates a new column in virus_data called makesCpG
  virus_data$makesCpG <- 0
  
  #loop that determines if a CpG could occur due to mutation at each spot in the list of nucleotides (WTnt)
  #loops from row 1 to the the last row of the column WTnt
  for(x in 1:(data_length-1)){
    
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
        virus_data$makesCpG[x] = 1
      }
      #otherwise, nothing is changed
      else{}
    }
    #repeat the same loop as before, but modified slightly to look for CA pairs
    if(current_nuc == "C"){
      if(current_neighbor == "A"){
        #here, instead of changing the position of the current nuc to a 1, we change the neighbor's position because that is where the mutation would be, hence x+1
        virus_data$makesCpG[x] = 1
      }
      else{}
    }
    else{}
    #print(c("Data length is ", nrow(virus_data), "List length is ", length(makesCpG)))
  }
  
  #when the entire function is done, return the amended data file
  return(virus_data)
}