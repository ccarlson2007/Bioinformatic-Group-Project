CpG_finder <- function(virus_data){
  WTnt <- virus_data$WTnt
  
  data_length <- nrow(virus_data)
  
  virus_data$makesCpG <- 0
  
  for(x in 1:data_length-1){
    current_nuc <- toupper(WTnt[x])
    current_neighbor <- toupper(WTnt[x + 1])
    
    if(current_nuc == "T"){
      if(current_neighbor == "G"){
        virus_data$makesCpG[x] = 1
      }
      else{}
    }
    if(current_nuc == "C"){
      if(current_neighbor == "A"){
        virus_data$makesCpG[x] = 1
      }
      else{}
    }
    else{}
  }
  
  return(virus_data)
}