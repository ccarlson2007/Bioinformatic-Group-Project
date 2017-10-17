#Importing data, selecting desired WTnt column, creating a new variable of the same length
HIV_data <- read.csv("HIV_Data.csv")
WT <- HIV_data$WTnt
WT_length <- length(WT)

#Making a vector, and using it to create a new data frame
makesCpG <- vector(mode = "numeric", length = WT_length)
makesCpG
new_CpG <- data.frame(WT, makesCpG)
new_CpG
length(makesCpG)

#Function
CpG_finder <- function(input_seq) {
  for(x in 1:WT_length) {
    nuc_pos_1 <- WT[x]
    nuc_pos_2 <- WT[x+1]
    if(nuc_pos_1 == "t") {
      if(nuc_pos_2 == "g") {
        makesCpG[x] <- 1
      }
      else {}
    }
    if(nuc_pos_1 == "c") {
      if(nuc_pos_2 == "a") {
        makesCpG[x+1] <- 1
      }
      else {}
    }
  }
  new_CpG$makesCpG <- makesCpG
  return(new_CpG)
}

CpG_finder(HIV_data)
