###############################################################################
# Goal: Filtering MGX/MTX/MPX abundance data    
# Name: Suyeon Kim                             
# Creation Date: Mar-06-2022                   
# Last updated: Mar-06-2022                    
###############################################################################

# In this function, we want to remove any rows that do have species information
# and unclassified OR ungrouped. 
filt_Input<- function(raw_ab, entity){
  colnames(raw_ab)[1]<- entity
  print(paste("Original size of data", dim(raw_ab)))
  if(length(grep("g__", raw_ab[,1])) || 
     length(grep("unclassified", raw_ab[,1])) >0){
    raw_ab<- raw_ab[-grep("g__", raw_ab[,1]),]
    raw_ab<- raw_ab[-grep("unclassified", raw_ab[,1]),]
    print(paste("After 1st filetering ", dim(raw_ab)))
  }else if(length(grep("UNGROUPED", raw_ab[,1])) > 0){
    raw_ab<- raw_ab[-grep("UNGROUPED", raw_ab[,1]),]
    print(paste("After 2nd filtering", dim(raw_ab)))
  }
  else{
    raw_ab<- raw_ab
    print(paste("No change in size of data", dim(raw_ab)))
  }
  return(raw_ab)
}
