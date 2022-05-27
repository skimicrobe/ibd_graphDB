#----------------------------------------------#
# Goal: Formatting dataframe in Neo4j format   #
# Name: Suyeon Kim                             #
# Creation Date: Dec-05-2021                   #
# Last updated: Feb-25-2021                    #
#----------------------------------------------#

create_Entity_uniqueIds<- function(entity, prefix, final_propertyTab){
  print(paste("Generating unique", entity, "Ids"))
  preFix<- prefix
  suFFix<- seq(1:nrow(final_propertyTab))
  entity_ids<- paste(preFix, suFFix, sep="")
  return(entity_ids)
}

## Make consistent column names (special characters)
edit_specialCharacters<- function(participant_property) {
  print("Rename the special chatercter in columnNames")
  
  ## Filtering column names for Participant Node dataframe. 
  colnames(participant_property)<- gsub("X","",colnames(participant_property))
  colnames(participant_property)<- gsub("\\.+","\\_", colnames(participant_property))
  
  # gsub special characters only at the end of the string!
  colnames(participant_property)<- gsub('\\_$','',colnames(participant_property))
  colnames(participant_property)<- gsub('^\\_',"",colnames(participant_property))
  
  return(participant_property)
}

## Reformat Node data as importable format in the Neo4j      
#reformat_Node_file<- function(participant.node, col_removed, old_col, new_col, node_label) {
   
#  format_dat<- participant.node %>% dplyr::select(-col_removed)
#  names(format_dat)[names(format_dat) == old_col]<- new_col
#  print(dim(format_dat))
#  format_dat$`:LABEL`<- node_label
#  print(dim(format_dat))
  
#  return(format_dat)
#}
reformat_Node_file<- function(participant_node, old_col, new_col, node_label){
  format_dat<- participant_node
  names(format_dat)[names(format_dat) == old_col]<- new_col
  print(dim(format_dat))
  format_dat$`:LABEL`<- node_label
  print(dim(format_dat))
  return(format_dat)
}

## Reformat Edge data as importable format in the Neo4j 
reformat_Edge_file<- function(edge_input, start_col, end_col, edge_label) {
  print("begin reformatting edge file")
  
  ## reorder columns in neo4j format
  move_col<- c(start_col, end_col)
  reorder_edgeDat<- edge_input[c(move_col, setdiff(names(edge_input), move_col))]
  
  ## rename columns 
  names(reorder_edgeDat)[names(reorder_edgeDat) == start_col]<- ":START_ID"
  names(reorder_edgeDat)[names(reorder_edgeDat) == end_col]<- ":END_ID"
  
  ## Add Label column for representation in neo4j 
  reorder_edgeDat$`:TYPE`<- edge_label
  return(reorder_edgeDat)
}
