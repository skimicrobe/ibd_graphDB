#######################################################
# Goal: Formatting dataframe in Neo4j-import format   #
# Author: Suyeon Kim                                  #
# Creation Date: Dec-05-2021                          #
# Last updated: Feb-25-2021                           #
#######################################################

# In this function, we create unique Ids for entity(Feature)
create_Entity_uniqueIds<- function(entity, prefix, final_propertyTab){
  print(paste("Generating unique", entity, "Ids"))
  preFix<- prefix
  suFFix<- seq(1:nrow(final_propertyTab))
  entity_ids<- paste(preFix, suFFix, sep="")
  return(entity_ids)
}

# In this function, we want to make consistent column names, 
# containing special characters
edit_specialCharacters<- function(participant_property) {
  print("Rename the special chatercter in columnNames")
  
  # Filtering column names for Participant Node dataframe. 
  colnames(participant_property)<- gsub("X","",colnames(participant_property))
  colnames(participant_property)<- gsub("\\.+","\\_", colnames(participant_property))
  
  # Use `gsub` function: special characters only at the end of the string!
  colnames(participant_property)<- gsub('\\_$','',colnames(participant_property))
  colnames(participant_property)<- gsub('^\\_',"",colnames(participant_property))
  
  return(participant_property)
}

# In this function, we want to reformat 'Node' dataframe as importable format 
# in the Neo4j      
reformat_Node_file<- function(participant_node, old_col, new_col, node_label){
  print("begin reformatting node datafile")
  format_dat<- participant_node
  names(format_dat)[names(format_dat) == old_col]<- new_col
  print(dim(format_dat))
  # Add label column for labeling node in neo4j 
  format_dat$`:LABEL`<- node_label
  print(dim(format_dat))
  return(format_dat)
}

# In this function, we want to reformat 'Edge' dataframe as importable format 
# in the Neo4j 
reformat_Edge_file<- function(edge_input, start_col, end_col, edge_label) {
  print("begin reformatting edge datafile")
  # Reorder columns in target and source relationship
  move_col<- c(start_col, end_col)
  reorder_edgeDat<- edge_input[c(move_col, setdiff(names(edge_input), move_col))]
  
  # Rename columns 
  names(reorder_edgeDat)[names(reorder_edgeDat) == start_col]<- ":START_ID"
  names(reorder_edgeDat)[names(reorder_edgeDat) == end_col]<- ":END_ID"
  
  # Add Label column for labeling edge in neo4j  
  reorder_edgeDat$`:TYPE`<- edge_label
  return(reorder_edgeDat)
}
