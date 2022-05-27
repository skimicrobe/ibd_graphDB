################################################################################
# Script: serum.Node_AND_Edge.R                                                #
# Author: Suyeon Kim                                                           #
# Created: April 17, 2022                                                      #
# Last edited: April 17, 2022                                                  #
#                                                                              #
# Goal: The purpose of this script is to create a Serum nodes and sample-serum #
#       edges file in neo4j bulk importer format.                              #
# Input: 1) hmp2_serology_Compiled_ELISA_Data.tsv                              #
#        2) sample.node.[output argument].csv                                  #
# Output: a Serum node and sample_serum edge list                              #
################################################################################
library(stringr)
library(dplyr)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)
print(args)

sourceDir <- args[1]
InFile<- args[2]
sample_node<- args[3]
outName<- args[4]

# This source file is to check the global environment to be 
# "stringsAsFactors=FALSE".
source(paste0(sourceDir,"/","getOption_setOption.R"))

# In this function, we want to replace the serum marker name and sample ID with 
# the created sample Ids and serum marker Ids. 
get_Ids<- function(abundance, entity_Node, entity, sample_info){
  print(dim(abundance))
  # Merge abudance data with serum marker Node dataframe based on entity
  merged_serumId<- merge(abundance, entity_Node, by=c(entity))
  filter_cols<- merged_serumId %>% select(-c(entity, "Serum"))
  matched_sampleIndex<- 
    sample_info[match(colnames(filter_cols[,-ncol(filter_cols)]), 
                      sample_info$External_ID),1]
  colnames(filter_cols)[1:ncol(filter_cols)-1]<- matched_sampleIndex
  # Create abundance table only with created gene-ids and sample-ids
  final_ab<- cbind("serum_ID"= filter_cols[,ncol(filter_cols)], 
                   filter_cols[,-ncol(filter_cols)])
  print("**********Before Reshape final_ab *******")
  print(dim(final_ab))
  
  return(final_ab)
}

# In this function, we convert dataframe from wide to long format.
reshape_dat<- function(new_abundance, by_melt){
  reshape_dat<- melt(new_abundance, id.vars=c(by_melt))
  final_ab<- data.frame(cbind("s_ID:ID"=as.character(reshape_dat$variable), 
                              "serum_ID:ID"= reshape_dat$serum_ID, 
                              "Abundance:double"=reshape_dat$value), 
                        check.names=FALSE)
  print("**********After reshape final_ab and entity_Node dimensions *******")
  print(dim(final_ab))
  
  return(final_ab)
}

###############################################################################
# Read the Input file into R.                                                 #
###############################################################################
# This should result in a dataframe called 'serum_dat' that contains 
# 14 observations of 212 variables. 
serum_dat<- read.csv(InFile, header=TRUE, sep="\t", check.names=FALSE)
colnames(serum_dat)[1]<- "Serum"

# Lets read the 'sample node' file into R  
sample_N<- read.csv(sample_node, header=T, sep=",", check.names=FALSE)
sample_info<- sample_N %>% select(c("s_ID:ID", "External_ID", "data_type"))

###############################################################################
# Main Program: Create a Serum Node                                           #
# Input: serum_dat                                                            #
# Return: serum_node data                                                     #
# Note: Required source files - formatting_data_neo4j.R                       #    
###############################################################################
source(paste0(sourceDir,"/","formatting_data_neo4j.R"))

# We want to remove the rows containing 'collection site', 'plate' and 
# 'SampleNum'. This cuts our file, now named 'filt_sero'.  
filt_sero<- subset(serum_dat, Serum != c("Site", "Plate", "Sample"))
filt_sero$Serum<- gsub("\\..*", "", filt_sero$Serum)

# We want to look at unique Serum_ID for the serum marker. 
serum_info<- rename(data.frame(unique(filt_sero$Serum)), 
                    Serum = colnames(data.frame(unique(filt_sero$Serum))))

# We want to create and add 'unique_ids' to a Serum Node. 
se_uniqIds<- create_Entity_uniqueIds(entity="Serum node", prefix="se", 
                                     final_propertyTab=serum_info)

# Add created above unique-Ids to the Serum node and remove part after . using 
# gsub() function. 
serum_node<- cbind("serum_ID"=se_uniqIds, serum_info)
serum_node$Serum<- gsub("\\..*", "", serum_node$Serum)

################################################################################
# Main Program: Create a Sample-SerumMarker edge/relationship.                 #
# Input: serum_node, sample_info, filt_sero                                    #
# Return: serum_reshp                                                          #
################################################################################
# Subset Serology samples in sample_info dataframe. 
serum_sample<- sample_info[sample_info$data_type == "serology",]
# In serum data, we identified missing information for one sample. Since we 
# cannot find their subject Id (participant Id) in metadata, this sample will be 
# removed from our input data.
new_filt_serum<- filt_sero %>% select(-c("220948"))

# Let's create a 'Sample-Serum Edge'. Reshape the 'filt_mGene' 
# dataframe in Neo4j format. 
serum_abundance<- get_Ids(new_filt_serum, entity_Node = serum_node, 
                          entity ="Serum", sample_info = serum_sample) 

serum_reshp<- reshape_dat(serum_abundance, by_melt = "serum_ID")

################################################################################
# Main Program: Reformat data in neo4j-admin import format and Save Output     #
# Input: serum_node, serum_reshp                                               #
# Return: edited_serum_Node, edited_serum_Edge                                 #
# Note: Required source file called "formatting_data_neo4j.R"                  #
################################################################################
# Process serum_node in neo4j-import format 
neo4j_serum_Node<- reformat_Node_file(serum_node, old_col="serum_ID", 
                                      new_col="serum_ID:ID", 
                                      node_label="Serum Marker")

# Process Sample-Serum Marker edge/relationship in neo4j-import format 
neo4j_serum_Edge<- reformat_Edge_file(serum_reshp, start_col="s_ID:ID", 
                                      end_col="serum_ID:ID", 
                                      edge_label="HAS_SERUM_ABUNDANCE")

# We need to substitute NA with empty as Neo4 recognize NA as a string data type 
neo4j_serum_Node[is.na(neo4j_serum_Node)]<-""
neo4j_serum_Edge[is.na(neo4j_serum_Edge)]<- ""

# Write the Serum node Output file to a .csv file 
fileN1= paste("Serum", "node", outName, "csv", sep=".")
write.csv(neo4j_serum_Node, fileN1, row.names=F, quote=T, sep="\t")

# Write the sample-Serum file to a .csv file 
fileN2= paste("sample_serum","edge", outName, "csv", sep=".")
write.csv(neo4j_serum_Edge, fileN2, row.names=F, quote=F, sep=",")
