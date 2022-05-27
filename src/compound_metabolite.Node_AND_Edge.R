################################################################################
# Script: compound_metabolite.Node_AND_Edge.R                                  #
# Author: Suyeon Kim                                                           #
# Created: April 16, 2022                                                      #
# Last edited: April 16, 2022                                                  #
#                                                                              #
# Goal: The purpose of this script is to create a nodes and edges file in      #
# neo4j bulk importer format (neo4j-admin import tool).                        #
# Input: 1) HMP2_metabolomics.csv                                              # 
#        2) sample.node.[output argument].csv                                  #  
#                                                                              #
# Output: compound & metabolite nodes, four sample_compound edge list          #
#         (by technogloy), and  compound-metabolite edge                       #
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

# In this function, we want to remove if the columns have all NAs. 
Identify_missingColumns <- function(metabolite_dat) {
  all_missingCol <-
    colnames(metabolite_dat)[colSums(is.na(metabolite_dat)) == 
                               nrow(metabolite_dat)]
  print(length(all_missingCol))
  # Remove the rows only if it has all NAs.
  if (length(all_missingCol) != 0) {
    filt_metabolite <- metabolite_dat %>% select(-all_missingCol)
  }
  else{
    filt_metabolite <- metabolite_dat
  }
  return(filt_metabolite)
}

# In this function, we want to subset abundance data based on the 
# metabolomics technology
subset_data<- function(methodType=method_type, filt_metabolite){
  print(paste("Subsetting data - technology name:", print(methodType)))
  subset_dat<- filt_metabolite[grep(methodType, filt_metabolite$Method),]
  return(subset_dat)
}

# In this function, we want to replace the compound name and sample ID with 
# the created sample Ids and compound Ids. 
get_CompoundIds<- function(abundance, entity_Node, entity, sample_info){
  print(dim(abundance))
  # Merge abudance data with Compound Node dataframe based on entity
  merged_cpId<- merge(abundance, entity_Node, by=c(entity))
  filter_cols<- merged_cpId %>% select(-c(entity))
  matched_sampleIndex<- 
    sample_info[match(colnames(filter_cols)[7:ncol(filter_cols)-1], 
                      sample_info$External_ID),1]
  colnames(filter_cols)[7:ncol(filter_cols)-1]<- matched_sampleIndex
  
  # Create abundance table only with created compound Ids and sample Ids
  final_ab<- cbind("compound_ID"=filter_cols$compound_ID, 
                   filter_cols[,-ncol(filter_cols)])
  print("**********Before Reshape final_ab *******")
  print(dim(final_ab))
  
  return(final_ab)
}

# In this function, we convert dataframe from wide to long format. 
reshape_CompoundInfo<- function(new_abundance, by_melt){
  reshape_dat<- melt(new_abundance, id.vars=c(by_melt))
  final_ab<- data.frame(cbind("s_ID:ID"=as.character(reshape_dat$variable), 
                              "compound_ID:ID"= reshape_dat$compound_ID, 
                              "Abundance:double"=reshape_dat$value,
                              "Method"=reshape_dat$Method, 
                  "Pooled_QC_sample_CV:double"=reshape_dat$Pooled.QC.sample.CV,
                              "mz:double"= reshape_dat$m.z, 
                              "RT:double"=reshape_dat$RT,
              "HMDB_Representative_ID"= reshape_dat$HMDB...Representative.ID.),
                        check.names = FALSE)
  print("**********After reshape final_ab and entity_Node dimensions *******")
  print(dim(final_ab))
  
  return(final_ab)
}

################################################################################
# Read the Input file into R.                                                  #
################################################################################
# This should result in a dataframe called 'metabolite_dat' that contains 
# 81867 observations of 553 variables. 
metabolite_dat<- read.csv(paste0(InFile), quote='"')

# Remove If columns have all NAs. 
filt_metabolite<- Identify_missingColumns(metabolite_dat)

# Lets read the 'sample node' file into R  
sample_N<- read.csv(sample_node, header=T, sep=",", check.names=FALSE)
sample_info<- sample_N %>% select(c("s_ID:ID", "External_ID", "data_type"))

################################################################################
# Main Program: Create a Compound Node                                         #
# Input: metabolite_dat and sample_info                                        #
# Return: Compound_node                                                        #
# Note: Required source file called "formatting_data_neo4j.R"                  #
################################################################################
source(paste0(sourceDir,"/","formatting_data_neo4j.R"))

# Subset a column called 'Compound' in filt_metabolite dataframe. 
compound_info<- filt_metabolite %>% select("Compound")

# We want to create and add 'unique_ids' to a Compound node. 
compound_uniqIds<- create_Entity_uniqueIds(entity="Compound", prefix="compound", 
                                           final_propertyTab=compound_info)

# Add created above unique-Ids to the Compound node.
Compound_node<- cbind("compound_ID"=compound_uniqIds, compound_info)

################################################################################
# Main Program: Create a Metabolite Node                                       #
# Input: metabolite_dat and sample_info                                        #
# Return: Metabolite node                                                      #                       
################################################################################

# Subset a column called 'Metabolite' in filt_metabolite dataframe. 
metabolite_info<- filt_metabolite %>% select("Metabolite")
uniq_metabolite<- rename(data.frame(unique(metabolite_info$Metabolite)),
            Metabolite=colnames(data.frame(unique(metabolite_info$Metabolite))))

# We want to create and add 'unique_ids' to a Metabolite node. 
metabolite_uniqIds<- create_Entity_uniqueIds(entity="Metabolite", 
                        prefix="metabolite", final_propertyTab=uniq_metabolite)

# Add created above unique-Ids to the Metabolite node.
Metabolite_node<- cbind("metabolite_ID"=metabolite_uniqIds, uniq_metabolite)

################################################################################
# Main Program: Create a Sample-Compound edge/relationship.                    #
# Input: Compound_node, sample_info, and filt_metabolite                       #
# Return: Four sample_compound edges                                           #
#    - HILIC-pos : Positive ion mode of polar metabolites                      #
#    - HILIC-neg : Negative ion mode of polar metabolites                      #
#    - C18-neg: Intermediate metabolites                                       #
#    - C8-pos: Lipid metabolites                                               #
################################################################################
# Classify 'filt_metabolite' dataframe based on the different technology used 
# for the metabolite measurements. 
method_Type<- c("C18-neg", "C8-pos", "HILIC-neg", "HILIC-pos")
interM_mb<- subset_data(methodType=method_Type[1], filt_metabolite)
lipid_mb<- subset_data(methodType=method_Type[2], filt_metabolite)
polar_neg<- subset_data(methodType = method_Type[3], filt_metabolite)
polar_pos<- subset_data(methodType = method_Type[4], filt_metabolite)

# Let's create a 'Sample-Compound Edge'. Reshape the 'filt_metabolite' 
# dataframe in Neo4j format. 
# Subset metabolomics samples in sample_info dataframe. 
mb_sample<- sample_info[sample_info$data_type == "metabolomics",]

# Technology1 : Negative ion mode of interM metabolites (C18-neg)
filt_interM<- interM_mb %>% select(-"Metabolite")
new_interM_ab<- get_CompoundIds(filt_interM, entity_Node = Compound_node, 
                        entity="Compound", sample_info = mb_sample) 
interM_reshp<- reshape_CompoundInfo(new_interM_ab,
          by_melt=c("compound_ID", "Method", "Pooled.QC.sample.CV", "m.z", "RT",
                    "HMDB...Representative.ID."))

# Technology2: lipids metabolites (C8-pos)
filt_lipid<- lipid_mb %>% select(-"Metabolite")
new_lipid_ab<- get_CompoundIds(filt_lipid, entity_Node = Compound_node, 
                                entity="Compound", sample_info = mb_sample) 
lipid_reshp<- reshape_CompoundInfo(new_lipid_ab,
          by_melt=c("compound_ID", "Method", "Pooled.QC.sample.CV", "m.z", "RT",
                                              "HMDB...Representative.ID."))

# Technology3: negative ion mode polar metabolites (HILIC-neg)
filt_polarNeg<- polar_neg %>% select(-"Metabolite")
new_polarNeg_ab<- get_CompoundIds(filt_polarNeg, entity_Node = Compound_node, 
                               entity="Compound", sample_info = mb_sample) 
polarNeg_reshp<- reshape_CompoundInfo(new_polarNeg_ab,
          by_melt=c("compound_ID", "Method", "Pooled.QC.sample.CV", "m.z", "RT",
                                             "HMDB...Representative.ID."))

# Technology4: positive ion mode polar metabolites (HILIC-pos)
filt_polarPos<- polar_pos %>% select(-"Metabolite")
new_polarPos_ab<- get_CompoundIds(filt_polarPos, entity_Node = Compound_node, 
                                  entity="Compound", sample_info = mb_sample) 
polarPos_reshp<- reshape_CompoundInfo(new_polarPos_ab,
          by_melt=c("compound_ID", "Method", "Pooled.QC.sample.CV", "m.z", "RT",
                                                "HMDB...Representative.ID."))

################################################################################
# Main Program: Create a Compound-Metabolite edge/relationship.                #
# Input: Metabolite_node, Compound_node                                        #
################################################################################
# Let's subset "Compound" and "Metabolite" columns in filt_metabolite dataframe.
compound_metabolite<- filt_metabolite %>% select(c("Compound", "Metabolite"))

# Merge "compound_metabolite and Compound_node" dataframes by 'Compound' column. 
matched_compoundId<- merge(compound_metabolite,Compound_node,by=c("Compound"))
# Merge a created previous dataframe called 'matched_compoundId' and 
# 'Metabolite_node' by ' Metabolite' column. 
matched_metaboliteId<- merge(matched_compoundId,Metabolite_node,
                             by=c("Metabolite"))
# We want to create a 'compound_metabolite' edge. 
cp_mb_Edge<- matched_metaboliteId %>% select(c("compound_ID", "metabolite_ID"))

###############################################################################
# Main Program: Reformat data in neo4j-admin import format and Save Output.   #
# Input: compound_node, metabolite_node, cp_mb_Edge                           #
# Return: neo4j_comPound_Node, neo4j_metabolite_Node, neo4j_[4tech]_Edge, and #
#         neo4j_cp_Mb_Edg                                                     #
# Note: Required source file called "formatting_data_neo4j.R"                 #
###############################################################################
# Process Compound_node in neo4j-import format 
neo4j_comPound_Node<- reformat_Node_file(Compound_node, old_col="compound_ID", 
                                         new_col="compound_ID:ID", 
                                         node_label="Compound")

# Process Metabolite_node in neo4j-import format 
neo4j_metabolite_Node<- reformat_Node_file(Metabolite_node, 
                                           old_col="metabolite_ID", 
                                           new_col="metabolite_ID:ID", 
                                           node_label="Metabolite")

# Process Sample-Compound (interM) edge/relationship in neo4j-import format 
neo4j_interM_Edge<- reformat_Edge_file(interM_reshp, start_col="s_ID:ID", 
                                       end_col="compound_ID:ID", 
                       edge_label="HAS_NEGATIVE_ION_MODE_OF_INTERM_METABOLITES")

# Process Sample-Compound (lipid) edge/relationship in neo4j-import format 
neo4j_lipid_Edge<- reformat_Edge_file(lipid_reshp, start_col="s_ID:ID", 
                                      end_col="compound_ID:ID",
                                      edge_label="HAS_LIPIDS_METABOLITES")

# Process Sample-Compound (Polar-Neg) edge/relationship in neo4j-import format
neo4j_polarNeg_Edge<- reformat_Edge_file(polarNeg_reshp, start_col="s_ID:ID", 
                                         end_col="compound_ID:ID",
                        edge_label="HAS_NEGATIVE_ION_MODE_OF_POLAR_METABOLITES")

# Process Sample-Compound (Polar-Pos) edge/relationship in neo4j-import format 
neo4j_polarPos_Edge<- reformat_Edge_file(polarPos_reshp, start_col="s_ID:ID", 
                                         end_col="compound_ID:ID",
                        edge_label="HAS_POSITIVE_ION_MODE_OF_POLAR_METABOLITES")

# Process Compound-Metabolite edge/relationship in neo4j-import format 
neo4j_cp_Mb_Edge<- reformat_Edge_file(cp_mb_Edge, start_col="compound_ID", 
                                   end_col="metabolite_ID",edge_label="IS_FROM")

# We need to substitute NA with empty as Neo4 recognize NA as a string data type 
neo4j_comPound_Node[is.na(neo4j_comPound_Node)]<-""
neo4j_metabolite_Node[is.na(neo4j_metabolite_Node)]<- ""
neo4j_interM_Edge[is.na(neo4j_interM_Edge)]<- ""
neo4j_lipid_Edge[is.na(neo4j_lipid_Edge)]<- ""
neo4j_polarNeg_Edge[is.na(neo4j_polarNeg_Edge)]<- ""
neo4j_polarPos_Edge[is.na(neo4j_polarPos_Edge)]<- ""
neo4j_cp_Mb_Edge[is.na(neo4j_cp_Mb_Edge)]<- ""

# Write the Compound node Output file to a .csv file 
fileN1= paste("Compound", "node", outName, "csv", sep=".")
write.csv(neo4j_comPound_Node, fileN1, row.names=F, quote=T, sep="\t")

# Write the Compound node Output file to a .csv file 
fileN2= paste("Metabolite", "node", outName, "csv", sep=".")
write.csv(neo4j_metabolite_Node, fileN2, row.names=F, quote=T, sep="\t")

# Write the sample-InterMediate Compound file to a .csv file 
fileN3= paste("sample_InterMediate_cp","edge", outName, "csv", sep=".")
write.csv(neo4j_interM_Edge, fileN3, row.names=F, quote=F, sep=",")

# Write the sample-lipid Compound file to a .csv file 
fileN4= paste("sample_lipid_cp","edge", outName, "csv", sep=".")
write.csv(neo4j_lipid_Edge, fileN4, row.names=F, quote=F, sep=",")

# Write the sample-PolarNegative Compound file to a .csv file 
fileN5= paste("sample_polarNeg_cp","edge", outName, "csv", sep=".")
write.csv(neo4j_polarNeg_Edge, fileN5, row.names=F, quote=F, sep=",")

# Write the sample-PolarPositive Compound file to a .csv file 
fileN6= paste("sample_polarPos_cp","edge", outName, "csv", sep=".")
write.csv(neo4j_polarPos_Edge, fileN6, row.names=F, quote=F, sep=",")

# Write the Compound-Metabolite file to a .csv file 
fileN7= paste("compound_metabolite","edge", outName, "csv", sep=".")
write.csv(neo4j_cp_Mb_Edge, fileN7, row.names=F, quote=F, sep=",")
