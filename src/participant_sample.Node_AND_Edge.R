################################################################################
# Script: participant_sample.Node_AND_Edge.R                                   #  
# Author: Suyeon Kim                                                           #
# Created: Feb 22, 2022                                                        #  
# Last edited: Feb 25, 2022                                                    #
#                                                                              #
# Goal: The purpose of this script is to create a Participant, sample nodes and# 
#     particiapnt_sample edges file in neo4j bulk importer format.             #
# Input: HMP2 metadata from ibdmdb database. This input file is required to be #
#        processed by 'merged_clinicalData.R' script.                          # 
# Output: a participant node, a sample node, and participant_sample edge list  #  
################################################################################

# Install the package if you don't installed yet. 
library(dplyr)
library(DataExplorer)

args = commandArgs(trailingOnly=TRUE)
print(args)

sourceDir <- args[1]
InFile<- args[2]
outName<- args[3]

# This source file is to check the global environment to be 
# "stringsAsFactors=FALSE".
source(paste0(sourceDir,"/","getOption_setOption.R"))

# In this function, we will identify the properties for a Participant node based
# on the list of discrete variables. 
get_discreteVal_for_ParticipantNode<- function(discrete_val, total_participant){
  # Get the current Participant.ID 
  participant_info<- discrete_val[discrete_val$Participant.ID == 
                                    total_participant,]
  participant_properties<- c() 
  
  # We are looking for discrete variables which has constant values for 
  # each of the participants. To do this, each of the discrete variables have to
  # check the following:
  #  condition 1) Get variables have an unique value (one-level). 
  #  condition 2) Get variables can be unique value but it wasn't unique in 
  #               the first place due to the empty character.  
  for( i in 1:ncol(participant_info)){
    if( length(unique(participant_info[,i])) == 1 || 
        length(unique(participant_info[,i][participant_info[,i] != ""])) == 1) {
      x<- colnames(participant_info)[i]
      participant_properties<- c(participant_properties, x)
    }
  }
  return(participant_properties)
}

# In this function, we will identify the properties for a Participant node based 
# on the list of continuous variables. 
get_contiVal_for_ParticipantNode<- function(continuous_val, total_participant){
  # Get the current Participant.ID
  participant_info<- continuous_val[continuous_val$Participant.ID == 
                                      total_participant,]
  participant_properties<- c()
  
  # Here, we are looking for continuous variables which has dynamic values for 
  # each of the participants. To do this, each of the continuous variables have 
  # to check the following: 
  #   condition 1) Get variables have all NA values for each of the participants.
  #                If values are all NA, we keep this as a candidate property 
  #                for a Participant node. 
  #   condition 2) Get variables have some values, if the minimum and
  #                maximum values are the same. 
  # We want to exclude the Participant.ID (column 2) and External.ID (column 1) 
  # in the 'continuous_val' dataframe. We note that a dataframe called 
  # 'continuous_val' that contains 4337 observations of 57 variables. 
  for( i in 3:ncol(participant_info)) {
    if(min(participant_info[,i],na.rm=T) == max(participant_info[,i],na.rm=T) || 
       length(which(is.na(participant_info[,i])) == "TRUE") == 
       nrow(participant_info)){    
      x<- colnames(participant_info)[i]
      participant_properties<- c(participant_properties, x)
    }
  }
  return(participant_properties)
}

# In this function, we will identify the properties for a Sample node based on 
# the list of discrete variables. 
get_discVal_for_SampleNode<- function(disVal_after_pNode, total_sample){
  
  # Get the current External.ID 
  sample_info<- disVal_after_pNode[disVal_after_pNode$External.ID == 
                                     total_sample,]
  sample_properties<- c()
  
  # We are looking for discrete variables which has constant values for 
  # each of the samples. To do this, each of the discrete variables have to
  # check the following:
  #  condition 1) Get variables have an unique value (one-level). 
  #  condition 2) Get variables can be unique value but it wasn't unique in 
  #               the first place due to the empty character.  
  for( i in 1:ncol(sample_info)){
    if(length(unique(sample_info[,i])) == 1 || 
       length(unique(sample_info[,i][sample_info[,i] != ""])) == 1) {
      x<- c(colnames(sample_info)[i]) 
      sample_properties<- cbind(sample_properties, x)
    }
  }
  return(sample_properties)
}

# In this function, we will identify the properties for a Sample node based on 
# the list of continuous variables. 
get_contiVal_for_SampleNode<- function(contiVal_after_pNode, total_sample){
  
  # Grep info based on current Participant ID
  sample_info<- contiVal_after_pNode[contiVal_after_pNode$External.ID == 
                                       total_sample,]
  sample_properties<- c()
  
  # Here, we are looking for continuous variables which has dynamic values for 
  # each of the participants. To do this, each of the continuous variables have 
  # to check the following: 
  #   condition 1) Get variables have all NA values for each of the participants.
  #                If values are all NA, we keep this as a candidate property 
  #                for a Participant node. 
  #   condition 2) Get variables have some values, if the minimum and
  #                maximum values are the same. 
  # We want to exclude the Participant.ID (column 2) and External.ID (column 1) 
  # in the 'contiVal_after_pNode' dataframe. We note that a dataframe called 
  # 'contiVal_after_pNode' that contains 4337 observations of 39 variables. 
  for( i in 3:ncol(sample_info)) {
    if(min(sample_info[,i], na.rm=T) == max(sample_info[,i], na.rm=T) || 
       length(which(is.na(sample_info[,i])) == "TRUE")==nrow(sample_info)){  
      x<- c(colnames(sample_info)[i]) 
      sample_properties<- c(sample_properties, x)
    }
  }
  return(sample_properties)
}

################################################################################
# Read the Input file into R.                                                  #
################################################################################
# This should result in a dataframe called 'meta_dat' that contains 
# 4337 observations of 491 variables. 
meta_dat<- read.csv(InFile, header=TRUE, sep=",", stringsAsFactors = FALSE, 
                    check.names=FALSE)

# Next, we remove columns for which contain all missing values. 
# This should result in a dataframe called 'filt_meta' that contains 
# 4337 observations of 449 variables. 
all_missingCol<- colnames(meta_dat)[colSums(is.na(meta_dat)) == nrow(meta_dat)]
filt_meta<- dplyr::select(meta_dat, -all_missingCol) 

################################################################################
# Data Pre-processing                                                          #
################################################################################

# Given our metadata, we separate the metadata('filt_meta') into discrete and 
# continuous data. This step will identify the columns from metadata as the 
# properties for a participant or sample nodes. 
two_data<- split_columns(filt_meta)
discrete_val<- two_data$discrete
continuous_val<- two_data$continuous

# As the columns called 'Participant.ID' and 'External.ID(Sample.ID)' were 
# removed in the previous step, we need to append these two columns on to the 
# dataframe called 'continuous_val'.
# Once this complete, we can begin to create a participant and sample node.
continuous_val<- cbind("Participant.ID"= filt_meta$Participant.ID, 
                       continuous_val)
continuous_val<- cbind("External.ID" = filt_meta$External.ID, continuous_val)

################################################################################
# Main Program: Create a Participant Node                                      #
# Input: filt_meta, continuous_val, and discrete_val                           #
# Return: This should result in a dataframe called 'participant_Node_prop'     #
#         that contains 131 observations of 269 variables.                     #
################################################################################

# First, let's identify which discrete variables can be the candidate properties
# for a Participant node. 
total_participant<- unique(filt_meta$Participant.ID)

# Create new empty lists to hold discrete variables properties for 
# 131 participants. 
list_disc_par_total<- c()

# This for loop will run until we identify the discrete variables that can be
# the candidate properties for 131 participants. 
for( sb in 1:length(total_participant)){
  participant_discreteVal<- get_discreteVal_for_ParticipantNode(discrete_val, 
                                      total_participant=total_participant[sb])
  length(participant_discreteVal)
  list_disc_par_total[[sb]]<- participant_discreteVal
}

# Get variables which appear across all of the participants in the list. 
colName.x.discrete.participant<- list_disc_par_total[[1]]
for( i in 2:length(list_disc_par_total)){
  colName.y.discrete.participant<- list_disc_par_total[[i]]
  colName.x.discrete.participant<- intersect(colName.x.discrete.participant, 
                                             colName.y.discrete.participant)
}

# Once this is complete, we can begin to identify the candidate properties for 
# a Participant node from continuous variables. To do this, let's create new 
# empty lists to hold continuous variables properties for 131 participants. 
list_conti_par_total<- c() 

# Similarly, let's identify the continuous variables that can be
# the candidate properties for 131 participants. 
for( sb in 1:length(total_participant)){
  participant_contiVal<- get_contiVal_for_ParticipantNode(continuous_val, 
                                      total_participant=total_participant[sb])
  length(participant_contiVal)
  list_conti_par_total[[sb]]<- participant_contiVal
}

# Get variables which appear across all of the participants in the list. 
colName.x.continuous.participant<- list_conti_par_total[[1]]
for( j in 2:length(list_conti_par_total)){
  colName.y.continuous.participant<- list_conti_par_total[[j]]
  colName.x.continuous.participant<- intersect(colName.x.continuous.participant, 
                                               colName.y.continuous.participant)
}

# Combine discrete and continuous variables into one. 
property_P_list<- c(colName.x.discrete.participant, 
                  colName.x.continuous.participant)

# Get the subset of the dataset called 'filt_meta' used above the 
# candidate list of properties for a Participant node.  
participant_Property<- filt_meta[,match(property_P_list, colnames(filt_meta))]

# Check if there are duplicates in the dataset. 
participant_Node_prop<- 
  participant_Property[which(duplicated(participant_Property$Participant.ID)==
                               FALSE),]

################################################################################
# Main Program: Create a Sample Node                                           #
# Input: filt_meta, continuous_val, and discrete_val                           #
# Return: This should result in a dataframe called 'sample_Node_prop'          #
#         that contains 4337 observations of 168 variables.                    #
################################################################################

# Let's identify which discrete variables can be the candidate properties
# for a Sample node. As we know the variable (property) for a Participant node, 
# we exclude these variables from a dataframe called 'discrete_val'. 
disVal_after_pNode<- discrete_val %>% 
                     dplyr::select(!(colName.x.discrete.participant))

# First, let's identify which discrete variables can be the candidate properties
# for a Sample node. 
total_sample<- unique(filt_meta$External.ID)

# Create new empty lists to hold discrete variables properties for 2505 samples. 
list_disc_samp_total<- c()

# This for loop will run until we identify the discrete variables that can be 
# the candidate properties for 2505 samples. 
for( sp in 1:length(total_sample)){
  samp_discreteVal<- get_discVal_for_SampleNode(disVal_after_pNode, 
                                                total_sample=total_sample[sp])
  length(samp_discreteVal)
  list_disc_samp_total[[sp]]<- samp_discreteVal
}

# Get variables which appear across all of the samples in the list. 
colName.x.discrete.sample<- list_disc_samp_total[[1]]
for( j in 2:length(list_disc_samp_total)){
  colName.y.discrete.sample<- list_disc_samp_total[[j]]
  colName.x.discrete.sample<- intersect(colName.x.discrete.sample, 
                                        colName.y.discrete.sample)
}

# Once this is complete, let's identify which continuous variables can be the 
# candidate properties for a Sample node. 
contiVal_after_pNode<- continuous_val %>% 
                       dplyr::select(!(colName.x.continuous.participant))

# Create new empty lists to hold discrete variables properties for 2505 samples. 
list_conti_samp_total<- c()

# This for loop will run until we identify the continuous variables that can be 
# the candidate properties for 2505 samples. 
for( sb in 1:length(total_sample)){
  sample_contiVal<- get_contiVal_for_SampleNode(contiVal_after_pNode, 
                                                total_sample=total_sample[sb])
  length(sample_contiVal)
  list_conti_samp_total[[sb]]<- sample_contiVal
}

# Get variables which appear across all of the samples in the list. 
colName.x.continuous.sample<- list_conti_samp_total[[1]]
for( j in 2:length(list_conti_samp_total)){
  colName.y.continuous.sample<- list_conti_samp_total[[j]]
  colName.x.continuous.sample<- intersect(colName.x.continuous.sample, 
                                          colName.y.continuous.sample)
}

# Combine discrete and continuous variables into one. 
property_S_list<- c(colName.x.discrete.sample, colName.x.continuous.sample)

# Get the subset of the dataset called 'filt_meta' used above the candidate list 
# of properties for a Sample node. 
sample_Node_prop<- filt_meta[,match(property_S_list, colnames(filt_meta))]

# Note: we want to add 'data_type' as a property for a Sample node. We know that 
# the 'data_type' is not the candidate property for neither a Participant nor 
# Sample node. In the metadata, each of the participants can have same multiple 
# sample Ids (External.ID) with different types of the data. 
#
# Including the data_type property for a Sample node will help users
# to query in an efficient manner. It is also useful to have for purpose of 
# assigning unique sample-ids for a sample node table. 
new_samp_Node<- filt_meta %>% 
  dplyr::select(c(colnames(sample_Node_prop),"data_type"))

################################################################################
# Main Program: Create and add 'unique_Ids' for a Participant and Sample node  #
#             This step would useful for querying onto our populated database. #
#                                                                              #
# Input: participant_Node_prop, new_sampl_Node,                                #
# Return: participant_node, sample_node                                        #
################################################################################
source(paste0(sourceDir,"/","formatting_data_neo4j.R"))
# We want to create and add 'unique_ids' to our a Participant and Sample node. 
p_uniqIds<- create_Entity_uniqueIds(entity="participant node", prefix="p", 
                                    final_propertyTab=participant_Node_prop)

s_uniqIds<- create_Entity_uniqueIds(entity="sample node", prefix="s",
                                    final_propertyTab=new_samp_Node)

# Add created above unique-Ids to each of the Nodes. 
participant_node<- cbind("p_ID" = p_uniqIds, participant_Node_prop)
sample_node<- cbind("s_ID" = s_uniqIds, new_samp_Node)

################################################################################
# Main Program: Create a Participant-Sample Edge                               #
# Input: participant_node, sample_node, and 'filt_meta'                        #
# Return: This should result in a dataframe called 'ps_Edge' that contains     #
#         4337 observations of 14 variables.                                   #
################################################################################
# Now let's create a participant-sample edge (Hereafter we will use 'p-s-E').
# To identify the candidate property for a 'p-s-E', first we will need and use 
# the unique IDs created above for each of the nodes and merge with metadata. 
# Once this is complete, we can begin to extract the candidate property for a 
# 'p-s-E' relationship. 
list_pIds<- participant_node[,c("p_ID", "Participant.ID")]
first_merge<- merge(list_pIds, filt_meta)
second_merge<- merge(first_merge, sample_node)

# First, we identify and retain the properties for a Participant Node except 
# participant unique ID created above. 
property_P<- setdiff(colnames(participant_node), c("p_ID"))

# Next up, we identify and retain the properties for a Sample Node except sample
# unique ID and data_type created above. 
property_S<- setdiff(colnames(sample_node), c("s_ID", "data_type"))

# Extract the properties for a Participant and Sample node from the created 
# above merged data. 
ps_Edge<- second_merge %>% 
          dplyr::select(-property_P) %>% dplyr::select(-property_S)

################################################################################
# Main Program: Reformat data in neo4j-admin import format                     #
# Input: participant_node, sample_node, and 'ps_Edge'                          #
# Return:                                                                      #
# Note: Required source file called "formatting_data_neo4j.R"                  #
################################################################################
# Process Participant Node data in neo4j-import format 
edited_P_Node<- reformat_Node_file(participant_node, old_col="p_ID", 
                                   new_col="p_ID:ID", 
                                   node_label="Participant")

# Process Sample Node data in neo4j-import format  
edited_S_Node<- reformat_Node_file(sample_node, 
                                   old_col="s_ID", 
                                   new_col="s_ID:ID", 
                                   node_label="Sample")

# Process p-s-E in neo4j-import format 
edited_ps_Edge<- reformat_Edge_file(ps_Edge, start_col="p_ID", end_col="s_ID",
                                    edge_label="HAS_SAMPLE_AT_TIME_POINTS")

# Remove special letters and characters. 
neo4j_participantNode<- edit_specialCharacters(edited_P_Node)
neo4j_sampleNode<- edit_specialCharacters(edited_S_Node)
neo4j_psEdge<- edit_specialCharacters(edited_ps_Edge)

# Note: The next few steps of type conversion by adding suffixing the name with 
# property value. At the present, we manually checked each of the variables to 
# assign a property value. 

# Participant Node Property 
print("Rename columns for participant node dataframe")

names(neo4j_participantNode)[names(neo4j_participantNode) == 
                                   "Age_at_diagnosis"]<- "Age_at_diagnosis:long"
names(neo4j_participantNode)[names(neo4j_participantNode) == 
                           "Other_inflamed_flora"]<- "Other_inflamed_flora:long"
names(neo4j_participantNode)[names(neo4j_participantNode) == 
                               "Rectum_cell_biopsy"]<- "Rectum_cell_biopsy:long"
names(neo4j_participantNode)[names(neo4j_participantNode) == 
                                 "Ileum_cell_biopsy"]<- "Ileum_cell_biopsy:long"
names(neo4j_participantNode)[names(neo4j_participantNode) == 
       "Other_non_inflamed_cell_biopsy"]<- "Other_non_inflamed_cell_biopsy:long"
names(neo4j_participantNode)[names(neo4j_participantNode) == 
                                             "consent_age"]<- "consent_age:long"
names(neo4j_participantNode)[names(neo4j_participantNode) == 
          "If_yes_how_many_times_per_month_do_you_smoke_marijuana_on_average"]<- 
        "If_yes_how_many_times_per_month_do_you_smoke_marijuana_on_average:long"
names(neo4j_participantNode)[names(neo4j_participantNode) == 
                       "If_no_in_what_year_did_you_come_to_the_United_States"]<- 
                     "If_no_in_what_year_did_you_come_to_the_United_States:long"
names(neo4j_participantNode)[names(neo4j_participantNode) == 
                         "If_yes_how_many_times"]<- "If_yes_how_many_times:long"
names(neo4j_participantNode)[names(neo4j_participantNode) == "BMI"]<- 
                                                                    "BMI:double"
names(neo4j_participantNode)[names(neo4j_participantNode) == "Height"]<- 
                                                                 "Height:double"
names(neo4j_participantNode)[names(neo4j_participantNode) == "Weight_1"]<- 
                                                               "Weight_1:double"
names(neo4j_participantNode)[names(neo4j_participantNode) == "SIBDQ_Score"]<- 
                                                              "SIBDQ_Score:long"

# Sample Node Property  
print("Rename columns for sample node dataframe")

names(neo4j_sampleNode)[names(neo4j_sampleNode) == "week_num"]<- "week_num:long"
names(neo4j_sampleNode)[names(neo4j_sampleNode) == "interval_days"]<-
                                                            "interval_days:long"
names(neo4j_sampleNode)[names(neo4j_sampleNode) == "reads_ribosomal"]<- 
                                                          "reads_ribosomal:long"
names(neo4j_sampleNode)[names(neo4j_sampleNode) == "reads_viral"]<- 
                                                              "reads_viral:long"
names(neo4j_sampleNode)[names(neo4j_sampleNode) == 
     "Number_of_flora_tubes_collected"]<- "Number_of_flora_tubes_collected:long"
names(neo4j_sampleNode)[names(neo4j_sampleNode) == "Rectum_Flora"]<- 
                                                             "Rectum_Flora:long"
names(neo4j_sampleNode)[names(neo4j_sampleNode) == "Ileum_flora"]<- 
                                                              "Ileum_flora:long"
names(neo4j_sampleNode)[names(neo4j_sampleNode) == "Non_inflamed_flora"]<- 
                                                       "Non_inflamed_flora:long"
names(neo4j_sampleNode)[names(neo4j_sampleNode) == 
                    "Number_of_tubes_collected_for_epithelial_cell_biopsies"]<-
                   "Number_of_tubes_collected_for_epithelial_cell_biopsies:long"
names(neo4j_sampleNode)[names(neo4j_sampleNode) == 
                                         "Number_of_DNA_RNA_tubes_collected"]<- 
                                        "Number_of_DNA_RNA_tubes_collected:long"
names(neo4j_sampleNode)[names(neo4j_sampleNode) == "Weight"]<- "Weight:double"
names(neo4j_sampleNode)[names(neo4j_sampleNode) == "fecalcal_ng_ml"]<-
                                                         "fecalcal_ng_ml:double"
names(neo4j_sampleNode)[names(neo4j_sampleNode) == 
                 "Number_of_liquid_or_very_soft_stools_in_the_past_24_hours"]<- 
              "Number_of_liquid_or_very_soft_stools_in_the_past_24_hours:double"
names(neo4j_sampleNode)[names(neo4j_sampleNode) == "hbi"]<- "hbi:double"
names(neo4j_sampleNode)[names(neo4j_sampleNode) == "CRP_mg_L"]<- 
                                                               "CRP_mg_L:double"
names(neo4j_sampleNode)[names(neo4j_sampleNode) == "ESR_mm_hr"]<- 
                                                              "ESR_mm_hr:double"
names(neo4j_sampleNode)[names(neo4j_sampleNode) == "Modified_Baron_s_Score"]<- 
                                                   "Modified_Baron_s_Score:long"
names(neo4j_sampleNode)[names(neo4j_sampleNode) == 
                                      "Highest_dose_of_Prednisone_taken_mg"]<- 
                                      "Highest_dose_of_Prednisone_taken_mg:long"
names(neo4j_sampleNode)[names(neo4j_sampleNode) == 
                                        "Highest_dose_of_Entocort_taken_mg"]<-
                                        "Highest_dose_of_Entocort_taken_mg:long"
names(neo4j_sampleNode)[names(neo4j_sampleNode) == 
         "Duration_of_Humira_use_months"]<- "Duration_of_Humira_use_months:long"
names(neo4j_sampleNode)[names(neo4j_sampleNode) == 
         "Duration_of_Cimzia_use_months"]<- "Duration_of_Cimzia_use_months:long"
names(neo4j_sampleNode)[names(neo4j_sampleNode) == "sccai"]<- "sccai:long"
names(neo4j_sampleNode)[names(neo4j_sampleNode) == "Total"]<- "Total:long"
names(neo4j_sampleNode)[names(neo4j_sampleNode) == "Total_1"]<- "Total_1:long"
names(neo4j_sampleNode)[names(neo4j_sampleNode) == "Total_2"]<- "Total_2:long"
names(neo4j_sampleNode)[names(neo4j_sampleNode) == "Total_3"]<- "Total_3:long"
names(neo4j_sampleNode)[names(neo4j_sampleNode) == "SES_CD_Score"]<- 
                                                            "SES_CD_Score:long"
names(neo4j_sampleNode)[names(neo4j_sampleNode) == "fecalcal"]<- 
                                                              "fecalcal:double"
names(neo4j_sampleNode)[names(neo4j_sampleNode) == "GSSR_IDs"]<- "GSSR_IDs:long"
names(neo4j_sampleNode)[names(neo4j_sampleNode) == "WR_ID"]<- "WR_ID:long"
names(neo4j_sampleNode)[names(neo4j_sampleNode) == "Lanes_in_Aggregation"]<- 
                                                     "Lanes_in_Aggregation:long"

# Participant-Sample Edge Property 
print("Rename columns for participant-sample edge/relationship dataframe")
names(neo4j_psEdge)[names(neo4j_psEdge) == "reads_raw"]<- "reads_raw:long"
names(neo4j_psEdge)[names(neo4j_psEdge) == "reads_filtered"]<- 
                                                           "reads_filtered:long"
names(neo4j_psEdge)[names(neo4j_psEdge) == "reads_qc_fail"]<- 
                                                            "reads_qc_fail:long"
names(neo4j_psEdge)[names(neo4j_psEdge) == "reads_human"]<- "reads_human:long"


# Note: we want to decode columns 113 and 114 for more meaningful values.(How to break this code before the line)
neo4j_sampleNode$Bowel_frequency_during_the_day[neo4j_sampleNode$Bowel_frequency_during_the_day == "6-Apr"]<- "4 - 6"
neo4j_sampleNode$Bowel_frequency_during_the_day[neo4j_sampleNode$Bowel_frequency_during_the_day == "9-Jul"]<- "7 - 9"
neo4j_sampleNode$Bowel_frequency_during_the_night[neo4j_sampleNode$Bowel_frequency_during_the_night == "3-Jan"]<- "1 - 3"
neo4j_sampleNode$Bowel_frequency_during_the_night[neo4j_sampleNode$Bowel_frequency_during_the_night == "6-Apr"]<- "4 - 6"

# We need to substitute NA with empty for neo4j to recognize NA as a string 
# data type 
neo4j_participantNode[is.na(neo4j_participantNode)]<-""
neo4j_sampleNode[is.na(neo4j_sampleNode)]<- ""
neo4j_psEdge[is.na(neo4j_psEdge)]<- ""

# Write the Participant node Output files to a .csv file 
fileN= paste("participant", "node", outName, "csv", sep=".")
write.csv(neo4j_participantNode, fileN, row.names=F, quote=T)

# Write the Sample node Output files to a .csv file
fileN= paste("sample", "node", outName, "csv", sep=".")
write.csv(neo4j_sampleNode, fileN, row.names=F, quote=T)

# Write the Participant-Sample edge Output files to a .csv file
fileN= paste("participant_sample", "edge", outName, "csv", sep=".")
write.csv(neo4j_psEdge, fileN, row.names=F, quote=T)
