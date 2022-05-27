################################################################################  
# Script: gene.Node_AND_Edge.new.R                                             #
# Author: Suyeon Kim                                                           #
# Created: April 13, 2022                                                      #
# Last edited: April 13, 2022                                                  #
#                                                                              #
# Goal: The purpose of this script is to create a Gene node and sample-hostgene# 
#       edges file and sample-microbialgene edges in neo4j bulk importer format#
# Input: 1) Three processed Gene data (metatranscriptomic and host             #
#           metranscriptome folder from ibdmdb database are needed.            #
#        2) sample.node.[output argument].csv                                  #
# Output: a GENE node and sample_microbialGene/ hostGene edge list             #                         
################################################################################
library(stringr)
library(dplyr)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)
print(args)

sourceDir <- args[1]
InFile1<- args[2]
InFile2<- args[3]
sample_node<- args[4]
outName<- args[5]

# This source file is to check the global environment to be 
# "stringsAsFactors=FALSE".
source(paste0(sourceDir,"/","getOption_setOption.R"))

# In this function, we want to replace the gene name and sample ID with 
# the created sample Ids and gene Ids. 
get_Ids<- function(abundance, entity_Node, entity_1, entity_2, sample_info){
  print(dim(abundance))
  # Merge abudance data with Gene Node dataframe based on entity_1 
  # entity_1 in microbial gene abundance is "UniRef90" 
  # entity_1 in host gene abundance is "Gene"
  merged_GeneId<- merge(abundance, entity_Node, by=c(entity_1))
  filter_cols<- merged_GeneId %>% select(-c(entity_1, entity_2))
  matched_sampleIndex<- 
    sample_info[match(colnames(filter_cols[,-ncol(filter_cols)]), 
                      sample_info$External_ID),1]
  colnames(filter_cols)[1:ncol(filter_cols)-1]<- matched_sampleIndex
  # Create abundance table only with created gene-ids and sample-ids
  final_ab<- cbind("gene_ID"=filter_cols$gene_ID, 
                   filter_cols[,-ncol(filter_cols)])
  print("**********Before Reshape final_ab *******")
  print(dim(final_ab))
  
  return(final_ab)
}

# In this function, we convert dataframe from wide to long format.
reshape_dat<- function(new_abundance, by_melt){
  reshape_dat<- melt(new_abundance, id.vars=c(by_melt))
  final_ab<- data.frame(cbind("s_ID:ID"=as.character(reshape_dat$variable), 
                              "gene_ID:ID"= reshape_dat$gene_ID, 
                              "Abundance:double"=reshape_dat$value), 
                        check.names=FALSE)
  print("**********After reshape final_ab and entity_Node dimensions *******")
  print(dim(final_ab))
 
  return(final_ab)
}
################################################################################
# Read the Input file into R.                                                  #
################################################################################
# This should result in a dataframe called 'hGene_dat' that contains 
# 55765 observations of 255 variables. 
hGene_dat<- read.csv(InFile1, header=TRUE, sep="\t", check.names=FALSE)
print("The size of the host gene abundance")
print(dim(hGene_dat))

# The TSV file downloaded was 3.7 GB. The file itself contained 2,164,740 lines 
# and is tab-delimited. 
mGene_dat<- read.csv(InFile2, header=TRUE, comment.char = "!", check.names=FALSE, 
                     sep="\t")

# Lets read the 'sample node' dataframe.
sample_N<- read.csv(sample_node, header=T, sep=",",check.names=F)
sample_info<- sample_N %>% select(c("s_ID:ID", "External_ID", "data_type"))

################################################################################
# Main Program: Create a GENE Node                                             #
# Input: hGene_dat and mGene_dat                                               #
# Return: gene_node                                                            #
# Note: Required source files - formatting_data_neo4j.R, filtering_Input.R     #
################################################################################
source(paste0(sourceDir,"/","formatting_data_neo4j.R"))
source(paste0(sourceDir,"/","filtering_Input.R"))

# The dataframe that has been filtered by 'filt_Input' function. 
filt_mGene<- filt_Input(mGene_dat, entity="UniRef90")
# Rename for sample names in columns. Let's Keep technical recpliation(TR) 
# samples 
colnames(filt_mGene)[2:ncol(filt_mGene)]<- gsub("_Abundance-RPKs.*", "", 
                                      colnames(filt_mGene)[2:ncol(filt_mGene)])

filt_hGene<- filt_Input(hGene_dat, entity="Gene")

# First, let's combine each column containing gene-relevant information from the 
# all inputs. 
bac_gene<- data.frame(filt_mGene$UniRef90)
colnames(bac_gene)<- "UniRef90"
bac_gene$Gene<- ""

host_gene<- data.frame(filt_hGene$Gene)
colnames(host_gene)<- "Gene"
host_gene$UniRef90<- ""

col_order<- c("UniRef90", "Gene")
new_hostGene<- host_gene[,col_order]

# Merge 'bac_gene' and 'new_hostGene' dataframes into a dataframe called 
# 'all_gene'. 
all_gene<- rbind(bac_gene, new_hostGene)

# We want to create and add 'unique_ids' to a Gene Node. 
gene_uniqIds<- create_Entity_uniqueIds(entity="Gene node", prefix="gene", 
                                       final_propertyTab=all_gene)

# Add created above unique-Ids to the Gene Node. 
gene_node<- cbind("gene_ID"=gene_uniqIds, all_gene)

################################################################################
# Main Program: Create a Sample-Microbial-Gene Edge                            #
# Input: gene_node, sample_info, and filt_mGene                                #
# Return: reshape_mGene dataframe                                              #
################################################################################
# Subset MTX samples in sample_info dataframe.
mtx_sample<- sample_info[sample_info$data_type == "metatranscriptomics",]
# Let's get Gene-Ids and Sample-Ids onto the abundance dataframe. 
new_mGene_ab<- get_Ids(filt_mGene, entity_Node = gene_node, 
                       entity_1 = "UniRef90", entity_2 = "Gene", 
                       sample_info = mtx_sample)

# Once this complete, we can change the dataframe from long format to wide 
# format
reshape_mGene<- reshape_dat(new_mGene_ab, by_melt = "gene_ID")

################################################################################
# Main Program: Create a Sample-Host-Gene Edge.                                #
# Input: gene_node, sample_info, and filt_hGene                                #
# Return: reshape_hGene dataframe                                              #
################################################################################
# Subset MTX samples in sample_info dataframe.
htx_sample<- sample_info[sample_info$data_type == "host_transcriptomics",]
# In human transcriptome data, we identified missing information for few samples. 
# Since we cannot find their subject Id (participant Id) in metadata, these 
# samples will be removed from our input data.
new_filt_hGene<- filt_hGene %>% select(-c("CSMDRVY8","CSMDRVY9","HSm9JTBX"))
# Let's get Gene-Ids and Sample-Ids onto the abundance dataframe. 
new_hGene_ab<- get_Ids(new_filt_hGene, entity_Node = gene_node, 
                       entity_1 = "Gene", entity_2 = "UniRef90", 
                       sample_info = htx_sample)

# Once this complete, we can change the dataframe from long format to wide 
# format
reshape_hGene<- reshape_dat(new_hGene_ab, by_melt = "gene_ID")

################################################################################
# Main Program: Reformat data in neo4j-admin import format and Save Output     #
# Input: gene_node, reshape_mGene, and reshape_hGene                           #
# Return: neo4j_gene_Node, neo4j_mGene_edge, and neo4j_hGene_edge              #
# Note: Required source file called "formatting_data_neo4j.R"                  #
################################################################################
# Process Gene_node in neo4j-import format 
neo4j_gene_Node<- reformat_Node_file(gene_node, old_col="gene_ID", 
                                     new_col="gene_ID:ID", 
                                     node_label="Gene")

# Process Sample-Microbial GENE edge/relationship in neo4j-import format 
neo4j_mGene_Edge<- reformat_Edge_file(reshape_mGene, start_col="s_ID:ID", 
                                      end_col="gene_ID:ID", 
                                      edge_label="HAS_MICROBIAL_GENE_ABUNDANCE")

# The size of dataframe called 'edited_mGene_Edge' contains a lot of information
# , we want to remove rows if the values in the column('Abundance') is >0.
filt_neo4j_mGene_Edge<- filter(neo4j_mGene_Edge, `Abundance:double`>0)

# Process Sample-Host GENE edge/relationship in neo4j-import format 
neo4j_hGene_Edge<- reformat_Edge_file(reshape_hGene, start_col="s_ID:ID", 
                                      end_col="gene_ID:ID", 
                                      edge_label="HAS_HOST_GENE_ABUNDANCE")

# We need to substitute NA with empty as Neo4 recognize NA as a string data type 
neo4j_gene_Node[is.na(neo4j_gene_Node)]<-""
filt_neo4j_mGene_Edge[is.na(filt_neo4j_mGene_Edge)]<- ""
neo4j_hGene_Edge[is.na(neo4j_hGene_Edge)]<- ""

# Write the Gene node Output file to a .csv file 
fileN1= paste("Gene", "node", outName, "csv", sep=".")
write.csv(neo4j_gene_Node, fileN1, row.names=F, quote=T, sep="\t")

# Write the sample-Microbial Gene file to a .csv file 
fileN2= paste("sample_microbialGene","edge", outName, "csv", sep=".")
write.csv(filt_neo4j_mGene_Edge, fileN2, row.names=F, quote=F, sep=",")

# Write the sample-Host Gene file to a .csv file 
fileN3= paste("sample_hostGene","edge", outName, "csv", sep=".")
write.csv(neo4j_hGene_Edge, fileN3, row.names=F, quote=F, sep=",")
