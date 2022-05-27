################################################################################
# Script: taxon.Node_AND_Edge.R                                                #
# Author: Suyeon Kim                                                           #
# Created: April 07, 2022                                                      #
# Last edited: April 11, 2022                                                  # 
#                                                                              #
# Goal: The purpose of this script is to create a nodes and edges file in      #
# neo4j bulk importer format (neo4j-admin import tool).                        #
# Input: three 'taxonomic_profiles.tsv' from three sources: MGX/2018-05-04/,   #
#        16S/2018-01-07/, and MVX/                                             #  
#                                                                              #
# Output:                                                                      #
# (in neo4j bulk importer format).                                             #
################################################################################
library(stringr)
library(dplyr)
library(reshape2)

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
print(args)
sourceDir <- args[1]
InFile1<- args[2]
InFile2<- args[3]
InFile3<- args[4]
sample_node<- args[5]
outName<- args[6]

# This source file is to check the global environment to be 
# "stringsAsFactors=FALSE".
source(paste0(sourceDir,"/","getOption_setOption.R"))

create_taxonLineage<- function(rank, list_taxonomy){
  # Let's identify the lineage for each taxonomy rank
  print(paste("This is the taxonomy rank", rank))
  parent_Info<- unlist(lapply(list_taxonomy, 
                function(x) {paste(x[1:as.numeric(rank-1)], collapse=";")}))
  child_Info<- unlist(lapply(list_taxonomy, function(x) {x[rank]}))
  parent_child<- data.frame(cbind(parent_Info, child_Info))
  # we will check if there are NAs, and remove NAs. 
  new_parent_child<- parent_child[complete.cases(parent_child),]
  # Next, we remove duplicates by removing rows for which child_Info is 
  # duplicates. 
  unique_parent_child<- unique(new_parent_child)
  unique_parent_child$`lineage`<- paste(unique_parent_child[,1], 
                                        unique_parent_child[,2], sep="|")
  return(unique_parent_child)
}

format_TaxonNode<- function(rank, taxon_Inf, taxon_name, 
                            prev_Taxon) {
  # In order to assigning taxon-Id index, we need multiple resources to create 
  # index. We need to know current 'list_taxonomy' number of rows along with 
  # previous taxon's very last taxon-Id. 
  print(paste("This is the taxonomy rank", rank))
  print(dim(taxon_Inf))
  if(rank == 1) {
    start_num<-1
    end_num<- nrow(taxon_Inf)
  }
  else{
    last_taxId<- prev_Taxon[nrow(prev_Taxon),"tx_Id:ID"]
    raw_num<- gsub("tx","", last_taxId)
    prev_num<- as.numeric(raw_num)
    start_num<-prev_num + 1 
    end_num<- start_num + (nrow(taxon_Inf)-1)  
  }
  # Once we know start number and end number for a 'list_taxonomy' dataframe in 
  # current rank, we can create taxon_Ids. 
  prefix<- "tx" 
  suffix<- seq(start_num, end_num)
  taxon_Id<- paste(prefix, suffix, sep="")
  
  # Reorder columns the following orders. Make sure taxon_Id is placed in first
  # place. 
  print(taxon_name)
  taxon_Node<- data.frame(cbind(taxon_Id, taxon_Inf$child_Info, 
                  taxon_Inf$parent_Info, taxon_Inf$lineage))
  print(dim(taxon_Node))
  print(class(taxon_Node))
  colnames(taxon_Node)<-c("tx_Id:ID",taxon_name,"ancestral_lineage","lineage")
  taxon_Node$ancestral_lineage<-gsub("\\;","\\|", taxon_Node$ancestral_lineage)
  taxon_Node$lineage<- gsub("\\;","\\|",taxon_Node$lineage)
  # We also want to add the taxonomy level information. 
  taxon_Node$level<- paste0("level",rank)
  # Reformat dataframe based on the Neo4j format
  taxon_Node$`:LABEL`<- taxon_name
  
  # Write the Taxon node Output file to a .csv file 
  #fileN= paste(taxon_name, "node", outName, "csv", sep=".")
  #write.csv(taxon_Node, fileN, row.names=F, quote=T, sep="\t")
  
  return(taxon_Node)
}

get_abundance<- function(mgx_tax, rank){
  # We want to create abundance dataframes based on the taxonomic rank.
  taxon_dat<- mgx_tax[which(str_count(rownames(mgx_tax), "\\|") == rank),]
  tax_ab<- data.frame("lineage"=rownames(taxon_dat), taxon_dat)
  rownames(tax_ab)<- NULL
  return(tax_ab)
}

get_16s_Ids<- function(taxon_dat, entity_Node, entity, sample_info){
  print(dim(taxon_dat))
  merged_Taxon<-merge(taxon_dat, entity_Node, by.x=c(entity,'ancestral_lineage'), 
                      by.y=c(entity,'ancestral_lineage'))
  filter_cols<- merged_Taxon %>% 
    select(-c(entity,"ancestral_lineage","lineage","level",":LABEL"))
  matched_sampIndex<-
    sample_info[match(colnames(filter_cols[,-ncol(filter_cols)]), 
                      sample_info$External_ID),1]
  colnames(filter_cols)[1:ncol(filter_cols)-1]<- matched_sampIndex
  final_ab<- cbind("tx_Id:ID"=filter_cols$`tx_Id:ID`, 
                   filter_cols[,-ncol(filter_cols)]) 
  print(dim(final_ab))
  return(final_ab)
}

get_mgx_Ids<- function(taxon_dat, entity_Node, entity, sample_info){
  print(dim(taxon_dat))
  merged_mgx<- merge(taxon_dat, entity_Node, by=c("lineage"))
  filter_cols<- merged_mgx %>% 
    select(-c("lineage",entity,"ancestral_lineage","level",":LABEL"))
  matched_sampIndex<-
    sample_info[match(colnames(filter_cols[,-ncol(filter_cols)]), 
                      sample_info$External_ID),1]
  colnames(filter_cols)[1:ncol(filter_cols)-1]<- matched_sampIndex
  final_ab<- cbind("tx_Id:ID"=filter_cols$`tx_Id:ID`, 
                   filter_cols[,-ncol(filter_cols)]) 
  print(dim(final_ab))
  return(final_ab)
}

reshape_taxonDat<- function(abundance_dat, by_melt){
  reshape_dat<- melt(abundance_dat, id.vars=c(by_melt))
  final_ab<- data.frame(cbind("s_ID"=as.character(reshape_dat$variable),
                              "tx_Id:ID"=as.character(reshape_dat$`tx_Id:ID`), 
                              "Abundance:double"=as.numeric(reshape_dat$value)), 
                        check.names = F)
  return(final_ab)
}

edge_Neo4j<- function(reshape_dat, edge_label){
  names(reshape_dat)[names(reshape_dat) == "s_ID"]<- ":START_ID"
  names(reshape_dat)[names(reshape_dat) == "tx_Id:ID"]<- ":END_ID"
  colnames(reshape_dat)[3]<- "Abundance:double"
  reshape_dat$`:TYPE`<- edge_label
  
  return(reshape_dat)
}

################################################################################
# Read the Input file into R.                                                  #
################################################################################
# This should result in a dataframe called 'biopsy_tax' that contains 
# 982 observations of 179 variables. 
biopsy_tax<- read.csv(InFile1, header=TRUE, comment.char = "$", sep="\t", 
                      row.names=1,stringsAsFactors = FALSE, check.names=FALSE)
print(dim(biopsy_tax))
# This should result in a dataframe called 'mgx_tax' that contains 
# 1479 observations of 1638 variables. 
mgx_tax<- read.csv(InFile2, header=TRUE, comment.char = "$", sep="\t", 
                   row.names=1,stringsAsFactors = FALSE, check.names=FALSE)
print(dim(mgx_tax))
# This should result in a dataframe called 'mvx_tax' that contains 
# 56 observations of 329 variables. 
mvx_tax<- read.csv(InFile3, header=TRUE, comment.char = "$", sep="\t", 
                   row.names=1,stringsAsFactors = FALSE, check.names=FALSE)
print(dim(mvx_tax))
# Lets read the 'sample node' file into R  
sample_N<- read.csv(sample_node, header=T, sep=",", check.names=FALSE)
sample_info<- sample_N %>% select(c("s_ID:ID", "External_ID", "data_type"))
print(dim(sample_info))
################################################################################
# Main Program: Create "Kingdom, Phylum, Class, Order, Family, Genus, Species  #
#               Strain" Nodes                                                  #
# Input: biopsy_tax, mgx_tax, mvx_tax                                          #
# Return:                                                                      #
################################################################################
#
### Data pre-processing
### As each abundance data were generated in different sequencing technologies,
### Taxonomy columns in created previous dataframes need to be generalized. 
#

# Subset taxonomy rownames from 'mgx_tax' dataframe
taxonomy_MGX<- data.frame("taxonomy"=rownames(mgx_tax))
# Subset taxonomy rownames from 'mvx_tax' dataframe
taxonomy_MVX<- data.frame("taxonomy"=rownames(mvx_tax))
# Subset taxonomy id with taxonomy column from 'biopsy_tax' dataframe
taxonomy_16s<- data.frame(cbind("otu_ID"=rownames(biopsy_tax), 
                                "taxonomy"=biopsy_tax$taxonomy))
# Split 'taxonomy' column into number of taxons in taxonomy lineage.
split_taxon<- data.frame(str_split_fixed(taxonomy_16s$taxonomy, ";", n=6))
colnames(split_taxon)<- c("Kingdom","Phylum","Class","Order","Family","Genus")
# We want to format of 16s taxonomy name to MGX taxonomy name to make 
# comparability of taxonomy. 
split_taxon$Genus<- gsub("\\___","\\__", split_taxon$Genus)
split_taxon$Genus<- gsub("__","g__", split_taxon$Genus)
split_taxon$Family<- gsub("__","f__", split_taxon$Family)
split_taxon$Order<- gsub("__","o__", split_taxon$Order)
split_taxon$Class<- gsub("__","c__", split_taxon$Class)
split_taxon$Phylum<- gsub("__","p__", split_taxon$Phylum)
split_taxon$Kingdom<- gsub("Bacteria", "k__Bacteria", split_taxon$Kingdom)
split_taxon$Kingdom<- gsub("Archaea", "k__Archaea", split_taxon$Kingdom)

new_16sTaxonomy<- paste0(split_taxon$Kingdom, "|", split_taxon$Phylum, "|", 
                         split_taxon$Class, "|", split_taxon$Order,
                         "|" ,split_taxon$Family,"|", split_taxon$Genus)
new_16sTaxonomy<- gsub("\\| ","|", new_16sTaxonomy)
new_16sTaxonomy<- data.frame("taxonomy"=(new_16sTaxonomy))

# Unique taxonomy information from 16s biopsy
unique_16sTax<- unique(new_16sTaxonomy)

# Once this is complete, we can begin to combine taxonomy information from three
# different sources. 
print("Union of all taxonomy")
union_Taxonomy<- union(union(unique_16sTax,taxonomy_MGX),taxonomy_MVX)
print(dim(union_Taxonomy))
list_taxonomy<- strsplit(as.character(union_Taxonomy$taxonomy), split="\\|")

kingdom_Inf<- create_taxonLineage(rank=1, list_taxonomy)
print(dim(kingdom_Inf))
kingdom_Node<- format_TaxonNode(rank=1, kingdom_Inf,
                                taxon_name = "Kingdom",prev_Taxon=0)
kingdom_Node$lineage<- gsub("\\|.*", "", kingdom_Node$lineage)
print(paste("This is the size of ", dim(kingdom_Node)))

# Create a Phylum Node 
phylum_Inf<- create_taxonLineage(rank=2, list_taxonomy)
phylum_Node<- format_TaxonNode(rank=2, phylum_Inf, 
                               taxon_name = "Phylum", prev_Taxon=kingdom_Node)
print(paste("This is the size of ", dim(phylum_Node)))

# Create a Class Node 
class_Inf<- create_taxonLineage(rank=3, list_taxonomy)
class_Node<- format_TaxonNode(rank=3, class_Inf, 
                              taxon_name = "Class", prev_Taxon=phylum_Node)
print(paste("This is the size of ", dim(class_Node)))
# Create a Order Node 
order_Inf<- create_taxonLineage(rank=4, list_taxonomy)
order_Node<- format_TaxonNode(rank=4, order_Inf, 
                              taxon_name = "Order", prev_Taxon=class_Node)
print(paste("This is the size of ", dim(order_Node)))
# Create a Family Node 
family_Inf<- create_taxonLineage(rank=5, list_taxonomy)
family_Node<- format_TaxonNode(rank=5, family_Inf, 
                               taxon_name = "Family", prev_Taxon=order_Node)
print(paste("This is the size of ", dim(family_Node)))
# Create a Genus Node 
genus_Inf<- create_taxonLineage(rank=6, list_taxonomy)
genus_Node<- format_TaxonNode(rank=6, genus_Inf, 
                              taxon_name = "Genus", prev_Taxon = family_Node)
print(paste("This is the size of ", dim(genus_Node)))
# Create a Species Node
species_Inf<- create_taxonLineage(rank=7, list_taxonomy)
species_Node<- format_TaxonNode(rank=7, species_Inf,
                                taxon_name = "Species", prev_Taxon=genus_Node)
print(paste("This is the size of ", dim(species_Node)))
# Create a Strain Node
strain_Inf<- create_taxonLineage(rank=8, list_taxonomy)
strain_Node<- format_TaxonNode(rank=8, strain_Inf, 
                               taxon_name = "Strain", prev_Taxon=species_Node)
print(paste("This is the size of ", dim(strain_Node)))

################################################################################
# Main Program: Create a Sample-16sBiopsy-Taxon Edge                           #
# Input: Kingdom, Phylum, Class, Order, Family, Genus, Species, Strain Nodes,  #
#        sample_info, and biopsy_tax                                           #
# Return:                                                                      #
################################################################################
source(paste0(sourceDir,"/","taxonAbundance_generator_Aim1.R"))

# Drop the taxonomy columns to aggregate the columns based on the taxonomic rank
biopsy_ab<- biopsy_tax %>% select(-"taxonomy")
# We want to aggreagte the table based on the taxonomic rank. 
genus_dat<- aggregate_dat(biopsy_ab, split_taxon, columnName = "Genus")
family_dat<- aggregate_dat(biopsy_ab, split_taxon, columnName = "Family")
order_dat<- aggregate_dat(biopsy_ab, split_taxon, columnName = "Order")
class_dat<- aggregate_dat(biopsy_ab, split_taxon, columnName = "Class")
phylum_dat<- aggregate_dat(biopsy_ab, split_taxon, columnName = "Phylum")
kingdom_dat<- aggregate_dat(biopsy_ab, split_taxon, columnName = "Kingdom")

#save_output(genera_dat, "Genera", "16sBiopsy")
#save_output(family_dat, "Family", "16sBiopsy")
#save_output(order_dat, "Order", "16sBiopsy")
#save_output(class_dat, "Class", "16sBiopsy")
#save_output(phylum_dat, "Phylum", "16sBiopsy")
#save_output(kingdom_dat, "Kingdom", "16sBiopsy")

# Replace sample-ids and taxonomy info with our created sample-ids and tax-ids 
# in Taxon node. [Better Codelines]
# Subset 16sbiopsy samples in sample_info dataframe. 
biopsy_sample<- sample_info[sample_info$data_type == "biopsy_16S",]
genus_16sAb<- get_16s_Ids(genus_dat, genus_Node, entity="Genus", 
                          sample_info = biopsy_sample)
family_16sAb<- get_16s_Ids(family_dat, family_Node, entity ="Family", 
                           sample_info = biopsy_sample)
order_16sAb<- get_16s_Ids(order_dat, order_Node, entity="Order", 
                          sample_info = biopsy_sample)
class_16sAb<- get_16s_Ids(class_dat, class_Node, entity="Class", 
                          sample_info = biopsy_sample)
phylum_16sAb<- get_16s_Ids(phylum_dat, phylum_Node, entity="Phylum", 
                           sample_info = biopsy_sample)
kingdom_16sAb<- get_16s_Ids(kingdom_dat, kingdom_Node, entity="Kingdom", 
                            sample_info = biopsy_sample)
                        
# Reshape the previous created dataframes from wide to long foramt. 
reshape_genus<- reshape_taxonDat(genus_16sAb, by_melt="tx_Id:ID")
reshape_family<- reshape_taxonDat(family_16sAb, by_melt="tx_Id:ID")
reshape_order<- reshape_taxonDat(order_16sAb, by_melt="tx_Id:ID")
reshape_class<- reshape_taxonDat(class_16sAb, by_melt="tx_Id:ID")
reshape_phylum<- reshape_taxonDat(phylum_16sAb, by_melt="tx_Id:ID")
reshape_kingdom<- reshape_taxonDat(kingdom_16sAb, by_melt="tx_Id:ID")

# Once this complete, let's reformat the dataframe in Neo4j format. 
neo4j_G_edge<- edge_Neo4j(reshape_genus,
                          edge_label="HAS_16s_BIOPSY_G_ABUNDANCE")
neo4j_F_edge<- edge_Neo4j(reshape_family,
                          edge_label="HAS_16s_BIOPSY_F_ABUNDANCE")
neo4j_O_edge<- edge_Neo4j(reshape_order,
                          edge_label="HAS_16s_BIOPSY_O_ABUNDANCE")
neo4j_C_edge<- edge_Neo4j(reshape_class,
                          edge_label="HAS_16s_BIOPSY_C_ABUNDANCE")
neo4j_P_edge<- edge_Neo4j(reshape_phylum,
                          edge_label="HAS_16s_BIOPSY_P_ABUNDANCE")
neo4j_K_edge<- edge_Neo4j(reshape_kingdom,
                          edge_label="HAS_16s_BIOPSY_K_ABUNDANCE")

print("Compeleted sample-16sBiopsy data")
################################################################################
# Main Program: Create a Sample-Metagenomics-Taxon Edge                        #
# Input: Kingdom, Phylum, Class, Order, Family, Genus, Species, Strain Nodes,  #
#        sample_info, and mgx_tax                                              #
# Return:                                                                      #
################################################################################
k_mgx_ab<- get_abundance(mgx_tax, rank=0)
p_mgx_ab<- get_abundance(mgx_tax, rank=1)
c_mgx_ab<- get_abundance(mgx_tax, rank=2)
o_mgx_ab<- get_abundance(mgx_tax, rank=3)
f_mgx_ab<- get_abundance(mgx_tax, rank=4)
g_mgx_ab<- get_abundance(mgx_tax, rank=5)
sp_mgx_ab<- get_abundance(mgx_tax, rank=6)
st_mgx_ab<- get_abundance(mgx_tax, rank=7)

# Subset MGX samples in sample_info dataframe. 
mgx_sample<- sample_info[sample_info$data_type == "metagenomics",]
mgx_k<- get_mgx_Ids(k_mgx_ab, entity_Node = kingdom_Node, entity = "Kingdom", 
                    sample_info = mgx_sample)
mgx_p<- get_mgx_Ids(p_mgx_ab, entity_Node = phylum_Node, entity = "Phylum", 
                    sample_info = mgx_sample)
mgx_c<- get_mgx_Ids(c_mgx_ab, entity_Node = class_Node, entity = "Class", 
                    sample_info = mgx_sample)
mgx_o<- get_mgx_Ids(o_mgx_ab, entity_Node = order_Node, entity = "Order", 
                    sample_info = mgx_sample)
mgx_f<- get_mgx_Ids(f_mgx_ab, entity_Node = family_Node, entity = "Family", 
                    sample_info = mgx_sample)
mgx_g<- get_mgx_Ids(g_mgx_ab, entity_Node = genus_Node, entity = "Genus", 
                    sample_info = mgx_sample)  
mgx_sp<- get_mgx_Ids(sp_mgx_ab, entity_Node = species_Node, entity="Species", 
                     sample_info = mgx_sample)
mgx_st<- get_mgx_Ids(st_mgx_ab, entity_Node = strain_Node, entity="Strain", 
                     sample_info = mgx_sample)

# Reshape the previous created dataframes from wide to long foramt. 
reshp_mgx_k<- reshape_taxonDat(mgx_k, by_melt="tx_Id:ID")
reshp_mgx_p<- reshape_taxonDat(mgx_p, by_melt="tx_Id:ID")
reshp_mgx_c<- reshape_taxonDat(mgx_c, by_melt="tx_Id:ID")
reshp_mgx_o<- reshape_taxonDat(mgx_o, by_melt="tx_Id:ID")
reshp_mgx_f<- reshape_taxonDat(mgx_f, by_melt="tx_Id:ID")
reshp_mgx_g<- reshape_taxonDat(mgx_g, by_melt="tx_Id:ID")
reshp_mgx_sp<- reshape_taxonDat(mgx_sp, by_melt="tx_Id:ID")
reshp_mgx_st<- reshape_taxonDat(mgx_st, by_melt="tx_Id:ID")

# Once this complete, let's reformat the dataframe in Neo4j format. 
neo4j_K_mgx<- edge_Neo4j(reshp_mgx_k,
                          edge_label="HAS_MGX_K_ABUNDANCE")
neo4j_P_mgx<- edge_Neo4j(reshp_mgx_p,
                          edge_label="HAS_MGX_P_ABUNDANCE")
neo4j_C_mgx<- edge_Neo4j(reshp_mgx_c,
                          edge_label="HAS_MGX_C_ABUNDANCE")
neo4j_O_mgx<- edge_Neo4j(reshp_mgx_o,
                          edge_label="HAS_MGX_O_ABUNDANCE")
neo4j_F_mgx<- edge_Neo4j(reshp_mgx_f,
                          edge_label="HAS_MGX_F_ABUNDANCE")
neo4j_G_mgx<- edge_Neo4j(reshp_mgx_g,
                          edge_label="HAS_MGX_G_ABUNDANCE")
neo4j_Sp_mgx<- edge_Neo4j(reshp_mgx_sp,
                         edge_label="HAS_MGX_Sp_ABUNDANCE")
neo4j_St_mgx<- edge_Neo4j(reshp_mgx_st,
                         edge_label="HAS_MGX_St_ABUNDANCE")

################################################################################
# Main Program: Create a Sample-MVX-Taxon Edge                                 #
# Input: Kingdom, Phylum, Class, Order, Family, Genus, Species, Strain Nodes,  #
#        sample_info, and mvx_tax                                              #
# Return:                                                                      #
################################################################################
# Split 'taxonomy' column into number of taxons in taxonomy lineage.
split_taxon<- data.frame(str_split_fixed(rownames(mvx_tax), "\\|", n=7))
colnames(split_taxon)<- c("Kingdom","Phylum","Class","Order","Family","Genus",
                          "Species")

# Make it Null for rownames in mvx_ab dataframe. 
mvx_ab<- mvx_tax 
rownames(mvx_ab)<- NULL
# We want to aggregate the table based on the taxonomic rank. 
species_mvx<- aggregate_dat(mvx_ab, split_taxon, columnName = "Species")
genus_mvx<- aggregate_dat(mvx_ab, split_taxon, columnName = "Genus")
family_mvx<- aggregate_dat(mvx_ab, split_taxon, columnName = "Family")
order_mvx<- aggregate_dat(mvx_ab, split_taxon, columnName = "Order")
class_mvx<- aggregate_dat(mvx_ab, split_taxon, columnName = "Class")
phylum_mvx<- aggregate_dat(mvx_ab, split_taxon, columnName = "Phylum")
kingdom_mvx<- aggregate_dat(mvx_ab, split_taxon, columnName="Kingdom")

# Subset MVX samples in sample_info dataframe. 
mvx_sample<- sample_info[sample_info$data_type == "viromics",]
# Replace sample/taxon information with Ids. 
k_mvx_ids<- get_16s_Ids(kingdom_mvx, kingdom_Node, entity="Kingdom", 
                        sample_info = mvx_sample)
p_mvx_ids<- get_16s_Ids(phylum_mvx, phylum_Node, entity="Phylum", 
                        sample_info = mvx_sample)
c_mvx_ids<- get_16s_Ids(class_mvx, class_Node, entity="Class", 
                        sample_info = mvx_sample)
o_mvx_ids<- get_16s_Ids(order_mvx, order_Node, entity="Order", 
                        sample_info = mvx_sample)
f_mvx_ids<- get_16s_Ids(family_mvx, family_Node, entity="Family", 
                        sample_info = mvx_sample)
g_mvx_ids<- get_16s_Ids(genus_mvx, genus_Node, entity="Genus", 
                        sample_info = mvx_sample)
sp_mvx_ids<- get_16s_Ids(species_mvx, species_Node, entity="Species", 
                         sample_info = mvx_sample)

# Reshape the previous created dataframes from wide to long foramt. 
reshp_mvx_k<- reshape_taxonDat(k_mvx_ids, by_melt="tx_Id:ID")
reshp_mvx_p<- reshape_taxonDat(p_mvx_ids, by_melt="tx_Id:ID")
reshp_mvx_c<- reshape_taxonDat(c_mvx_ids, by_melt="tx_Id:ID")
reshp_mvx_o<- reshape_taxonDat(o_mvx_ids, by_melt="tx_Id:ID")
reshp_mvx_f<- reshape_taxonDat(f_mvx_ids, by_melt="tx_Id:ID")
reshp_mvx_g<- reshape_taxonDat(g_mvx_ids, by_melt="tx_Id:ID")
reshp_mvx_sp<- reshape_taxonDat(sp_mvx_ids, by_melt="tx_Id:ID")

# Once this complete, let's reformat the dataframe in Neo4j format. 
neo4j_K_mvx<- edge_Neo4j(reshp_mvx_k,
                         edge_label="HAS_MVX_K_ABUNDANCE")
neo4j_P_mvx<- edge_Neo4j(reshp_mvx_p,
                         edge_label="HAS_MVX_P_ABUNDANCE")
neo4j_C_mvx<- edge_Neo4j(reshp_mvx_c,
                         edge_label="HAS_MVX_C_ABUNDANCE")
neo4j_O_mvx<- edge_Neo4j(reshp_mvx_o,
                         edge_label="HAS_MVX_O_ABUNDANCE")
neo4j_F_mvx<- edge_Neo4j(reshp_mvx_f,
                         edge_label="HAS_MVX_F_ABUNDANCE")
neo4j_G_mvx<- edge_Neo4j(reshp_mvx_g,
                         edge_label="HAS_MVX_G_ABUNDANCE")
neo4j_Sp_mvx<- edge_Neo4j(reshp_mvx_sp,
                          edge_label="HAS_MVX_Sp_ABUNDANCE")

# Write the Taxon Nodes Output file to a .csv file 
fileN1= paste("Kingdom", "node", outName, "csv", sep=".")
write.csv(kingdom_Node, fileN1, row.names=F, quote=T, sep="\t")
fileN2= paste("Phylum", "node", outName, "csv", sep=".")
write.csv(phylum_Node, fileN2, row.names=F, quote=T, sep="\t")
fileN3= paste("Class", "node", outName, "csv", sep=".")
write.csv(class_Node, fileN3, row.names=F, quote=T, sep="\t")
fileN4= paste("Order", "node", outName, "csv", sep=".")
write.csv(order_Node, fileN4, row.names=F, quote=T, sep="\t")
fileN5= paste("Family", "node", outName, "csv", sep=".")
write.csv(family_Node, fileN5, row.names=F, quote=T, sep="\t")
fileN6= paste("Genus", "node", outName, "csv", sep=".")
write.csv(genus_Node, fileN6, row.names=F, quote=T, sep="\t")
fileN7= paste("Species", "node", outName, "csv", sep=".")
write.csv(species_Node, fileN7, row.names=F, quote=T, sep="\t")
fileN8= paste("Strain", "node", outName, "csv", sep=".")
write.csv(strain_Node, fileN8, row.names=F, quote=T, sep="\t")

# Write the sample-16sTaxon file to a .csv file 
fileN9= paste("sample_16sBiopsy_K","edge", outName, "csv", sep=".")
write.csv(neo4j_K_edge, fileN9, row.names=F, quote=T, sep=",")
fileN10= paste("sample_16sBiopsy_P","edge", outName, "csv", sep=".")
write.csv(neo4j_P_edge, fileN10, row.names=F, quote=T, sep=",")
fileN11= paste("sample_16sBiopsy_C","edge", outName, "csv", sep=".")
write.csv(neo4j_C_edge, fileN11, row.names=F, quote=T, sep=",")
fileN12= paste("sample_16sBiopsy_O","edge", outName, "csv", sep=".")
write.csv(neo4j_O_edge, fileN12, row.names=F, quote=T, sep=",")
fileN13= paste("sample_16sBiopsy_F","edge", outName, "csv", sep=".")
write.csv(neo4j_F_edge, fileN13, row.names=F, quote=T, sep=",")
fileN14= paste("sample_16sBiopsy_G","edge", outName, "csv", sep=".")
write.csv(neo4j_G_edge, fileN14, row.names=F, quote=T, sep=",")

# Write the sample-MGX file to a .csv file
fileN15= paste("sample_MGX_K","edge", outName, "csv", sep=".")
write.csv(neo4j_K_mgx, fileN15, row.names=F, quote=T, sep=",")
fileN16= paste("sample_MGX_P","edge", outName, "csv", sep=".")
write.csv(neo4j_P_mgx, fileN16, row.names=F, quote=T, sep=",")
fileN17= paste("sample_MGX_C","edge", outName, "csv", sep=".")
write.csv(neo4j_C_mgx, fileN17, row.names=F, quote=T, sep=",")
fileN18= paste("sample_MGX_O","edge", outName, "csv", sep=".")
write.csv(neo4j_O_mgx, fileN18, row.names=F, quote=T, sep=",")
fileN19= paste("sample_MGX_F","edge", outName, "csv", sep=".")
write.csv(neo4j_F_mgx, fileN19, row.names=F, quote=T, sep=",")
fileN20= paste("sample_MGX_G","edge", outName, "csv", sep=".")
write.csv(neo4j_G_mgx, fileN20, row.names=F, quote=T, sep=",")
fileN21= paste("sample_MGX_Sp","edge", outName, "csv", sep=".")
write.csv(neo4j_Sp_mgx, fileN21, row.names=F, quote=T, sep=",")
fileN22= paste("sample_MGX_St","edge", outName, "csv", sep=".")
write.csv(neo4j_St_mgx, fileN22, row.names=F, quote=T, sep=",")

# Write the sample-MVX file to a .csv file
fileN23= paste("sample_MVX_K","edge", outName, "csv", sep=".")
write.csv(neo4j_K_mvx, fileN23, row.names=F, quote=T, sep=",")
fileN24= paste("sample_MVX_P","edge", outName, "csv", sep=".")
write.csv(neo4j_P_mvx, fileN24, row.names=F, quote=T, sep=",")
fileN25= paste("sample_MVX_C","edge", outName, "csv", sep=".")
write.csv(neo4j_C_mvx, fileN25, row.names=F, quote=T, sep=",")
fileN26= paste("sample_MVX_O","edge", outName, "csv", sep=".")
write.csv(neo4j_O_mvx, fileN26, row.names=F, quote=T, sep=",")
fileN27= paste("sample_MVX_F","edge", outName, "csv", sep=".")
write.csv(neo4j_F_mvx, fileN27, row.names=F, quote=T, sep=",")
fileN28= paste("sample_MVX_G","edge", outName, "csv", sep=".")
write.csv(neo4j_G_mvx, fileN28, row.names=F, quote=T, sep=",")
fileN29= paste("sample_MVX_Sp","edge", outName, "csv", sep=".")
write.csv(neo4j_Sp_mvx, fileN29, row.names=F, quote=T, sep=",")