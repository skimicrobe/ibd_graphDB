################################################################################
# Goal: This script is for generating abundance matrix according to different  #
#       level of taxonomy.                                                     #
# Author: Suyeon Kim                                                           #
# Date: Nov-01-2021                                                            #
# Last updated: Apr-08-2022                                                    #
################################################################################
aggregate_dat<- function(mvx_ab, split_taxon, columnName) {
  # Get the index of taxon of interest. 
  columnIndex<- which(colnames(split_taxon) == columnName)
  # Grep a lineage info (from 'Kingdom [index 1]' upto before taxon of interest. 
  columns_to_select <- colnames(split_taxon)[1:(columnIndex-1)]
  # Create a lineage info table. 
  lineage<- subset(split_taxon, select=columns_to_select)
  # Format the lineage info. 
  alineage<- apply(lineage, 1, paste, collapse='|')
  # Combine lineage info table with original abundance data 
  new_table<- cbind("ancestor" = alineage, "leaf" = split_taxon[,columnIndex], 
                    mvx_ab)
  colnames(new_table)
  # Aggregate the abundance matrix. 
  agg_table<- aggregate(new_table[,c(-1, -2)], by=list(new_table$ancestor, 
                       new_table$leaf), FUN=sum, na.rm=T)
  print(colnames(agg_table)[1:2])
  colnames(agg_table)[1:2]<- c("ancestral_lineage", columnName)
  print(agg_table[,1:5])
  agg_table[,1]<- gsub("\\| ", "|", agg_table[,1])
  agg_table[,2]<- gsub(" g__*","g__", agg_table[,2])
  agg_table[,2]<- gsub(" f__*","f__", agg_table[,2])
  agg_table[,2]<- gsub(" o__*","o__", agg_table[,2])
  agg_table[,2]<- gsub(" c__*","c__", agg_table[,2])
  agg_table[,2]<- gsub(" p__*","p__", agg_table[,2])
  print(colnames(agg_table)[1:2])
  return(agg_table)
} 

save_output<- function(abundance_dat, taxonName, data_type){
  
  ## Drop the lineage column 
  #final_ab<- abundance_dat[,-1]
  ## Save the file 
  fileN= paste(taxonName, data_type, "abundance", "csv", sep=".")
  write.csv(abundance_dat, fileN, row.names=F, quote=F)
}
