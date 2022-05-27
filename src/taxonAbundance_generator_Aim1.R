#-----------------------------------------------------------------------------------------#
# Goal: This script is for generating abundance matrix according to different             #
#       level of taxonomy                                                                 #
#                                                                                         #
# Date: Nov-01-2021                                                                       #
# Last updated: Apr-08-2022                                                               #
# Name: Suyeon Kim                                                                        #
#-----------------------------------------------------------------------------------------#
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


## Test Aggregate function. 
#aggregate_ab<- function(mvx_ab, split_taxon, columnName) {

## Create a lineage 
## ASK Ishwor [Character vs Numeric] 
#  columnIndex <- which(colnames(split_taxon)== columnName)

#  columns_to_select <- colnames(split_taxon)[1:(columnIndex-1)]

#  lineage<- subset(split_taxon, select=columns_to_select)

#  alineage <- apply(lineage, 1, paste, collapse='|')

## Combine lineage with original abundance data 
#  new_table<- cbind("ancestor"=alineage, "leaf"=split_taxon[,columnIndex],mvx_ab)
#  colnames(new_table)    

## Generate abundance matrix 
#  aggre_dat<- aggregate( new_table[,!(names(new_table) %in% colnames(lineage))], by=list(new_table$Kingdom, new_table$Phylum, new_table$Class, new_table$Order, new_table$Family, new_table$Genera), FUN=sum, na.rm=T)
#filt_table<- aggregate(.~Kingdom+Phylum+Class+Order+Family+Genera, data= new_table, sum, na.rm=T)  
#  agg_table <- aggregate(new_table[,c(-1,-2)], by=list(new_table$ancestor, new_table$leaf),FUN=sum, na.rm=T)

#  return(filt_table)
#}