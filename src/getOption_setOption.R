#-------------------------------------------------------------#
# Goal: Getoption seption                                     #
# - Make the globalenvironment 'stringsAsFactors' == FALSE    #
#                                                             #
# Name: Suyeon Kim                                            #
# Creation Date: Oct-08-2021                                  #
# Last updated: Oct-08-2021                                   #
#-------------------------------------------------------------#

## If Global environment for "stringsAsFactors: is equal to FALSE, it is good. 
## Otherwise, "stringsAsFactors" is equal to TRUE, change the environment to FALSE
if(getOption("stringsAsFactors")) options(stringsAsFactors=F)

getOption("stringsAsFactors")

