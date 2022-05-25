# Multi-Omics Graph Database Model  

Multi-Omics integrated graph database. Our graph database provides a useful analytical platform to help a better understanding of integrated data. The current version of model is specifically used the multi-omcis datasets from 'The Inflammatory Bowel Disease Multi'omics Database (IBDMDB)' (Lloyd-Price et al. 2019). Neo4j is a graph database management system that uses a graph database model. Our multi-omics graph database is highly queryalbe and provides several case studies that utilize simple as well as complex queries for studying IBD. Each step for populating our graph DB will be deatiled below. 

Contact the author at suyeonkim [at] unomaha.edu. This version has been tested for OSX. 
## Steps for populating the multi-omics graph database

### Pre-requisites
Please make sure you have Conda installed. 
neo4r is an R package that can be run as an R function. 

#### Installation N4o4j and running a local database 
To populate our graph database model, it needs to able to access a Neo4j database. Here, the instracutions below are to run the Neo4j via a server on macOS compouters. 

From command line 
  1. Open up your terminal and create a 'demo' folder.
  ```
  mkdir demo
  ```
  2. Download the macOS version of Neo4j listed under community server 
  ```
  wget https://neo4j.com/artifact.php?name=neo4j-community-4.3.13-unix.tar.gz
  tar -xf neo4j-community-4.3.13-unix.tar.gz
  ```
  3. 

#### Download IBDMDB Data 
1. Go visit the [The Inflammatory Bowel Disease Multi'omics Database](https://ibdmdb.org/) website. 
2. On the `Download Data` tab, we use batch downloading using FTP transfer tool.  
3. Copy `ftp.broadinstitute.org` address and paste it in new tab, with adding 'ftp` instead 'http'.
4. Once open the folder, use given `username/password` to access the datasets. 


Download the macOS version of Neo4j listed under Community server
Unzip the Neo4j download in a folder of your choice
In command line, navigate to the Neo4j folder where you unzipped the download
Run ./bin/neo4j console to start the server
Run CTRL+C to stop the server


##### From command line 







### Running Neo4j on a server 

### Running Neo4j in R 




Steps for populating the database and running the server
### raw files

