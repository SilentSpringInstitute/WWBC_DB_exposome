library(tidyverse)

####### POSITIVE MODE #######

#import 2021 frag data

node2021_pos <- read.csv("~/Box/SHE lab/WWBC/Data analysis/Non-targeted analysis/2021Fragmentation matching/node_attributes_table_pos1.tsv", header = T, sep = "\t") #Postive Mode


node2021_pos$cluster.index<-as.factor(node2021_pos$cluster.index)

##Unlist chemicals ID
node2021_pos$ConsensusID<-as.character(node2021_pos$ConsensusID)
node2021_pos$ConsensusScore<-as.character(node2021_pos$ConsensusScore)
node2021_pos$ConsensusSMILES<-as.character(node2021_pos$ConsensusSMILES)

#this code is to transform the shape of the dataset to have one line per consensus ID  
node2021_pos1<- node2021_pos %>%   
  transform(ConsensusID=strsplit(ConsensusID,","),
            ConsensusScore=strsplit(ConsensusScore,","),
            ConsensusSMILES=strsplit(ConsensusSMILES,",")) %>%
  unnest(c(ConsensusID, ConsensusScore, ConsensusSMILES))

#in addition this line of code keeps only those with a consensus ID. it removes fragmentation peaks that do not match. 

node2021_pos1<-node2021_pos1[,c(1,2,4,5,17:19)]  

node2021_pos1 <- node2021_pos1 %>%
  mutate(parent.mass= round(parent.mass, 1))

#add year to the colname
colnames_temp <- colnames(node2021_pos1)
cor_names_paste <- paste0(colnames_temp, ".21")
node2021_posmerge <- node2021_pos1
colnames(node2021_posmerge) <- cor_names_paste

#import 2020 fragmentation data
node2020_pos <- read.csv("~/Box/SHE lab/WWBC/Data analysis/Non-targeted analysis/DoNotUse_2020Fragmentation/node_attributes_table_positive.tsv", header = T, sep = "\t") #Postive Mode


node2020_pos$cluster.index<-as.factor(node2020_pos$cluster.index)

##Unlist chemicals ID
node2020_pos$ConsensusID<-as.character(node2020_pos$ConsensusID)
node2020_pos$ConsensusScore<-as.character(node2020_pos$ConsensusScore)
node2020_pos$ConsensusSMILES<-as.character(node2020_pos$ConsensusSMILES)

#this code is to transform the shape of the dataset to have one line per consensus ID  
node2020_pos1<- node2020_pos %>%   
  transform(ConsensusID=strsplit(ConsensusID,","),
            ConsensusScore=strsplit(ConsensusScore,","),
            ConsensusSMILES=strsplit(ConsensusSMILES,",")) %>%
  unnest(c(ConsensusID, ConsensusScore, ConsensusSMILES))

#in addition this line of code keeps only those with a consensus ID. it removes fragmentation peaks that do not match. 

node2020_pos1<-node2020_pos1[,c(1,2,4,5,17:19)]  

node2020_pos1 <- node2020_pos1 %>%
  mutate(parent.mass= round(parent.mass, 1))

colnames_temp <- colnames(node2020_pos1)
cor_names_paste <- paste0(colnames_temp, ".20")
node2020_posmerge <- node2020_pos1
colnames(node2020_posmerge) <- cor_names_paste


overlap_pos <- merge(node2021_posmerge, node2020_posmerge, by.x   = "parent.mass.21", by.y= "parent.mass.20")
 write.csv(overlap_pos, "Box/SHE lab/WWBC/Data analysis/Non-targeted analysis/2021Fragmentation matching/QA_QC_fragmentation_runs/fragmentationdata_merged_2020_2021.csv")

overlap_pos %>%
   group_by(parent.mass.21)%>%
   summarise(count=n_distinct(parent.mass.21))

remove(list=ls())
####### NEGATIVE MODE #######

#import 2021 frag data

node2021_neg <- read.csv("~/Box/SHE lab/WWBC/Data analysis/Non-targeted analysis/2021Fragmentation matching/node_attributes_table_neg.tsv", header = T, sep = "\t") #negtive Mode


node2021_neg$cluster.index<-as.factor(node2021_neg$cluster.index)

##Unlist chemicals ID
node2021_neg$ConsensusID<-as.character(node2021_neg$ConsensusID)
node2021_neg$ConsensusScore<-as.character(node2021_neg$ConsensusScore)
node2021_neg$ConsensusSMILES<-as.character(node2021_neg$ConsensusSMILES)

#this code is to transform the shape of the dataset to have one line per consensus ID  
node2021_neg1<- node2021_neg %>%   
  transform(ConsensusID=strsplit(ConsensusID,","),
            ConsensusScore=strsplit(ConsensusScore,","),
            ConsensusSMILES=strsplit(ConsensusSMILES,",")) %>%
  unnest(c(ConsensusID, ConsensusScore, ConsensusSMILES))

#in addition this line of code keeps only those with a consensus ID. it removes fragmentation peaks that do not match. 

node2021_neg1<-node2021_neg1[,c(1,2,4,5,17:19)]  

node2021_neg1 <- node2021_neg1 %>%
  mutate(parent.mass= abs(round(parent.mass, 2)))

#add year to the colname
colnames_temp <- colnames(node2021_neg1)
cor_names_paste <- paste0(colnames_temp, ".21")
node2021_negmerge <- node2021_neg1
colnames(node2021_negmerge) <- cor_names_paste

#import 2020 fragmentation data
node2020_neg <- read.csv("~/Box/SHE lab/WWBC/Data analysis/Non-targeted analysis/DoNotUse_2020Fragmentation/node_attributes_table_negative.tsv", header = T, sep = "\t") #negtive Mode


node2020_neg$cluster.index<-as.factor(node2020_neg$cluster.index)

##Unlist chemicals ID
node2020_neg$ConsensusID<-as.character(node2020_neg$ConsensusID)
node2020_neg$ConsensusScore<-as.character(node2020_neg$ConsensusScore)
node2020_neg$ConsensusSMILES<-as.character(node2020_neg$ConsensusSMILES)

#this code is to transform the shape of the dataset to have one line per consensus ID  
node2020_neg1<- node2020_neg %>%   
  transform(ConsensusID=strsplit(ConsensusID,","),
            ConsensusScore=strsplit(ConsensusScore,","),
            ConsensusSMILES=strsplit(ConsensusSMILES,",")) %>%
  unnest(c(ConsensusID, ConsensusScore, ConsensusSMILES))

#in addition this line of code keeps only those with a consensus ID. it removes fragmentation peaks that do not match. 

node2020_neg1<-node2020_neg1[,c(1,2,4,5,17:19)]  

node2020_neg1 <- node2020_neg1 %>%
  mutate(parent.mass= round(parent.mass, 1))

colnames_temp <- colnames(node2020_neg1)
cor_names_paste <- paste0(colnames_temp, ".20")
node2020_negmerge <- node2020_neg1
colnames(node2020_negmerge) <- cor_names_paste


overlap_neg <- merge(node2021_negmerge, node2020_negmerge, by.x   = "parent.mass.21", by.y= "parent.mass.20")
write.csv(overlap_neg, "Box/SHE lab/WWBC/Data analysis/Non-targeted analysis/2021Fragmentation matching/QA_QC_fragmentation_runs/Negative_fragmentationdata_merged_2020_2021.csv")

overlap_neg %>%
  group_by(parent.mass.21)%>%
  summarise(count=n_distinct(parent.mass.21))  
