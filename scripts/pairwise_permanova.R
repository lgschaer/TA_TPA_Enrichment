#PAIRWISE PERMANOVA

library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

# Adonis test, all samples
head(sdata3)

adonisData <- sequence_table %>%
  rownames_to_column(var = "SampleID") %>%                                          #change rownames to a column so there is a common variable to join by
  left_join(sdata3, by = "SampleID") %>%                                            #join sample data to the sequence table
  select(-c("SampleID", "Environment", "Replicate", "Transfer", "Sample_Type")) %>% #remove all metadata columns except the one to be used to compare
  select(Substrate, everything())                                                   #rearrange so substrate column is first
dim(adonisData)

#similarity euclidean from vegdist and holm correction
pairwise.adonis(x=adonisData[,2:12320],factors=adonisData$Substrate,sim.function='vegdist',
                sim.method='bray',p.adjust.m='holm')
