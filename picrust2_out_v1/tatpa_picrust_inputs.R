### Making inputs for picrust

#load libraries
library(phyloseq)
library(tidyverse)
#install.packages("phylotools", dependencies = TRUE)
library(phylotools)
#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install("biomformat")
library(biomformat)


#load sequence table
sequence_table <- readRDS("/Users/lauraschaerer/Desktop/TATPAEnrichment/seqtab.rds")
sequence_table[1:5,1:5]
dim(sequence_table)

pre_fna_summary <- sequence_table %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "seq.text") %>%
  mutate(seq.name = paste("sp", 1:12003, sep="")) %>%
  select(c("seq.name", "seq.text"))
head(pre_fna_summary)
dim(pre_fna_summary)


phylotools::dat2fasta(pre_fna_summary, outfile = "/Users/lauraschaerer/Desktop/TATPAEnrichment/arsenic.fna")

## Making the biom table

ps <- readRDS("/Users/lauraschaerer/Desktop/TATPAEnrichment/ps_nomito_nochloro_noinocula_T4.rds")
ps
#View(sample_data(ps))

# making otu.biom
count_table <- (otu_table(ps))
count_table[1:5,1:5]

#transpose
t_count_table <- t(count_table)
t_count_table[,5:1][1:5,]
class(t_count_table)

#write table
write.table(t_count_table, "/Users/lauraschaerer/Desktop/TATPAEnrichment/t_count_table.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
write_biom(t_count_table, "/Users/lauraschaerer/Desktop/TATPAEnrichment/table.biom")

##next step: use putty to convert to .biom table

#Processing the picrust output

#load packages
library(tidyverse)
library(csv)
library(phyloseq)
library(DESeq2) #differential abundance
library("KEGGREST") #getting KO annotations



#load phyloseq object
ps2 <- read_rds("/Users/lauraschaerer/Desktop/TATPAEnrichment/ps_nomito_nochloro_noinocula_T4.rds") %>%
  subset_samples(Substrate != "Control")
ps2
#View(otu_table(ps)[1:5,1:5])
#ps2 <- subset_samples(ps, Timepoint == "T7")
unique(sample_data(ps2)$Substrate)


#load sample data and format for phyloseq
sdata <- as_tibble(sample_data(ps2)) %>%
  mutate(SampleID2 = SampleID) %>%
  column_to_rownames(var = "SampleID2") #%>%
  #filter(!is.na(Fern.type))
head(sdata)
#View(sdata)
#unique(sdata$Fern.type)
dim(sdata)


#load picrust THIS WILL BE THE SEQUENCE TABLE
pred_mg <- read_tsv("/Users/lauraschaerer/Desktop/TATPAEnrichment/picrust2_out_11142024/KO_metagenome_out/pred_metagenome_unstrat.tsv") 

pmg2 <- pred_mg %>%
  column_to_rownames(var = "function") %>%
  t()
pmg2[1:10,1:10]
rownames(pmg2)


# getting a list of gene names to match KO numbers THIS WILL BE THE TAXA TABLE
KO <- pred_mg %>%
  column_to_rownames(var = "function") %>%
  rownames_to_column(var = "KO") %>%
  select(KO)
head(KO)
dim(KO)



#make phyloseq object
samdata = sample_data(sdata)                                      #define sample data
colnames(pmg2) <- NULL                                          #remove column names from "nonzero"
seqtab = otu_table(pmg2, taxa_are_rows = FALSE)                 #define sequence table
taxtab = tax_table(KO)                                     #define taxa table
rownames(taxtab) <- NULL                                           #remove rownames from taxa table


sample_names(samdata)
sample_names(seqtab)

psGene = phyloseq(otu_table(seqtab), tax_table(taxtab), sample_data(samdata)) 

#View(sample_data(phyloseq_object_all))
unique(sample_data(psGene)$Substrate)

psGene <- subset_samples(psGene, sample_sums(psGene) > 0)
psGene <- subset_samples(psGene, taxa_sums(psGene) > 0)
psGene

####### Looking at interesting enzymes
library(KEGGREST)

#interesting enzymes list
test <- keggGet(c("ec00624",   #Polycyclic aromatic hydrocarbon degradation
                  "ec00362",   #Benzoate degradation
                  "ec00627"    #Aminobenzoate degradation
)) 
test[[1]]$ENZYME

e1 <- as.data.frame(test[[1]]$ENZYME) %>% 
  mutate(EC_number = test[[1]]$ENZYME, 
         General_Pathway = test[[1]]$NAME)

e2 <- as.data.frame(test[[2]]$ENZYME) %>% 
  mutate(EC_number = test[[2]]$ENZYME, 
         General_Pathway = test[[2]]$NAME)

e3 <- as.data.frame(test[[3]]$ENZYME) %>% 
  mutate(EC_number = test[[3]]$ENZYME, 
         General_Pathway = test[[3]]$NAME)
e4 <- data.frame(EC_number = c("3.5.1.4"), General_Pathway = c("Amidase"))
head(e1)


interesting_enzymes <- e1 %>%
  full_join(e2, by = c("EC_number", "General_Pathway")) %>%
  full_join(e3, by = c("EC_number", "General_Pathway")) %>%
  full_join(e4, by = c("EC_number", "General_Pathway")) %>%
  dplyr::select(c(EC_number, General_Pathway))
head(interesting_enzymes)
dim(interesting_enzymes)
interesting_EC <- interesting_enzymes$EC_number

ec1 <- interesting_EC[1:100] %>% keggList() %>% as.data.frame()
View(ec1)
KO1 <- KO[1:100] %>% keggList() %>% as.data.frame()

KO_details2 <- KO %>%
  rownames_to_column(var = "KO") %>%
  separate(details, into = c("Gene_Abbreviation", "Gene_Name"), sep = "; ") %>%
  separate(Gene_Name, into = c("Gene_Name", "EC_Number"), sep = " \\[EC:") %>%
  separate(EC_Number, into = c("EC_Number", NA), sep = "\\]") %>%
  mutate(Gene_Prefix = str_extract(Gene_Abbreviation, "^..."),
         Variable = "x")
#separate(Gene_Abbreviation, into = c("Gene_Prefix", "NA"), sep = "\\]")
head(KO_details2)

genes <- psGene %>%
  tax_glom(taxrank = "ta1") %>%                          # Set to smallest taxonomic level you are interested in
  transform_sample_counts(function(x) {x/sum(x)} ) %>%   # Transform to rel. abundance
  psmelt()                                               # Melt to long format
head(genes)

KO_of_interest <- c("K01426", "K18074", "K18075", "K00448", "K00449", "K18076")

KO1 <- KO_of_interest %>% keggList() %>% as.data.frame()
colnames(KO1) <- c("details")
KO_details <- KO1 %>%
  rownames_to_column(var = "KO") %>%
  separate(details, into = c("Gene_Abbreviation", "Gene_Name"), sep = "; ") %>%
  separate(Gene_Name, into = c("Gene_Name", "EC_Number"), sep = " \\[EC:") %>%
  separate(EC_Number, into = c("EC_Number", NA), sep = "\\]") %>%
  mutate(Gene_Prefix = str_extract(Gene_Abbreviation, "^...")) %>%
  select(KO, Gene_Name, EC_Number)
head(KO_details)

genes_filt <- genes %>%
  filter(Sample_Type == "Enrichment" & Transfer == "T4") %>%
  mutate(KO = ta1,
         #KO = ifelse(Abundance < 0.001, NA, KO)
         ) %>%
  filter(KO %in% KO_of_interest) %>%
  group_by(Environment, Substrate, KO) %>%
  summarise(
    Abundance = sum(Abundance)
  ) %>%
  left_join(KO_details)
head(genes_filt)
#length(unique(genes_filt$KO))

ggplot(genes_filt, aes(x = Substrate, y = Abundance, fill = Gene_Name))+
  facet_grid(cols = vars(Environment))+
  geom_col(color = "black")+
  xlab("Substrate")+
  ylab("Relative Abundance")+
  theme_linedraw()+ 
  theme(axis.text.x = element_text(size = 20, angle = 0, vjust = 1, hjust = 0.5, color = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14, angle = 0, hjust = 1, vjust = 0.5, color = "black"),
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        legend.title = element_blank(),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text.x = element_text(size = 18, face = "bold", angle = 0),
        strip.text.y = element_text(size = 18, face = "bold", angle = 0))+
  guides(fill=guide_legend(nrow=5,byrow=TRUE))


