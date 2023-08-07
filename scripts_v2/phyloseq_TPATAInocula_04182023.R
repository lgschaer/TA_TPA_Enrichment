
### Create a phyloseq object with the new inocula sequencing data from Zymo to look at the decrease in diversity during the enrichment process.

library(phyloseq)
library(tidyverse)
library(csv)

###load the original phyloseq object and subset it to only include inocula samples
ps_full <- readRDS("/home/lgschaer/old/Plastic_Deg/TA_TPA_Enrichment_09062021/ps_unrarefied_09072021.rds")
head(sample_data(ps_full))

ps_full_inocula <- subset_samples(ps_full, XXXX == "Inocula")



###create a phyloseq object with the new sequencing data from Zymo
#load metadata

sdata <- as.csv("/home/lgschaer/old/Plastic_Deg/TA_TPA_Enrichment_09062021/Inocula_Seq_Zymo/Inocula_Metadata_04132023.csv", row.names = 1, header = TRUE, sep = ",", check.names = TRUE, stringsAsFactors = TRUE)
head(sdata)

#load sequence table
sequence_table <- readRDS("/home/lgschaer/old/Plastic_Deg/TA_TPA_Enrichment_09062021/Inocula_Seq_Zymo/dada2_out/seqtab.rds")    
sequence_table[1:3,1:3]

#1 remove column names
colnames(sequence_table) <- NULL
sequence_table[1:5,1:5]

#2 convert to data frame
sequence_table <- as.data.frame(sequence_table)
class(sequence_table)

#dim(seq_table)
dim(sdata)
dim(sequence_table)

#7 remove metadata columns
subset_all <- sequence_table

#8 remove column names from the sequence table
colnames(subset_all) <- NULL
subset_all[1:5,1:5]

#load taxa table
taxa <- readRDS("/home/lgschaer/old/Plastic_Deg/TA_TPA_Enrichment_09062021/Inocula_Seq_Zymo/dada2_out/taxa.rds")
taxa[1:5,1:5]

#make phyloseq object
samdata = sample_data(sdata3)
seqtab = otu_table(subset_all, taxa_are_rows = FALSE)
taxtab = tax_table(taxa)

sample_names(samdata)
sample_names(seqtab)
sample_names(taxtab)

taxa_names(seqtab)
sample_names(taxtab)

#combine all components into a phyloseq object
ps = phyloseq(otu_table(seqtab), tax_table(taxtab), sample_data(samdata))
ps

head(sample_data(ps))

write_rds(ps, "/home/lgschaer/old/Plastic_Deg/TA_TPA_Enrichment_09062021/Inocula_Seq_Zymo/phyloseq_out/ps_unrarefied_04182023.rds")




#now merge the two phyloseq objects together

ps_new <- merge_phyloseq(ps, ps_full_inocula)

#Filter out eukaryotes and mitochondria
ps2 <- ps_new %>%
  subset_taxa(
    Kingdom == "Bacteria" &                   #only bacteria
      Family  != "Mitochondria" &             #filter out mitochondria
      Class   != "Chloroplast"                #filter out chloroplasts
  )
ps2

### normalize

samplesover1000_all <- subset_samples(ps2, sample_sums(ps2) > 1000)

any(taxa_sums(samplesover1000_all) == 0)

sum(taxa_sums(samplesover1000_all) == 0)

prune_samplesover1000_all <- prune_taxa(taxa_sums(samplesover1000_all) > 0, samplesover1000_all)

any(taxa_sums(prune_samplesover1000_all) == 0)

#for reproducible data
set.seed(81)

rarefy_samplesover1000_all <- rarefy_even_depth(prune_samplesover1000_all, rngseed= 81, sample.size = min(sample_sums(prune_samplesover1000_all)))

pps<- rarefy_samplesover1000_all
pps

### making figures

#pps_inocula <- subset_samples(pps, Substrate == "Inocula")

#pps_inocula %>%                                                              #phyloseq object
 # plot_richness(
  #  x = "Environment",                                                 #compare diversity of datatype
   # measures = c("Observed", "Shannon")) +                           #choose diversity measures
#  geom_violin(aes(fill = Environment), show.legend = FALSE)+           #make violin plot, set fill aes to sampletype
#  geom_boxplot(width=0.1)# +                                          #add boxplot, set width
  #theme_classic()+                                                   #change theme to classic
  #xlab(NULL)+                                                        #no label on x-axis
  #theme(axis.text.y.left = element_text(size = 20),                  #adjust y-axis text
   #     axis.text.x = element_text(size = 20, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
    #    axis.title.y = element_text(size = 20))+                     #adjust y-axis title
  #theme(strip.text = element_text(face = "bold", size = 25))+        #adjust headings
  #scale_fill_manual(values = sample_colors)+                         #set fill colors
  #scale_x_discrete(                                                  #change x-axis labels
    #limits = sample_types)+                   
  #theme(plot.title=element_text(size = 25, face = "bold", hjust = 0.5)) #change title size, face and position


#sample colors
sample_colors <- c("Inocula" = "blue", "Control" = "lightblue", "TPA" = "maroon", "TA" = "olivedrab")
transfer_colors <- c("I" = "blue4", "T1" = "blue3", "T2" = "blue2", "T3" = "blue1", "T4" = "lightblue")

sample_types <- c("TPA", "TA", "Control")

#violin plot
pps %>%                                                              #phyloseq object
  plot_richness(
    x = "Transfer",                                                 #compare diversity of datatype
    measures = c("Observed", "Shannon")) +                           #choose diversity measures
  #facet_grid(cols = vars(Substrate_Label))+
  geom_boxplot(aes(fill = Substrate), show.legend = FALSE)+           #make violin plot, set fill aes to sampletype
  #geom_violin(width=0.1) +                                          #add boxplot, set width
  theme_linedraw()+                                                   #change theme to classic
  xlab(NULL)+                                                        #no label on x-axis
  theme(axis.text.y.left = element_text(size = 20),                  #adjust y-axis text
        axis.text.x = element_text(size = 20, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 20))+                     #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 25))+        #adjust headings
  scale_fill_manual(values = sample_colors)+                         #set fill colors
 # scale_x_discrete(                                                  #change x-axis labels
  #  limits = sample_types)+                   
  theme(plot.title=element_text(size = 25, face = "bold", hjust = 0.5)) #change title size, face and position

#Diversity bar plot with error bars
pps_filtered <- subset_samples(pps, Substrate != "Inocula" & Environment != "Inocula")

pps_diversity <- estimate_richness(pps, measures = c("Observed", "Shannon")) %>%
  rownames_to_column(var = "SampleID") %>%
  left_join(sdata3, by = "SampleID") %>%
  group_by(Substrate) %>%
  summarise(
    avgObs = mean(Observed),
    sdObs = sd(Observed),
    avgShan = mean(Shannon),
    sdShan = sd(Shannon))
head(pps_diversity)
#dim(pps_diversity)

pps_diversity <- estimate_richness(pps_filtered, measures = c("Observed", "Shannon")) %>%
  rownames_to_column(var = "SampleID") %>%
  left_join(sdata3, by = "SampleID") %>%
  group_by(Substrate, Substrate_Label, Environment, Transfer) %>%
  summarise(
    avgObs = mean(Observed),
    sdObs = sd(Observed),
    avgShan = mean(Shannon),
    sdShan = sd(Shannon))
head(pps_diversity)
#dim(pps_diversity)

colors <- c("lightblue", "maroon", "olivedrab")

ggplot(pps_diversity, aes(x = Substrate, y = avgObs))+
  #geom_errorbar(aes(ymin = minObs, ymax = maxObs), color = "black", position = "dodge") +
  geom_col(aes(fill = Substrate), show.legend = FALSE, position = "dodge", color = "black")+          #make violin plot, set fill aes to sampletype                                          #add boxplot, set width
  facet_grid(rows = vars(Environment), cols = vars(Transfer))+
  theme_linedraw()+                                                   #change theme to classic
  ylab("Observed ASVs")+ 
  xlab(NULL)+
  scale_fill_manual(values = colors)+  
  theme(axis.text.y.left = element_text(size = 12),                                        #adjust y-axis text
        axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 1, angle = 0),          #adjust x-axis label position
        axis.title.y = element_text(size = 20),
        strip.text = element_text(face = "bold", size = 20))                               #change title size, face and position

pps <- subset_samples(pps, Substrate != "Blank")
pps_filt <- subset_samples(pps, Substrate != "Inocula")

pps_diversity_filt <- estimate_richness(pps, measures = c("Observed", "Shannon")) %>%
  rownames_to_column(var = "SampleID") %>%
  left_join(sdata3, by = "SampleID") #%>%
  #pivot_longer(cols = c("Observed", "Shannon"), values_to = "Diversity", names_to = "Metric")
head(pps_diversity_filt)

ggplot(pps_diversity_filt, aes(x = Transfer, y = Shannon))+
  geom_boxplot(aes(fill = Substrate), show.legend = FALSE, position = "dodge", color = "black")+          #make violin plot, set fill aes to sampletype                                          #add boxplot, set width
  facet_nested(cols = vars(Substrate_Label), scales = "free_y", space = "free_y")+
  theme_linedraw()+                                                   #change theme to classic
  ylab("Shannon Diversity")+ 
  xlab(NULL)+
  scale_fill_manual(values = colors)+  
  theme(axis.text.y.left = element_text(size = 12),                                        #adjust y-axis text
        axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 1, angle = 0),          #adjust x-axis label position
        axis.title.y = element_text(size = 20),
        strip.text = element_text(face = "bold", size = 20))  

ggplot(pps_diversity_filt, aes(x = Transfer, y = Observed))+
  geom_boxplot(aes(fill = Substrate), show.legend = FALSE, position = "dodge", color = "black")+          #make violin plot, set fill aes to sampletype                                          #add boxplot, set width
  facet_nested(cols = vars(Substrate_Label), scales = "free_y", space = "free_y")+
  theme_linedraw()+                                                   #change theme to classic
  ylab("Number of ASVs")+ 
  xlab(NULL)+
  scale_fill_manual(values = colors)+  
  theme(axis.text.y.left = element_text(size = 12),                                        #adjust y-axis text
        axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 1, angle = 0),          #adjust x-axis label position
        axis.title.y = element_text(size = 20),
        strip.text = element_text(face = "bold", size = 20))  

##### ALPHA DIVERSITY STATS

alphadiv <- estimate_richness(pps_filtered, measures = c("Observed", "Shannon")) %>%
  rownames_to_column(var = "SampleID") %>%
  left_join(sdata3, by = "SampleID") 
head(alphadiv)

#Kruskal-Wallis Test
library(FSA)
set.seed(81)

##BY ENRICHMENT
#Observed
kruskal.test(Observed ~ Transfer, data = alphadiv) 

#Shannon
kruskal.test(Shannon ~ Transfer, data = alphadiv)

#Dunn test (post hoc)

##Observed
dunnO <- dunnTest(Observed ~ Transfer,
                  data=alphadiv,
                  method="bh")
dunnO

##Shannon
dunnS <- dunnTest(Shannon ~ Transfer,
                  data=alphadiv,
                  method="bh")
dunnS


#Making a PCoA plot
### UNIFRAC
library

#adding a phylogenetic tree to phyloseq object using ape library
random_tree = rtree(ntaxa(pps_filt), rooted=TRUE, tip.label=taxa_names(pps_filt))
plot(random_tree)

justbacteria3 = merge_phyloseq(pps_filt, samdata, random_tree)
justbacteria3

#ordination
distance <- ordinate(
  physeq = justbacteria3, 
  method = "PCoA", 
  distance = "unifrac"
)
#summary(distance)
#distance

#plot
colors <- c("lightblue", "maroon", "olivedrab")

PlotA <- plot_ordination(
  physeq = justbacteria3,                                                          #phyloseq object
  ordination = distance)+                                                #ordination
  geom_point(aes(fill = Substrate, shape = Replicate), size = 6) +                         #sets fill color to sampletype
  scale_fill_manual(values = colors) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25))+
  theme_linedraw() +                                                      #changes theme, removes grey background
  theme(                             
    legend.title = element_blank(),                                      #removes legend title
    #legend.background = element_rect(fill = "white", color = "black"),  #adds black boarder around legend
    legend.position = "bottom",
    legend.text = element_text(size = 20, face = "bold"),                                 
    axis.text.y.left = element_text(size = 10),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 20),
    strip.text = element_text(face = "bold", size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))          #fills legend points based on the fill command
PlotA

PlotB <- plot_ordination(
  physeq = justbacteria3,                                                          #phyloseq object
  ordination = distance)+                                                #ordination
  facet_grid(cols = vars(Transfer), rows = vars(Environment))+
  geom_point(aes(fill = Substrate, shape = Replicate), size = 6) +                         #sets fill color to sampletype
  scale_fill_manual(values = colors) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25))+
  theme_linedraw() +                                                      #changes theme, removes grey background
  theme(                             
    legend.title = element_blank(),                                      #removes legend title
    #legend.background = element_rect(fill = "white", color = "black"),  #adds black boarder around legend
    legend.position = "bottom",
    legend.text = element_text(size = 20, face = "bold"),                                 
    axis.text.y.left = element_text(size = 10),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 20),
    strip.text = element_text(face = "bold", size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))          #fills legend points based on the fill command
PlotB

ggarrange(PlotA, PlotB, labels = c("A", "B"), common.legend = TRUE, legend = "bottom")

#PAIRWISE PERMANOVA

library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
library(phyloseq)

# Adonis test

#Comparing substrate type (no inocula samples)
head(sdata3)

adonisData <- otu_table(pps_filtered) %>%                                      #use sequence table that has inocula removed
  as.data.frame() %>%
  rownames_to_column(var = "SampleID") %>%                                          #change rownames to a column so there is a common variable to join by
  left_join(sdata3, by = "SampleID") %>%                                            #join sample data to the sequence table
  select(-c("SampleID", "Environment", "Replicate", "Transfer", "Sample_Type")) %>% #remove all metadata columns except the one to be used to compare
  select(Substrate, everything())                                                   #rearrange so substrate column is first
dim(adonisData)
head(adonisData[,1:10])

#similarity euclidean from vegdist and holm correction
#pairwise.adonis(x=adonisData[,2:3722],factors=adonisData$Substrate,sim.function='vegdist',
#               sim.method='bray',p.adjust.m='holm')

pairwise.adonis2(adonisData[,2:3722]~Substrate,method="bray",data=adonisData)

#Comparing transfer
adonisData2 <- otu_table(pps_filtered) %>%                                      #use sequence table that has inocula removed
  as.data.frame() %>%
  rownames_to_column(var = "SampleID") %>%                                          #change rownames to a column so there is a common variable to join by
  left_join(sdata3, by = "SampleID") %>%                                            #join sample data to the sequence table
  select(-c("SampleID", "Environment", "Replicate", "Substrate", "Sample_Type")) %>% #remove all metadata columns except the one to be used to compare
  select(Transfer, everything())                                                   #rearrange so substrate column is first
dim(adonisData2)
head(adonisData2[,1:10])

pairwise.adonis2(adonisData2[,2:3722]~Transfer,method="bray",data=adonisData2)



#EXPLORING TAXA


#Summarize abundance of each class
genusabundance <- pps_filt %>%
  tax_glom(taxrank = "Species") %>%                      # agglomerate at class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Species) 
head(genusabundance)

#save color palatte
colors10 <- c(
  "black",   "darkcyan",     "orchid1",   "green",       "blue",   "olivedrab",
  "grey47",  "cyan",    "coral1",     "yellow",    "darkgreen",   "palegoldenrod",    
  "grey77",  "darkblue",     "orange",    "red",         "mediumpurple1", "tan4",   "purple4",
  "dodgerblue",    "white", "firebrick", "yellowgreen", "magenta", "blue", "green", "red", "orchid", "lightblue"
) 

#Select and summarize necessary variables
head(genusabundance)

all <- genusabundance %>%
  select(Phylum, Class, Family, Genus, Species, Sample, Abundance, Environment, Environment_Label, Substrate, Substrate_Label, Sample_Type, Transfer, Replicate) %>%
  filter(Abundance > 0) %>%
  mutate(
    Phylum = as.character(Phylum),
    Class = as.character(Class),
    Family = as.character(Family),
    Genus = as.character(Genus),
    Replicate2 = case_when(               #changing R2A to R2 and R2B to R3 for S4TA T2-T4
      Sample == "S4TAR2AT2" ~ "R1",
      Sample == "S4TAR2AT3" ~ "R1",
      Sample == "S4TAR2AT4" ~ "R1",
      Sample == "S4TAR2BT2" ~ "R2",
      Sample == "S4TAR2BT3" ~ "R2",
      Sample == "S4TAR2BT4" ~ "R2"
    ),
    Replicate = ifelse(is.na(Replicate2), Replicate, Replicate2)
  )
head(all)
#View(all)

class <- all %>%
  group_by(Sample_Type, Substrate, Transfer, Phylum, Class)%>%
  summarise(Abundance = sum(Abundance)/n()) %>%
#  mutate(Class = ifelse(Abundance < 0.001, "<0.1%", Class),
 #        Phylum = factor(Phylum, c("Acidobacteriota",   "Actinobacteriota",  "Bacteroidota",      "Bdellovibrionota",  "Deinococcota",     
  #                                 "Desulfobacterota",  "Firmicutes",        "Gemmatimonadota",   "Myxococcota",       "Nitrospirota",    
   #                                "Patescibacteria",   "Planctomycetota",   "Proteobacteria",    "Verrucomicrobiota")),
    #     Class = factor(Class, c("<0.1%",               "Blastocatellia",      "Vicinamibacteria",    "Acidimicrobiia",      "Actinobacteria",      "Thermoleophilia",    
     #                            "Bacteroidia",         "Rhodothermia",        "Bdellovibrionia",     "Deinococci",          "Syntrophia",         
      #                           "Bacilli",             "Clostridia",          "Gemmatimonadetes",    "Polyangia",          
       #                          "Nitrospiria",         "Saccharimonadia",     "Phycisphaerae",       "Planctomycetes",      "Alphaproteobacteria",
        #                         "Gammaproteobacteria", "Verrucomicrobiae",    "Acidobacteriae",      "Negativicutes",       "Myxococcia" ))) %>%
  filter(Transfer == "T4")
head(class)
unique(class$Phylum)
unique(class$Class)



#MAKING A TAXA PLOT 

#save color palatte

phylum_colors <- c(
  "grey22", "darkcyan", "orchid1", "green", "orange", "blue", "tomato2", "olivedrab", "grey47",
  "cyan", "coral3", "darkgreen", "magenta", "palegoldenrod", "dodgerblue", "firebrick", "yellow", "purple4",
  "lightblue", "grey77", "mediumpurple1", "tan4", "red", "darkblue", "yellowgreen")



class_colors <- c(
  "black",      "darkcyan",     "orchid1",          "blue",   "olivedrab",
  "grey47",     "cyan",         "coral1",      "darkgreen",   "palegoldenrod",    
  "dodgerblue", "firebrick",    "yellowgreen", "magenta",     "blue2", "red3", "yellow", "lightblue",
  "grey77",     "darkblue",     "orange",      "tomato2",         "mediumpurple1", "tan4",   "purple4")  



cla <- ggplot(class)+
  geom_col(mapping = aes(x = Substrate, y = Abundance, fill = Class), color = "black", position = "fill", show.legend = TRUE)+
 # facet_grid(cols = vars(Transfer), rows = vars(Environment))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = class_colors) +
  xlab(NULL)+
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 1, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 20),
        legend.position = "right",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text = element_text(size = 18, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill=guide_legend(ncol=1,byrow=TRUE))
cla


genus <- all %>%
  filter(Abundance > 0) %>%
  filter(Transfer == "T4") %>%
  #filter(Substrate != "Control") %>%
  dplyr::group_by(Environment, Substrate, Replicate, Phylum, Class, Family, Genus)%>%
  summarise(sumAbundance = sum(Abundance),
            count = n(),
            avg_abundance = sumAbundance/count) %>%
            #avgAbundance = Abundance/n()) %>%
  mutate(Genus = ifelse(avg_abundance < 0.1, "<10%", Genus))
head(genus)
#View(genus)

colFunc <- colorRampPalette(c("purple4", "lightpink", "firebrick", "orange", "lemonchiffon", "olivedrab4", "darkcyan", "lightblue", "darkblue"))
color_list <- colFunc(length(unique(genus$Genus)))
colors <- c("black", color_list)

gen <- ggplot(genus)+
  facet_grid(cols = vars(Substrate, Environment))+
  geom_col(mapping = aes(x = Replicate, y = avg_abundance, fill = Genus), color = "black", position = "stack", show.legend = TRUE)+
  # facet_grid(cols = vars(Transfer), rows = vars(Environment))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors) +
  xlab(NULL)+
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 1, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 20),
        legend.position = "right",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text = element_text(size = 18, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill=guide_legend(ncol=1,byrow=TRUE))
gen


#TA only
genusTA <- filter(genus, Substrate == "TA")
head(genusTA)

colFunc <- colorRampPalette(c("purple4", "lightpink", "firebrick", "orange", "lemonchiffon", "olivedrab4", "darkcyan", "lightblue", "darkblue"))
color_list <- colFunc(length(unique(genusTA$Genus)))
colors <- c("black", color_list)

ggplot(genusTA)+
  facet_grid(cols = vars(Substrate, Environment))+
  geom_col(mapping = aes(x = Replicate, y = avg_abundance, fill = Genus), color = "black", position = "stack", show.legend = TRUE)+
  # facet_grid(cols = vars(Transfer), rows = vars(Environment))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors) +
  xlab(NULL)+
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 1, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 20),
        legend.position = "right",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text = element_text(size = 18, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill=guide_legend(ncol=1,byrow=TRUE))



#TPA only
genusTPA <- filter(genus, Substrate == "TPA")

colFunc <- colorRampPalette(c("purple4", "lightpink", "firebrick", "orange", "lemonchiffon", "olivedrab4", "darkcyan", "lightblue", "darkblue"))
color_list <- colFunc(length(unique(genusTPA$Genus)))
colors <- c("black", color_list)

ggplot(genusTPA)+
  facet_grid(cols = vars(Substrate, Environment))+
  geom_col(mapping = aes(x = Replicate, y = avg_abundance, fill = Genus), color = "black", position = "stack", show.legend = TRUE)+
  # facet_grid(cols = vars(Transfer), rows = vars(Environment))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors) +
  xlab(NULL)+
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 1, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 20),
        legend.position = "right",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text = element_text(size = 18, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill=guide_legend(ncol=1,byrow=TRUE))


#Control only
genusCont <- filter(genus, Substrate == "Control")

colFunc <- colorRampPalette(c("purple4", "lightpink", "firebrick", "orange", "lemonchiffon", "olivedrab4", "darkcyan", "lightblue", "darkblue"))
color_list <- colFunc(length(unique(genusCont$Genus)))
colors <- c("black", color_list)

ggplot(genusCont)+
  facet_grid(cols = vars(Substrate, Environment))+
  geom_col(mapping = aes(x = Replicate, y = avg_abundance, fill = Genus), color = "black", position = "stack", show.legend = TRUE)+
  # facet_grid(cols = vars(Transfer), rows = vars(Environment))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors) +
  xlab(NULL)+
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 1, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 20),
        legend.position = "right",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text = element_text(size = 18, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill=guide_legend(ncol=1,byrow=TRUE))



#Looking at genera that are present in TA or TPA treatments @ > 10% relative abundance
#unique(genusCont$Genus)

#genera_list <- c(unique(genusTA$Genus), unique(genusTPA$Genus))
#genera_list <- genera_list[! genera_list %in% c("<10%")]
#genera_list

#View(all)

# same figure filtered even more (genus_list from phyloseq script)
genera_list <- c("Arthrobacter",       "Bradyrhizobium",     "Niabella",           #"Eoetvoesia",        
                "Hydrogenophaga",     "Ramlibacter",        "Variovorax",         "Pseudomonas",       
                "Hoeflea",            "Ellin6055",          "Pseudoxanthomonas",  #"Terrimonas",        
                "Taibaiella",         "UTBCD1",             "Rhodococcus",        "Bordetella",        
                "Noviherbaspirillum", "Micromonospora",     #"Pseudorhodoplanes",  
                "Thermomonas",       
                "Acidovorax",         "Paenibacillus",      "Mesorhizobium",      #"Luteimonas",        
                "Rhodococcus",        "Leadbetterella",     "SWB02",              "Mesorhizobium",     
                "Bordetella",         #"Bosea",              
                "Aliihoeflea",        "Verticiella",       
                "Aminobacter",        "Pseudomonas",        "Sphingobacterium",   "Hydrogenophaga",    
                "Flavobacterium",     "Achromobacter",      "Thermomonas",        "Bacillus",          
                "Micromonospora",     "Pseudoxanthomonas")

genus_subset <- all %>%
  filter(Abundance > 0) %>%
  mutate(Genus2 = ifelse(is.na(Genus), "", Genus),
         Species2 = ifelse(is.na(Species), "", Species)) %>%
  unite(Taxonomy, Genus2, Species2, sep = " ") %>%
  mutate(Taxonomy = ifelse(Taxonomy == " ", "Unclassified Organism", Taxonomy)) %>%
  dplyr::group_by(Environment, Environment_Label, Substrate, Substrate_Label, Transfer, Replicate, Phylum, Class, Family, Genus)%>%
  summarise(sumAbundance = sum(Abundance),
            count = n(),
            avg_abundance = sumAbundance/count) %>%
  filter(Genus %in% genera_list) %>%
  mutate(Treatment = case_when(
    Substrate == "Control" ~ "CT",
    Substrate == "TPA" ~ "TP", 
    Substrate == "TA" ~ "TA"
  ))
head(genus_subset)

#finding mean relative abundance and sd of select genera
genus_subset_filt <- all %>%
  ungroup() %>%
  filter(Genus == "Rhodococcus" | Genus == "Flavobacterium" | Genus == "Hydrogenophaga"|
           Genus == "Arthrobacter" | Genus == "Bradyrhizobium" | Genus == "Variovorax") %>%
  group_by(Genus, Substrate) %>%
  summarise(
    avg_abundance = sum(Abundance)/n(),
    sd_abundance = sd(Abundance)
  )
View(genus_subset_filt)

colFunc <- colorRampPalette(c("purple4", "lightpink", "maroon", "firebrick", "orange", "lemonchiffon", "olivedrab4", "darkgreen", "cyan", "lightblue", "darkblue"))
color_list <- colFunc(length(unique(genus_subset$Genus)))
colors <- c(#"black", 
            color_list)

#barplot of taxa that are present at > 10% in final transfer TA and TPA treatments
ggplot(genus_subset)+
  facet_nested(cols = vars(Environment_Label, Transfer), rows = vars(Substrate_Label))+
  geom_col(mapping = aes(x = Replicate, y = avg_abundance, fill = Genus), color = "black", position = "stack", show.legend = TRUE)+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors) +
  xlab(NULL)+
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.x = element_blank(),#size = 20, angle = 90, vjust = 1, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text = element_text(size = 18, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill=guide_legend(ncol=5,byrow=FALSE))

#making a heatmap

colors <- c("purple4", "chocolate", "orange")
max <- round(max(genus_subset$avg_abundance), digits = 1)
min <- round(min(genus_subset$avg_abundance), digits = 1)

ggplot(genus_subset)+
  facet_nested(cols = vars(Transfer, Substrate))+
  geom_tile(mapping = aes(x = Environment, y = Genus, fill = avg_abundance), color = "black", show.legend = TRUE)+
  scale_fill_gradientn(colors = colors, breaks = seq(min, max, by = 0.1),#)+#,
                       limits = c(min, max), labels = as.character(seq(min, max, by = 0.1))) +
  ylab(NULL)+
  xlab(NULL)+
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size = 13, angle = 90, vjust = 0.5, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        #legend.spacing.x = unit(0.1, 'mm'),
        #legend.spacing.y = unit(0.05, 'mm'),
        #plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text.x = element_text(size = 18, face = "bold", angle = 0),
        #strip.text.y = element_text(size = 18, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill = guide_colourbar(barwidth = 40, barheight = 1, title = "Relative Aundance"))

#Interesting genera - stats and boxplot
genus_subset2 <- genus_subset #filter(genus_subset, Transfer == "T4")
head(genus_subset2)

unique(genus_subset2$Genus)

Arthrobacter <- filter(genus_subset2, Genus == "Arthrobacter")
Rhodococcus <- filter(genus_subset2, Genus == "Rhodococcus")
Aminobacter <- filter(genus_subset2, Genus == "Aminobacter")
Mesorhizobium <- filter(genus_subset2, Genus == "Mesorhizobium")
Bradyrhizobium <- filter(genus_subset2, Genus == "Bradyrhizobium")
Acidovorax <- filter(genus_subset2, Genus == "Acidovorax")
Pseudomonas <- filter(genus_subset2, Genus == "Pseudomonas")
Pseudoxanthomonas <- filter(genus_subset2, Genus == "Pseudoxanthomonas")
Hydrogenophaga <- filter(genus_subset2, Genus == "Hydrogenophaga")
Variovorax <- filter(genus_subset2, Genus == "Variovorax")
Bordetella <- filter(genus_subset2, Genus == "Bordetella")
Flavobacterium <- filter(genus_subset2, Genus == "Flavobacterium")
Paenibacillus <- filter(genus_subset2, Genus == "Paenibacillus")
Sphingobacterium <- filter(genus_subset2, Genus == "Sphingobacterium")
Taibaiella <- filter(genus_subset2, Genus == "Taibaiella")
Micromonospora <- filter(genus_subset2, Genus == "Micromonospora")
Bacillus <- filter(genus_subset2, Genus == "Bacillus")

#Kruskal-Wallis test
kruskal.test(avg_abundance ~ Substrate, data = Arthrobacter)
kruskal.test(avg_abundance ~ Substrate, data = Rhodococcus)
kruskal.test(avg_abundance ~ Substrate, data = Aminobacter)
kruskal.test(avg_abundance ~ Substrate, data = Mesorhizobium)
kruskal.test(avg_abundance ~ Substrate, data = Bradyrhizobium)
kruskal.test(avg_abundance ~ Substrate, data = Acidovorax)
kruskal.test(avg_abundance ~ Substrate, data = Pseudomonas)
kruskal.test(avg_abundance ~ Substrate, data = Pseudoxanthomonas)
kruskal.test(avg_abundance ~ Substrate, data = Hydrogenophaga)
kruskal.test(avg_abundance ~ Substrate, data = Variovorax)
kruskal.test(avg_abundance ~ Substrate, data = Bordetella)
kruskal.test(avg_abundance ~ Substrate, data = Flavobacterium)
kruskal.test(avg_abundance ~ Substrate, data = Paenibacillus)
kruskal.test(avg_abundance ~ Substrate, data = Sphingobacterium)
kruskal.test(avg_abundance ~ Substrate, data = Taibaiella)
kruskal.test(avg_abundance ~ Substrate, data = Micromonospora)
kruskal.test(avg_abundance ~ Substrate, data = Bacillus)


#Dunn test (post hoc) for significant genera
rhodo_res <- dunnTest(avg_abundance ~ Substrate, data=Rhodococcus, method="bh")
rhodo_res

acido_res <- dunnTest(avg_abundance ~ Substrate, data=Acidovorax, method="bh")
acido_res

pseudx_res <- dunnTest(avg_abundance ~ Substrate, data=Pseudoxanthomonas, method="bh")
pseudx_res

vario_res <- dunnTest(avg_abundance ~ Substrate, data=Variovorax, method="bh")
vario_res

paeni_res <- dunnTest(avg_abundance ~ Substrate, data=Paenibacillus, method="bh")
paeni_res

#Making a boxplot

colors <- c("lightblue", "maroon", "olivedrab")

ggplot(genus_subset2)+
  facet_wrap(~Genus, nrow = 3)+#, rows = vars(Environment))+
  geom_boxplot(mapping = aes(x = Substrate, y = avg_abundance, fill = Substrate), color = "black", show.legend = TRUE)+
  ylab("Relative Abundance") +
  scale_fill_manual(values = colors) +
  xlab(NULL)+
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 22, color = "black"),
        axis.title.y = element_text(size = 25, color = "black"),
        axis.text.x = element_blank(),#size = 20, angle = 0, vjust = 1, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 25),
        legend.position = c(0.9, 0.2),
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text = element_text(size = 22, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 30))+
  guides(fill=guide_legend(nrow=3,byrow=TRUE))

#??theme

genus_subset3 <- genus_subset %>%
  filter(
    Species == "facilis" |
      Species == "hinzii" |
      Species == "elkanii" |
      Species == "japonicum" |
      Species == "koleovorans" |
      Species == "anguilliseptica" |
      Species == "yamanorum" |
      Species == "mexicana" | 
      Species == "erythropolis" |
      Species == "opacus" | 
      Species == "rhodochrous" |
      Species == "ruber" |
      Species == "boronicumulans"
  ) %>%
  filter(Transfer == "T4")
head(genus_subset3)


ggplot(genus_subset3)+
  facet_wrap(~Taxonomy, nrow = 2)+#, rows = vars(Environment))+
  geom_boxplot(mapping = aes(x = Substrate, y = avg_abundance, fill = Substrate), color = "black", show.legend = TRUE)+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors) +
  xlab(NULL)+
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.x = element_blank(),#size = 20, angle = 0, vjust = 1, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text = element_text(size = 10, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill=guide_legend(ncol=3,byrow=TRUE))


colors <- c("purple4", "firebrick", "orange", "lightgoldenrod", "olivedrab4", "darkgreen", "blue", "cyan", "lightblue")

ggplot(genus_subset3)+
  facet_nested(cols = vars(Environment, Substrate))+
  geom_col(mapping = aes(x = Transfer, y = avg_abundance, fill = Genus), color = "black", position = "stack", show.legend = TRUE)+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors) +
  xlab(NULL)+
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 1, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text = element_text(size = 18, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill=guide_legend(ncol=6,byrow=TRUE))




#REMOVE TAXA PRESENT IN CONTROLS FROM ALL SAMPLES
#devtools::install_github("donaldtmcknight/microDecon") #Installs microDecon
library(microDecon)

example <- cbind.data.frame(c("OTU1","OTU2","OTU3","OTU4","OTU5","OTU6"),
                            c(0,200,1000,50,0,25),
                            c(0,220,800,30,0,10),
                            c(0,180,1300,70,0,30),
                            c(60,660,1440,70,2400,30),
                            c(64,520,1000,48,1900,20),
                            c(40,480,700,35,2100,15),
                            c("K_Bacteria; P_Actinobacteria","K_Bacteria; P_Proteobacteria","K_Bacteria; P_Proteobacteria","K_Bacteria; P_Bacteroidetes","K_Bacteria","K_Bacteria"))
colnames(example) <- c("OTU_ID","Blank1","Blank2","Blank3","Pop1_Sample1","Pop1_Sample2","Pop2_Sample3","Taxa")
head(example)

decontaminated <- decon(data = example,numb.blanks=3,numb.ind=c(2,1),taxa=T)
head(decontaminated)



#Reorganize metadata so that the control samples from transfer 4 are listed first
#the following code formats a non-normalized phyloseq object for microDecon

sampledata <- as_tibble(sample_data(ps)) 

sdata <- sampledata %>%
  mutate(Substrate = factor(Substrate, levels = c("Control", "TPA", "TA", "Inocula")),
         Transfer = factor(Transfer, levels = c("T4", "T3", "T2", "T1"))) %>%
  group_by(Substrate, Transfer) %>%
  arrange(Substrate, .by_group = TRUE)
head(sdata)

#Reorder the sequence table to match the new order of the metadata, then join Genus classification
stable <- as.data.frame(otu_table(ps)) %>%
  rownames_to_column(var = "SampleID")
stable[1:4,1:4]

ttable <- as.data.frame(tax_table(ps)) %>%
  select(Genus) %>%
  rownames_to_column(var = "OTU_ID")
head(ttable)

stable2 <- sdata %>%
  ungroup() %>%
  select(SampleID) %>%
  left_join(stable) %>%
  column_to_rownames(var = "SampleID") %>%
  t() %>%                                        #transpose the sequence table to match the needed format for microDecon
  as.data.frame() %>%
  rownames_to_column(var = "OTU_ID") %>%
  left_join(ttable)
stable2[1:10,160:169]

#Make sure dimensions are acceptable
dim(sdata)
dim(stable)
dim(ttable)
dim(stable2)

#Will be removing all reads present in control samples from  transfer 4
#These lines of code determine how many samples are each of the categories to be used in numb.ind
sum(sampledata$Substrate == "Control" & sampledata$Transfer == "T4")     #15 "blank" samples
sum(sampledata$Substrate == "Control" & sampledata$Transfer != "T4")     #45 samples
sum(sampledata$Substrate == "TPA")                                       #60 samples
sum(sampledata$Substrate == "TA")                                        #58 samples
sum(sampledata$Substrate == "Inocula")                                   #14 samples

decontaminated <- decon(data = stable2,numb.blanks=15,numb.ind=c(45, 60, 58, 14),taxa=T)
summary(decontaminated)
decontaminated$decon.table  #this is the decontaminated otu data

decon <- decontaminated$decon.table
rownames(decon) <- NULL
decon[1:4,1:4]

#Now reformat to make a new decontaminated phyloseq object

decon_sdata <- sdata %>%
  filter(Substrate != "Control" | Transfer != "T4") %>%
  mutate(Sample_ID = SampleID) %>%
  column_to_rownames(var = "Sample_ID")
head(decon_sdata)

decon_stable <- as_tibble(decon) %>%
  select(-c("Mean.blank", "Genus")) %>%
  column_to_rownames(var = "OTU_ID") %>%
  t()
decon_stable[1:4,1:4]

original_ttable <- as.data.frame(tax_table(ps)) %>%
  rownames_to_column(var = "OTU_ID")
head(original_ttable)

decon_ttable <- decon %>%
  select(OTU_ID, Genus) %>%
  left_join(original_ttable) %>%
  select(OTU_ID, Kingdom, Phylum, Class, Order, Family, Genus) %>%
  column_to_rownames(var = "OTU_ID")
head(decon_ttable)


#make phyloseq object
samdata = sample_data(decon_sdata)
seqtab = otu_table(decon_stable, taxa_are_rows = FALSE)
taxtab = tax_table(decon_ttable)
colnames(taxtab) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
head(taxtab)

sample_names(samdata)
sample_names(seqtab)

taxa_names(seqtab)
taxa_names(taxtab)

#combine all components into a phyloseq object
ps_decon = phyloseq(otu_table(seqtab), tax_table(taxtab), sample_data(samdata))
ps_decon


#Filter out eukaryotes and mitochondria
ps_decon2 <- ps_decon %>%
  subset_taxa(
    Kingdom == "Bacteria" &                   #only bacteria
      Family  != "Mitochondria" &             #filter out mitochondria
      Class   != "Chloroplast"                #filter out chloroplasts
  )
ps_decon2

### Section 3: Normalize the data

samplesover1000_all <- subset_samples(ps_decon2, sample_sums(ps_decon2) > 1000)
any(taxa_sums(samplesover1000_all) == 0)
sum(taxa_sums(samplesover1000_all) == 0)
prune_samplesover1000_all <- prune_taxa(taxa_sums(samplesover1000_all) > 0, samplesover1000_all)
any(taxa_sums(prune_samplesover1000_all) == 0)
set.seed(81). #for reproducible data
rarefy_samplesover1000_all <- rarefy_even_depth(prune_samplesover1000_all, rngseed= 81, sample.size = min(sample_sums(prune_samplesover1000_all)))
psd<- rarefy_samplesover1000_all
psd

###Section 4: plotting

#sample colors
sample_colors <- c("Inocula" = "blue", "Control" = "lightblue", "TPA" = "maroon", "TA" = "olivedrab")
transfer_colors <- c("I" = "blue4", "T1" = "blue3", "T2" = "blue2", "T3" = "blue1", "T4" = "lightblue")

sample_types <- c("Inocula", "TPA", "TA", "Control")

#violin plot
psd %>%                                                              #phyloseq object
  plot_richness(
    x = "Substrate",                                                 #compare diversity of datatype
    measures = c("Observed", "Shannon")) +                           #choose diversity measures
  geom_violin(aes(fill = Substrate), show.legend = FALSE)+           #make violin plot, set fill aes to sampletype
  geom_boxplot(width=0.1) +                                          #add boxplot, set width
  theme_classic()+                                                   #change theme to classic
  xlab(NULL)+                                                        #no label on x-axis
  theme(axis.text.y.left = element_text(size = 20),                  #adjust y-axis text
        axis.text.x = element_text(size = 20, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 20))+                     #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 25))+        #adjust headings
  scale_fill_manual(values = sample_colors)+                         #set fill colors
  scale_x_discrete(                                                  #change x-axis labels
    limits = sample_types)+                   
  theme(plot.title=element_text(size = 25, face = "bold", hjust = 0.5)) #change title size, face and position


psd <- subset_samples(psd, Substrate != "Blank")
psd_filt <- subset_samples(psd, Substrate != "Inocula")


#Making a PCoA plot
### UNIFRAC
library(ape)

#only include the final time point
deconT4 <- subset_samples(psd_filt, Transfer == "T4")

#adding a phylogenetic tree to phyloseq object using ape library
random_tree = rtree(ntaxa(deconT4), rooted=TRUE, tip.label=taxa_names(deconT4))
plot(random_tree)

justbacteria4 = merge_phyloseq(deconT4, samdata, random_tree)
justbacteria4

#ordination
distance <- ordinate(
  physeq = justbacteria4, 
  method = "PCoA", 
  distance = "unifrac"
)
#summary(distance)
#distance

#plot

colors <- c("lightblue", "maroon", "olivedrab")

PlotA <- plot_ordination(
  physeq = justbacteria4,                                                          #phyloseq object
  ordination = distance)+                                                #ordination
  geom_point(aes(fill = Substrate, shape = Replicate), size = 6) +                         #sets fill color to sampletype
  scale_fill_manual(values = colors) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25))+
  theme_linedraw() +                                                      #changes theme, removes grey background
  theme(                             
    legend.title = element_blank(),                                      #removes legend title
    #legend.background = element_rect(fill = "white", color = "black"),  #adds black boarder around legend
    legend.position = "bottom",
    legend.text = element_text(size = 20, face = "bold"),                                 
    axis.text.y.left = element_text(size = 10),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 20),
    strip.text = element_text(face = "bold", size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))          #fills legend points based on the fill command
PlotA

#EXPLORING TAXA IN DECONTAMINATED DATA

#Summarize abundance of each class
genusabundance <- psd_filt %>%
  tax_glom(taxrank = "Genus") %>%                      # agglomerate at class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Genus) 
head(genusabundance)

#save color palatte
colors10 <- c(
  "black",   "darkcyan",     "orchid1",   "green",       "blue",   "olivedrab",
  "grey47",  "cyan",    "coral1",     "yellow",    "darkgreen",   "palegoldenrod",    
  "grey77",  "darkblue",     "orange",    "red",         "mediumpurple1", "tan4",   "purple4",
  "dodgerblue",    "white", "firebrick", "yellowgreen", "magenta", "blue", "green", "red", "orchid", "lightblue"
) 

#Select and summarize necessary variables
head(genusabundance)

all <- genusabundance %>%
  select(Phylum, Class, Family, Genus, Sample, Abundance, Environment, Substrate, Sample_Type, Transfer, Replicate) %>%
  filter(Abundance != 0) %>%
  mutate(
    Phylum = as.character(Phylum),
    Class = as.character(Class),
    Family = as.character(Family),
    Genus = as.character(Genus)
  )
head(all)

phylum <- all %>%
  dplyr::group_by(Environment, Sample_Type, Substrate, Transfer, Phylum)%>%
  summarise(Abundance = sum(Abundance)/n()) %>%
  mutate(Phylum.1p = ifelse(Abundance < 0.01, "<1%", Phylum))
head(phylum)

range(phylum$Abundance)

class <- all %>%
  group_by(Sample_Type, Substrate, Transfer, Phylum, Class)%>%
  summarise(Abundance = sum(Abundance)/n()) %>%
  mutate(Class = ifelse(Abundance < 0.01, "<1%", Class)) %>%
#  mutate(Class = ifelse(Abundance < 0.001, "<0.1%", Class),
 #        Phylum = factor(Phylum, c("Acidobacteriota",   "Actinobacteriota",  "Bacteroidota",      "Bdellovibrionota",  "Deinococcota",     
  #                                 "Desulfobacterota",  "Firmicutes",        "Gemmatimonadota",   "Myxococcota",       "Nitrospirota",    
   #                                "Patescibacteria",   "Planctomycetota",   "Proteobacteria",    "Verrucomicrobiota")),
    #     Class = factor(Class, c("<0.1%",               "Blastocatellia",      "Vicinamibacteria",    "Acidimicrobiia",      "Actinobacteria",      "Thermoleophilia",    
     #                            "Bacteroidia",         "Rhodothermia",        "Bdellovibrionia",     "Deinococci",          "Syntrophia",         
      #                           "Bacilli",             "Clostridia",          "Gemmatimonadetes",    "Polyangia",          
       #                          "Nitrospiria",         "Saccharimonadia",     "Phycisphaerae",       "Planctomycetes",      "Alphaproteobacteria",
        #                         "Gammaproteobacteria", "Verrucomicrobiae",    "Acidobacteriae",      "Negativicutes",       "Myxococcia" ))) %>%
  filter(Transfer == "T4")
head(class)
unique(class$Phylum)
unique(class$Class)

#MAKING A TAXA PLOT 

#save color palatte

phylum_colors <- c(
  "grey22", "darkcyan", "orchid1", "green", "orange", "blue", "tomato2", "olivedrab", "grey47",
  "cyan", "coral3", "darkgreen", "magenta", "palegoldenrod", "dodgerblue", "firebrick", "yellow", "purple4",
  "lightblue", "grey77", "mediumpurple1", "tan4", "red", "darkblue", "yellowgreen")



class_colors <- c(
  "black",      "darkcyan",     "orchid1",          "blue",   "olivedrab",
  "grey47",     "cyan",         "coral1",      "darkgreen",   "palegoldenrod",    
  "dodgerblue", "firebrick",    "yellowgreen", "magenta",     "blue2", "red3", "yellow", "lightblue",
  "grey77",     "darkblue",     "orange",      "tomato2",         "mediumpurple1", "tan4",   "purple4")  


length(phylum_colors)
length(unique(phylum$Phylum))
length(class_colors)
length(unique(class$Class))

#plot
head(phylum)

phy <- ggplot(phylum)+
  geom_col(mapping = aes(x = Substrate, y = Abundance, fill = Phylum), color = "black", position = "fill", show.legend = TRUE)+
  facet_grid(cols = vars(Transfer), rows = vars(Environment))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = phylum_colors) +
  xlab(NULL)+
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 18, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 18, angle = 90, hjust = 1, vjust = 0.5, color = "black"),
        legend.text = element_text(size = 20),
        legend.position = "right",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text = element_text(size = 15, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill=guide_legend(ncol=1,byrow=TRUE))
phy

View(class)

cla <- ggplot(class)+
  geom_col(mapping = aes(x = Substrate, y = Abundance, fill = Class), color = "black", position = "fill", show.legend = TRUE)+
  # facet_grid(cols = vars(Transfer), rows = vars(Environment))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = genus_colors) +
  xlab(NULL)+
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 1, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 20),
        legend.position = "right",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text = element_text(size = 18, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill=guide_legend(ncol=1,byrow=TRUE))
cla

#Plotting Genus

#genus_colors <- c(
 # "black",      "darkcyan",     "orchid1",          "blue",   "olivedrab",
  #"grey47",     "cyan",         "coral1",      "darkgreen",   "palegoldenrod",    
#  "dodgerblue", "firebrick",    "yellowgreen", "magenta",     "blue2", "red3", "yellow", "lightblue",
 # "grey77",     "darkblue",     "orange",      "tomato2",         "mediumpurple1", "tan4",   "purple4")  

genus <- all %>%
  group_by(Sample_Type, Substrate, Transfer, Phylum, Class, Family, Genus)%>%
  summarise(Abundance = sum(Abundance)/n()) %>%
  mutate(Genus = ifelse(Abundance < 0.1, "<10%", Genus)) %>%
  filter(Transfer == "T4")
head(genus)
#View(genus)

colfunc <- colorRampPalette(c("lemonchiffon", "orange", "green", "darkgreen", "darkcyan", "cyan",
                              "blue", "lightblue", "violet", "firebrick", "purple4"))
color_list <- colfunc(32)
genus_colors <- c("black", color_list)
length(genus_colors)

gen <- ggplot(genus)+
  geom_col(mapping = aes(x = Substrate, y = Abundance, fill = Genus), color = "black", position = "fill", show.legend = TRUE)+
  ylab("Proportion of Community") +
  scale_fill_manual(values = genus_colors) +
  xlab(NULL)+
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 1, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 20),
        legend.position = "right",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text = element_text(size = 18, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill=guide_legend(ncol=2,byrow=FALSE))
gen

#LOOKING AT CORE ASVs

#COMPOST
psT4_comp <- subset_samples(ps, Transfer == "T4" & Environment =="Compost")
psT4_comp

percent <- 0.8

compCore <- psT4_comp %>%
  filter_taxa(function(x) sum(x > 0) > (percent*length(x)), TRUE)  #change the decimal to desired percentage
compCore

tax_table(compCore)

compTPA <- psT4_comp %>%
  subset_samples(Substrate == "TPA") %>%
  filter_taxa(function(x) sum(x > 0) > (percent*length(x)), TRUE)  #change the decimal to desired percentage
compTPA

tax_table(compTPA)

compTA <- psT4_comp %>%
  subset_samples(Substrate == "TA") %>%
  filter_taxa(function(x) sum(x > 0) > (percent*length(x)), TRUE)  #change the decimal to desired percentage
compTA

tax_table(compTA)

compCont <- psT4 %>%
  subset_samples(Substrate == "Control") %>%
  filter_taxa(function(x) sum(x > 0) > (percent*length(x)), TRUE)  #change the decimal to desired percentage
compCont

tax_table(compCont)

#SEDIMENT 1
psT4_s1 <- subset_samples(ps, Transfer == "T4" & Environment =="Sediment1")
psT4_s1

percent <- 0.8

s1Core <- psT4_s1 %>%
  filter_taxa(function(x) sum(x > 0) > (percent*length(x)), TRUE)  #change the decimal to desired percentage
#s1Core
#tax_table(s1Core)

s1TPA <- psT4_s1 %>%
  subset_samples(Substrate == "TPA") %>%
  filter_taxa(function(x) sum(x > 0) > (percent*length(x)), TRUE)  #change the decimal to desired percentage
#s1TPA
#tax_table(s1TPA)

s1TA <- psT4_s1 %>%
  subset_samples(Substrate == "TA") %>%
  filter_taxa(function(x) sum(x > 0) > (percent*length(x)), TRUE)  #change the decimal to desired percentage
s1TA

tax_table(s1TA)

s1Cont <- psT4_s1 %>%
  subset_samples(Substrate == "Control") %>%
  filter_taxa(function(x) sum(x > 0) > (percent*length(x)), TRUE)  #change the decimal to desired percentage
s1Cont

#SEDIMENT 3
psT4_s3 <- subset_samples(ps, Transfer == "T4" & Environment =="Sediment3")
psT4_s3

percent <- 0.8

s3Core <- psT4_s3 %>%
  filter_taxa(function(x) sum(x > 0) > (percent*length(x)), TRUE)  #change the decimal to desired percentage
s3Core

tax_table(s3Core)

s3TPA <- psT4_s3 %>%
  subset_samples(Substrate == "TPA") %>%
  filter_taxa(function(x) sum(x > 0) > (percent*length(x)), TRUE)  #change the decimal to desired percentage
s3TPA

tax_table(s3TPA)

s3TA <- psT4_s3 %>%
  subset_samples(Substrate == "TA") %>%
  filter_taxa(function(x) sum(x > 0) > (percent*length(x)), TRUE)  #change the decimal to desired percentage
s3TA

tax_table(coreTA)

s3Cont <- psT4_s3 %>%
  subset_samples(Substrate == "Control") %>%
  filter_taxa(function(x) sum(x > 0) > (percent*length(x)), TRUE)  #change the decimal to desired percentage
s3Cont

tax_table(s3Cont)


#SEDIMENT 4
psT4_s4 <- subset_samples(ps, Transfer == "T4" & Environment =="Sediment4")
psT4_s4

percent <- 0.8

s4Core <- psT4_s4 %>%
  filter_taxa(function(x) sum(x > 0) > (percent*length(x)), TRUE)  #change the decimal to desired percentage
#s4Core

s4TPA <- psT4_s4 %>%
  subset_samples(Substrate == "TPA") %>%
  filter_taxa(function(x) sum(x > 0) > (percent*length(x)), TRUE)  #change the decimal to desired percentage
s4TPA

tax_table(s4TPA)

s4TA <- psT4_s4 %>%
  subset_samples(Substrate == "TA") %>%
  filter_taxa(function(x) sum(x > 0) > (percent*length(x)), TRUE)  #change the decimal to desired percentage
s4TA

tax_table(s4TA)

s4Cont <- psT4_s4 %>%
  subset_samples(Substrate == "Control") %>%
  filter_taxa(function(x) sum(x > 0) > (percent*length(x)), TRUE)  #change the decimal to desired percentage
s4Cont

tax_table(s4Cont)


#SEDIMENT 6
psT4_s6 <- subset_samples(ps, Transfer == "T4" & Environment =="Sediment6")
psT4_s6

percent <- 0.8

s6Core <- psT4_s6 %>%
  filter_taxa(function(x) sum(x > 0) > (percent*length(x)), TRUE)  #change the decimal to desired percentage
#s6Core
#tax_table(s6Core)

s6TPA <- psT4_s6 %>%
  subset_samples(Substrate == "TPA") %>%
  filter_taxa(function(x) sum(x > 0) > (percent*length(x)), TRUE)  #change the decimal to desired percentage
s6TPA

tax_table(s6TPA)

s6TA <- psT4_s6 %>%
  subset_samples(Substrate == "TA") %>%
  filter_taxa(function(x) sum(x > 0) > (percent*length(x)), TRUE)  #change the decimal to desired percentage
s6TA

tax_table(s6TA)

s6Cont <- psT4_s6 %>%
  subset_samples(Substrate == "Control") %>%
  filter_taxa(function(x) sum(x > 0) > (percent*length(x)), TRUE)  #change the decimal to desired percentage
#s6Cont
#tax_table(s6Cont)

#Making Summary of Core Taxa
compCore_list <- as.data.frame(tax_table(compCore)) %>% rownames_to_column(var = "OTU_ID") %>%
  mutate(description = "compost_core")
compTPA_list <- as.data.frame(tax_table(compTPA)) %>% rownames_to_column(var = "OTU_ID") %>%
  mutate(description = "compost_TPA")
compTA_list <- as.data.frame(tax_table(compTA)) %>% rownames_to_column(var = "OTU_ID") %>%
  mutate(description = "compost_TA")
compCont_list <- as.data.frame(tax_table(compCont)) %>% rownames_to_column(var = "OTU_ID") %>%
  mutate(description = "compost_control")

s1TA_list <- as.data.frame(tax_table(s1TA)) %>% rownames_to_column(var = "OTU_ID") %>%
  mutate(description = "sediment1_TA")
s1Cont_list <- as.data.frame(tax_table(s1Cont)) %>% rownames_to_column(var = "OTU_ID") %>%
  mutate(description = "sediment1_control")

s3Core_list <- as.data.frame(tax_table(s3Core)) %>% rownames_to_column(var = "OTU_ID") %>%
  mutate(description = "sediment3_core")
s3TPA_list <- as.data.frame(tax_table(s3TPA)) %>% rownames_to_column(var = "OTU_ID") %>%
  mutate(description = "sediment3_TPA")
s3TA_list <- as.data.frame(tax_table(s3TA)) %>% rownames_to_column(var = "OTU_ID") %>%
  mutate(description = "sediment3_TA")
s3Cont_list <- as.data.frame(tax_table(s3Cont)) %>% rownames_to_column(var = "OTU_ID") %>%
  mutate(description = "sediment3_control")

s4TPA_list <- as.data.frame(tax_table(s4TPA)) %>% rownames_to_column(var = "OTU_ID")  %>%
  mutate(description = "sediment4_TPA")
s4TA_list <- as.data.frame(tax_table(s4TA)) %>% rownames_to_column(var = "OTU_ID") %>%
  mutate(description = "sediment4_TA")
s4Cont_list <- as.data.frame(tax_table(s4Cont)) %>% rownames_to_column(var = "OTU_ID") %>%
  mutate(description = "sediment4_control")

s6TPA_list <- as.data.frame(tax_table(s6TPA)) %>% rownames_to_column(var = "OTU_ID") %>%
  mutate(description = "sediment6_TPA")
s6TA_list <- as.data.frame(tax_table(s6TA)) %>% rownames_to_column(var = "OTU_ID") %>%
  mutate(description = "sediment6_TA")

tax_list <- compCore_list %>%
  full_join(compTPA_list) %>%
  full_join(compTA_list) %>%
  full_join(compCont_list) %>%
  full_join(s1TA_list) %>%
  full_join(s1Cont_list) %>%
  full_join(s3Cont_list) %>%
  full_join(s3TPA_list) %>%
  full_join(s3TA_list) %>%
  full_join(s3Cont_list) %>%
  full_join(s4TPA_list)%>%
  full_join(s4TA_list) %>%
  full_join(s4Cont_list) %>%
  full_join(s6TPA_list) %>%
  full_join(s6TA_list) %>%
  group_by(Genus, description) %>%
  summarise(Count = n()) %>%
  separate(description, into = c("Environment", "Substrate"), sep = "_", remove = FALSE) %>%
  filter(Substrate != "core") %>%
  mutate(Substrate = factor(Substrate, levels = c("TA", "control", "TPA")))
head(tax_list)

library(ggh4x)

ggplot(tax_list, aes(x = description, y = Genus))+
  facet_nested(cols = vars(Environment, Substrate), space = "free_x", scales = "free_x")+
  geom_tile(aes(fill = Substrate), color = "black")+
  geom_text(aes(label = Count)) +
  xlab(NULL)+
  theme_bw()+
  theme(#axis.text.y = element_text(size = 18, color = "black"),
        #axis.text.x = element_text(size = 16, angle = 0, hjust = 0.5, vjust = 0.5, color = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 16),
        legend.position = "bottom",
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.background = element_rect(fill = "white"),
        #strip.text.y = element_text(size = 17, face = "bold", angle = 0, hjust = 0.5, vjust = 0.5),
        strip.text.x = element_text(size = 18, face = "bold", angle = 0),
        strip.text.y = element_blank(),
        #strip.text.x = element_blank(),
        title = element_text(size = 18))






merged_ps <- merge_phyloseq(compCore, compTPA, compTA, compCont,
                            s1TA, s1Cont, s3Core, s3TPA, s3TA, s3Cont,
                            s4TPA, s4TA, s4Cont, s6TPA, s6TA)

View(tax_table(merged_ps))

#tax_list <- as.data.frame(tax_table(coreASVs))$Genus
#length(tax_list)

genTax <- merged_ps %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Genus) %>%
  filter(Genus %in% tax_list$Genus) %>%
  filter(Abundance > 0) %>%
  select(Genus, Abundance, Substrate, Environment) %>%
  group_by(Genus, Substrate, Environment) %>%
  #summarize(avg_abundance = sum(Abundance)/n()) %>%
  mutate(Genus = as.character(Genus))
head(genTax)

colors <- c("lightblue", "orange", "olivedrab4", "firebrick", "cyan", "darkblue",
            "darkcyan", "lightgoldenrod", "tan4", "maroon")

ggplot(genTax)+
  geom_col(mapping = aes(x=Substrate, y = Abundance, fill = Genus), position = "stack", show.legend = TRUE, color = "black")+
  ylab("Relative Abundance") +
  #ylim(0,1)+
  #scale_fill_manual(values = colors) +
  xlab(NULL)+
  ggtitle("Genus")+
  theme_minimal()+
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 90, vjust = 0.5, hjust = 1),
        legend.text = element_text(size = 18),
        legend.position = "bottom",
        legend.title = element_blank(),
        title = element_text(size = 25))+
  guides(fill = guide_legend(nrow = 3))

