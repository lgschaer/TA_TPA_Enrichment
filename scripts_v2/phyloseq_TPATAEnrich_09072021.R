###Section 1

library(phyloseq)
library(tidyverse)
library(csv)

#load metadata
sample_names <- as.data.frame(sample.names)
head(sample_names)

write_csv(sample_names, "/home/lgschaer/old/Plastic_Deg/TA_TPA_Enrichment_09062021/sample_names.csv")

sdata <- as.csv("/home/lgschaer/old/Plastic_Deg/TA_TPA_Enrichment_09062021/TA_TPA_Metadata_09072021.csv", row.names = 1, header = TRUE, sep = ",", check.names = TRUE, stringsAsFactors = TRUE)
head(sdata)

sdata3 <- sdata %>%
  column_to_rownames(var = "SampleID")
View(sdata3)

#load sequence table
sequence_table <- readRDS("/home/lgschaer/old/Plastic_Deg/TA_TPA_Enrichment_09062021/seqtab.rds")    
sequence_table[1:3,1:3]

#1 remove column names
colnames(sequence_table) <- NULL
sequence_table[1:5,1:5]

#2 convert to data frame
sequence_table <- as.data.frame(sequence_table)
class(sequence_table)

#3 add "SampleID" header to rownames
#sequence_table2 <- rownames_to_column(sequence_table, var = "SampleID")
#sequence_table2[1:5,1:5]

#4 change sample names to match metadata (you won't normally have to do this)
#sequence_table <- sequence_table %>%
#  mutate(SampleID = str_replace_all(SampleID, "-", "_"))
#sequence_table[1:4,1:4]

#5 add metadata to sequence table, making sure dimensions and sample names line up
#seq_table <- sdata2 %>% 
#  left_join(sequence_table2, by = "SampleID")
#seq_table[1:5,1:5]

#dim(seq_table)
dim(sdata)
dim(sequence_table)

#7 remove metadata columns
subset_all <- sequence_table

#8 remove column names from the sequence table
colnames(subset_all) <- NULL
subset_all[1:5,1:5]

#load taxa table
taxa <- readRDS("/home/lgschaer/old/Plastic_Deg/TA_TPA_Enrichment_09062021/taxa.rds")
taxa[1:5,1:5]

### Section 2: making a phyloseq object

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

write_rds(ps, "/home/lgschaer/old/Plastic_Deg/TA_TPA_Enrichment_09062021/ps_unrarefied_09072021.rds")


### Section 3: Normalize the data

samplesover1000_all <- subset_samples(ps, sample_sums(ps) > 1000)

any(taxa_sums(samplesover1000_all) == 0)

sum(taxa_sums(samplesover1000_all) == 0)

prune_samplesover1000_all <- prune_taxa(taxa_sums(samplesover1000_all) > 0, samplesover1000_all)

any(taxa_sums(prune_samplesover1000_all) == 0)

#for reproducible data
set.seed(81)

rarefy_samplesover1000_all <- rarefy_even_depth(prune_samplesover1000_all, rngseed= 81, sample.size = min(sample_sums(prune_samplesover1000_all)))

pps<- rarefy_samplesover1000_all
pps

#
count_samples1 <- sample_data(ps) %>%
  group_by(Substrate, Transfer) %>%
  summarise(
    count = n()
  )
View(count_samples1)


count_samples2 <- sample_data(pps) %>%
  group_by(Substrate, Transfer) %>%
  summarise(
    count = n()
  )
View(count_samples2)

###Section 4: plotting

#sample colors
sample_colors <- c("Inocula" = "blue", "Control" = "indianred2", "TPA" = "orange", "TA" = "lightblue")

sample_types <- c("Inocula", "Control", "TPA", "TA")

#sample_labels <- c("bilge" ="Bilge", "water" = "Water", "boatback" = "Boat\nBack", "boatside" = "Boat\nSide", "boatrear" = "Boat\nRear", "dock" = "Dock", "air" = "Air")

#violin plot
pps %>%                                                              #phyloseq object
  plot_richness(
    x = "Substrate",                                                #compare diversity of datatype
    measures = c("Observed", "Shannon")) +                           #choose diversity measures
  geom_violin(aes(fill = Substrate), show.legend = FALSE)+          #make violin plot, set fill aes to sampletype
  geom_boxplot(width=0.1) +                                          #add boxplot, set width
  theme_classic()+                                                   #change theme to classic
  xlab(NULL)+                                                        #no label on x-axis
  theme(axis.text.y.left = element_text(size = 20),                  #adjust y-axis text
        axis.text.x = element_text(size = 20, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 20))+                     #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 25))+        #adjust headings
  scale_fill_manual(values = sample_colors)+   #set fill colors
 # scale_x_discrete(                                                  #change x-axis labels
   # breaks = sample_types)+                   
  ggtitle("Alpha Diversity") +                                       #add title
  theme(plot.title=element_text(size = 25, face = "bold", hjust = 0.5)) #change title size, face and position

#violin plot
pps %>%                                                              #phyloseq object
  plot_richness(
    x = "Transfer",                                                #compare diversity of datatype
    measures = c("Shannon")) +                           #choose diversity measures
  geom_violin(aes(fill = Substrate), show.legend = TRUE)+          #make violin plot, set fill aes to sampletype
  geom_boxplot(width=0.1) +                                          #add boxplot, set width
  theme_classic()+                                                   #change theme to classic
  xlab(NULL)+                                                        #no label on x-axis
  theme(axis.text.y.left = element_text(size = 20),                  #adjust y-axis text
        axis.text.x = element_text(size = 20, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 20))+                     #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 25))+        #adjust headings
  scale_fill_manual(values = colors)+   #set fill colors
  # scale_x_discrete(                                                  #change x-axis labels
  # breaks = sample_types)+                   
 # ggtitle("Alpha Diversity") +                                       #add title
  theme(plot.title=element_text(size = 25, face = "bold", hjust = 0.5)) #change title size, face and position

colors <- c("red", "blue", "orange", "green", "lightblue")

#violin plot
pps %>%                                                              #phyloseq object
  plot_richness(
    x = "Time_Point",                                                #compare diversity of datatype
    measures = c("Observed", "Shannon")) +                           #choose diversity measures
  geom_violin(aes(fill = Time_Point), show.legend = FALSE)+          #make violin plot, set fill aes to sampletype
  geom_boxplot(width=0.1) +                                          #add boxplot, set width
  theme_classic()+                                                   #change theme to classic
  xlab(NULL)+                                                        #no label on x-axis
  theme(axis.text.y.left = element_text(size = 20),                  #adjust y-axis text
        axis.text.x = element_text(size = 20, hjust = 1, vjust = 0.5, angle = 90),           #adjust x-axis label position
        axis.title.y = element_text(size = 20))+                     #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 25))+        #adjust headings
  scale_fill_manual(values = colors)+   #set fill colors
  # scale_x_discrete(                                                  #change x-axis labels
  # breaks = sample_types)+                   
  ggtitle("Alpha Diversity") +                                       #add title
  theme(plot.title=element_text(size = 25, face = "bold", hjust = 0.5)) #change title size, face and position


#Diversity bar plot with error bars
pps_filtered <- subset_samples(pps, Substrate != "Inocula" & Environment != "Inocula")

pps_diversity <- estimate_richness(pps_filtered, measures = "Observed") %>%
  rownames_to_column(var = "SampleID") %>%
  left_join(sdata, by = "SampleID") %>%
  group_by(Substrate, Environment, Transfer) %>%
  summarise(
    avgObs = mean(Observed),
    maxObs = max(Observed),
    minObs = min(Observed)
  )
head(pps_diversity)
#dim(pps_diversity)

ggplot(pps_diversity, aes(x = Transfer, y = avgObs))+
  geom_col(aes(fill = Environment), show.legend = TRUE, position = "dodge", color = "black")+          #make violin plot, set fill aes to sampletype
  geom_errorbar(aes(ymin = minObs, ymax = maxObs, color = Environment), position = "dodge") +                                          #add boxplot, set width
  facet_grid(rows = vars(Substrate))+
  theme_bw()+                                                   #change theme to classic
  xlab(NULL)+                                                        #no label on x-axis
  theme(axis.text.y.left = element_text(size = 20),                  #adjust y-axis text
        axis.text.x = element_text(size = 20, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 20))+                     #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 15))+        #adjust headings
  scale_fill_manual(values = colors)+   #set fill colors
  scale_color_manual(values = colors)+
  # scale_x_discrete(                                                  #change x-axis labels
  # breaks = sample_types)+                   
  ggtitle("Alpha Diversity") +                                       #add title
  theme(plot.title=element_text(size = 25, face = "bold", hjust = 0.5)) #change title size, face and position

pps <- subset_samples(pps, Substrate != "Blank")

#Making a PCoA plot
#ordination
all_pcoa <- ordinate(
  physeq = pps, 
  method = "PCoA", 
  distance = "bray"
)

head(sdata)

colors <- c("purple", "orange", "blue", "red", "white")

#plot
plot_ordination(
  physeq = pps,                                                         #phyloseq object
  ordination = all_pcoa)+                                                #ordination
  geom_point(aes(fill = Substrate, shape = Transfer), size = 5) +     #sets fill color to sampletype
  scale_shape_manual(values = c(21, 22, 23, 24, 25))+
  scale_fill_manual(values = colors) +
  theme_classic() +                                                      #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 20),                               #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command

#plot
plot_ordination(
  physeq = pps,                                                         #phyloseq object
  ordination = all_pcoa)+                                                #ordination
  geom_point(aes(fill = Transfer, shape = Substrate), size = 3) +     #sets fill color to sampletype
  scale_shape_manual(values = c(21, 22, 23, 24))+
  scale_fill_manual(values = colors) +
  theme_classic() +                                                      #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 20),                               #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command

# Dissimilarity Plot

vectors <- as.data.frame(all_pcoa$vectors) %>%
  rownames_to_column(var = "SampleID") %>%
  select(SampleID, Axis.1) %>%
  left_join(sdata, by = "SampleID") %>%
  filter(Substrate != "Inocula" & Environment != "Inocula")
head(vectors)

ggplot(vectors, aes(x = Transfer, y = Axis.1))+
  #facet_grid(rows = vars(Substrate))+
  geom_boxplot(aes(fill = Substrate), show.legend = TRUE)+          #make violin plot, set fill aes to sampletype
  #geom_boxplot(width=0.05) +                                          #add boxplot, set width
  xlab(NULL)+                                                        #no label on x-axis
  theme(axis.text.y.left = element_text(size = 20),                  #adjust y-axis text
        axis.text.x = element_text(size = 20, hjust = 1, vjust = 0.5, angle = 90),           #adjust x-axis label position
        axis.title.y = element_text(size = 20))+                     #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 25))+        #adjust headings
  scale_fill_manual(values = sample_colors)+ 
  theme_bw()

ggplot(vectors, aes(x = Transfer, y = Axis.1, fill = Substrate), color = "black")+
  #facet_grid(rows = vars(Substrate))+
  geom_col(position = "dodge") +                                          #add boxplot, set width
  xlab(NULL)+                                                        #no label on x-axis
  theme(axis.text.y.left = element_text(size = 20),                  #adjust y-axis text
        axis.text.x = element_text(size = 20, hjust = 1, vjust = 0.5, angle = 90),           #adjust x-axis label position
        axis.title.y = element_text(size = 20))+                     #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 25))+        #adjust headings
  scale_fill_manual(values = sample_colors)+ 
  theme_bw()

#plot
plot_ordination(
  physeq = pps,                                                         #phyloseq object
  ordination = all_pcoa)+                                                #ordination
  geom_point(aes(fill = Environment, shape = Substrate), size = 3) +     #sets fill color to sampletype
  scale_shape_manual(values = c(21, 22, 23, 24))+
  scale_fill_manual(values = colors) +
  theme_classic() +                                                      #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 20),                               #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command

head(sdata)

## Time PCOA

#subset phyloseq
T1ps <- subset_samples(pps, Transfer == "T1")
T2ps <- subset_samples(pps, Transfer == "T2")
T3ps <- subset_samples(pps, Transfer == "T3")
T4ps <- subset_samples(pps, Transfer == "T4")

#ordination
T1_pcoa <- ordinate(
  physeq = T1ps, 
  method = "PCoA", 
  distance = "bray"
)

T2_pcoa <- ordinate(
  physeq = T2ps, 
  method = "PCoA", 
  distance = "bray"
)

T3_pcoa <- ordinate(
  physeq = T3ps, 
  method = "PCoA", 
  distance = "bray"
)

T4_pcoa <- ordinate(
  physeq = T4ps, 
  method = "PCoA", 
  distance = "bray"
)


#plot
T1 <- plot_ordination(
  physeq = T1ps,                                                         #phyloseq object
  ordination = T1_pcoa)+                                                #ordination
  geom_point(aes(fill = Environment, shape = Substrate), size = 3) +     #sets fill color to sampletype
  scale_shape_manual(values = c(21, 22, 23, 24))+
  scale_fill_manual(values = colors) +
  theme_classic() +                                                      #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 20),                               #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command

T2 <- plot_ordination(
  physeq = T2ps,                                                         #phyloseq object
  ordination = T2_pcoa)+                                                #ordination
  geom_point(aes(fill = Environment, shape = Substrate), size = 3) +     #sets fill color to sampletype
  scale_shape_manual(values = c(21, 22, 23, 24))+
  scale_fill_manual(values = colors) +
  theme_classic() +                                                      #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 20),                               #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command

T3 <- plot_ordination(
  physeq = T3ps,                                                         #phyloseq object
  ordination = T3_pcoa)+                                                #ordination
  geom_point(aes(fill = Environment, shape = Substrate), size = 3) +     #sets fill color to sampletype
  scale_shape_manual(values = c(21, 22, 23, 24))+
  scale_fill_manual(values = colors) +
  theme_classic() +                                                      #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 20),                               #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command

T4 <- plot_ordination(
  physeq = T4ps,                                                         #phyloseq object
  ordination = T4_pcoa)+                                                #ordination
  geom_point(aes(fill = Environment, shape = Substrate), size = 3) +     #sets fill color to sampletype
  scale_shape_manual(values = c(21, 22, 23, 24))+
  scale_fill_manual(values = colors) +
  theme_classic() +                                                      #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 20),                               #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command

T1
T2
T3
T4

library(ggpubr)
ggarrange(T1, T2, T3, T4, common.legend = TRUE, labels = c("T1", "T2", "T3", "T4"))



## Substrate PCOAs

#subset phyloseq
TAps <- subset_samples(pps, Substrate == "TA")
TPAps <- subset_samples(pps, Substrate == "TPA")
Cps <- subset_samples(pps, Substrate == "Control")

#ordination
TA_pcoa <- ordinate(
  physeq = TAps, 
  method = "PCoA", 
  distance = "bray"
)

TPA_pcoa <- ordinate(
  physeq = TPAps, 
  method = "PCoA", 
  distance = "bray"
)

C_pcoa <- ordinate(
  physeq = Cps, 
  method = "PCoA", 
  distance = "bray"
)

#plot
TA <- plot_ordination(
  physeq = TAps,                                                         #phyloseq object
  ordination = TA_pcoa)+                                                #ordination
  geom_point(aes(fill = Environment, shape = Transfer), size = 3) +     #sets fill color to sampletype
  scale_shape_manual(values = c(21, 22, 23, 24))+
  scale_fill_manual(values = colors) +
  theme_classic() +                                                      #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 20),                               #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command

TPA <- plot_ordination(
  physeq = TPAps,                                                         #phyloseq object
  ordination = TPA_pcoa)+                                                #ordination
  geom_point(aes(fill = Environment, shape = Transfer), size = 3) +     #sets fill color to sampletype
  scale_shape_manual(values = c(21, 22, 23, 24))+
  scale_fill_manual(values = colors) +
  theme_classic() +                                                      #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 20),                               #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command

C <- plot_ordination(
  physeq = Cps,                                                         #phyloseq object
  ordination = C_pcoa)+                                                #ordination
  geom_point(aes(fill = Environment, shape = Transfer), size = 3) +     #sets fill color to sampletype
  scale_shape_manual(values = c(21, 22, 23, 24))+
  scale_fill_manual(values = colors) +
  theme_classic() +                                                      #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 20),                               #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command

ggarrange(TA, TPA, C, 
          common.legend = TRUE, 
          labels = c("TA", "TPA", "Control"), 
          nrow = 2,ncol=2,
          font.label = list(size = 15, face = "bold"))



#taxa plot

#1) Filter out eukaryotes and mitochondria
justbacteria <- ps %>%
  subset_taxa(
    Kingdom == "Bacteria" &                   #only bacteria
      Family  != "mitochondria" &             #filter out mitochondria
      Class   != "Chloroplast"                #filter out chloroplasts
  )
justbacteria

#2) Convert phyloseq object into a data frame with relative abundance counts
phylumabundance <- justbacteria %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Phylum) 
head(phylumabundance)

#Summarize the most abundant taxa, filter to only include taxa with > 1% abundance
all <- phylumabundance %>%
  select(Phylum, Environment, Abundance, Substrate, Transfer) %>%
  group_by(Phylum, Substrate, Environment, Transfer) %>%
  summarize(
    avg_abundance = mean(Abundance)
  ) #%>%
  #filter(avg_abundance > 0.01) #%>%
 # filter(Substrate != "Blank" & Environment != "Blank")
head(all)

#4) Save phylum colors
phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", "blue",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
)

#5) Create taxa plot!
ggplot(all)+
  geom_col(mapping = aes(x = Substrate, y = avg_abundance, fill = Phylum), position = "fill", show.legend = TRUE, color = "black")+
  facet_grid(rows = vars(Environment), cols = vars(Transfer))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = phylum_colors) +
  xlab(NULL)+
  theme_minimal()+
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 20),
        title = element_text(size = 25))

abundanceplot <- phylumabundance %>%
  select(Phylum, Abundance, Environment, Transfer) %>%
  group_by(Phylum, Environment, Transfer) %>%
  summarize(
    avg_abundance = mean(Abundance)
  ) %>%
  filter(avg_abundance > 0.10) %>% #filter everything less than 10% relative abundance
  filter(Environment != "Blank") 
head(abundanceplot)


ggplot(abundanceplot, mapping = aes(x = Transfer, y = avg_abundance, fill = Phylum), color = "black", show.legend = TRUE)+
  geom_point(position = "fill", shape = 21, size = 3, show.legend = TRUE)+
  facet_grid(rows = vars(Environment))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = phylum_colors) +
  xlab(NULL)+
  theme_minimal()+
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 1, hjust = 0.5),
        axis.title.y = element_text(size = 20),
        title = element_text(size = 25))


#2) Convert phyloseq object into a data frame with relative abundance counts
classabundance <- justbacteria %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Class) 
head(classabundance)

#Summarize the most abundant taxa, filter to only include taxa with > 1% abundance
all2 <- classabundance %>%
  select(Class, Substrate, Abundance) %>%
  mutate( Class = as.character(Class),
    Class = ifelse(Abundance >= 0.05, Class, "< 0.05%"))%>%
  group_by(Class, Substrate) %>%
  summarize(
    avg_abundance = mean(Abundance)
  ) %>%
  filter(Substrate != "Blank")
head(all2)

TA <- filter(all2, Substrate == "Terephthalamide")
TPA <- filter(all2, Substrate == "Terephthalate")
C <- filter(all2, Substrate == "Control")

#4) Save phylum colors
phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
   "#673770","#D14285",  "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "Black", "lightblue",
  "#CBD588", "#5F7FC7"
)

#5) Create taxa plot!
ggplot(all2)+
  geom_col(mapping = aes(x = Substrate, y = avg_abundance, fill = Class), position = "fill", show.legend = TRUE, color = "black")+
  #facet_grid(cols = vars(Inocula, Substrate), rows = vars(Time))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = phylum_colors) +
  xlab(NULL)+
  theme_minimal()+
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5),
        axis.title.y = element_text(size = 20),
        title = element_text(size = 25))

#2) Convert phyloseq object into a data frame with relative abundance counts
genusabundance <- justbacteria %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at genus level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Genus) 
head(genusabundance)

#Summarize the most abundant taxa, filter to only include taxa with > 1% abundance
all3 <- genusabundance %>%
  filter(Substrate != "Inocula") %>%
  select(Genus, Substrate, Abundance, Time) %>%
  mutate(Order = as.character(Order),
         Order = ifelse(Abundance >= 0.02, Order, "< 2%"))%>%
  group_by(Order, Substrate,Time) %>%
  summarize(
    avg_abundance = sum(Abundance)) #%>%
  #mutate()
head(all3)

#4) Save colors
genus_colors <- c(
  "#CBD588","orange", "grey","#DA5724", "blue", "#508578", 
   "#D14285",  "#C84248", 
  "#8569D5", "#D1A33D", "#8A7C64", "#599861",
  "purple4",      "orchid1",   "green",          "coral1", "blue", "yellow",
    "cyan",          "palegoldenrod",
   "darkblue",                "mediumpurple1",
  "grey",   "dodgerblue",   "firebrick", "yellowgreen", "magenta", "lightpink",
   "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "purple4", "darkcyan",     "orchid1",   "green",          "coral1", "blue", "yellow", "red",
  "grey47",  "cyan",         "darkgreen", "palegoldenrod",  "tan4",
  "grey77",  "darkblue",     "orange",    "red",         "mediumpurple1",
  "white",   "dodgerblue",   "firebrick", "yellowgreen", "magenta", "black", "lightpink",
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "purple4", "darkcyan",     "orchid1",   "green",          "coral1", "blue", "yellow", "red",
  "grey47",  "cyan",         "darkgreen", "palegoldenrod",  "tan4",
  "grey77",  "darkblue",     "orange",    "red",         "mediumpurple1",
  "white",   "dodgerblue",   "firebrick", "yellowgreen", "magenta", "black", "lightpink"
)

#5) Create taxa plot!
ggplot(all3)+
  geom_col(mapping = aes(x = Time, y = avg_abundance, fill = Order), color = "black", position = "fill", show.legend = TRUE)+
  #facet_grid(rows = vars(Inocula), cols = vars(Substrate))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = genus_colors) +
  xlab(NULL)+
  theme_minimal()+
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5),
        axis.title.y = element_text(size = 20),
        title = element_text(size = 25))


#Summarize the most abundant taxa, filter to only include taxa with > 1% abundance
all4 <- genusabundance %>%
  select(Genus, Substrate, Abundance, Inocula, Time_Point, Replicate) %>%
  filter(Substrate != "Blank" & Inocula != "Blank") %>%
  group_by(Genus, Substrate, Inocula, Time_Point, Replicate) %>%
  summarize(
    avg_abundance = mean(Abundance)
  ) %>%
  mutate(Genus = as.character(Genus),
         Genus = ifelse(mean(avg_abundance) >= 0.025, Genus, "< 0.025%"),
         Replicate = ifelse(Replicate == "R2A", "R2", Replicate),
         Replicate = ifelse(Replicate == "R2B", "R3", Replicate)) %>%
  filter(Replicate == "R1" | Replicate == "R2" | Replicate == "R3") %>%
  group_by(Genus, Substrate, Inocula, Time_Point, Replicate) %>%
  summarise(
    avg_abundance = sum(avg_abundance)
  ) %>%
  filter(!is.na(Genus))
head(all4)
dim(all4)

TA <- filter(all4, Substrate == "Terephthalamide")
TPA <- filter(all4, Substrate == "Terephthalate")
C <- filter(all4, Substrate == "Control")


#5) Create taxa plot!
ggplot(all4)+
  geom_col(mapping = aes(x = Time_Point, y = avg_abundance, fill = Genus), color = "black", position = "fill", show.legend = FALSE)+
  facet_grid(rows = vars(Inocula), cols = vars(Substrate, Replicate))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = genus_colors) +
  xlab(NULL)+
  theme_bw()+
  theme(axis.text.y.left = element_text(size = 12), legend.position = "bottom",
        axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5),
        axis.title.y = element_text(size = 12),
        title = element_text(size = 25))



#create plot showing number of genera in each sample

all5 <- genusabundance %>%
  select(Genus, Substrate, Abundance, Inocula, Time_Point, Replicate) %>%
  filter(Substrate != "Blank" & Inocula != "Blank") %>%
  filter(Abundance > 0) %>%
  mutate(Genus = as.character(Genus),
         Time = as.factor(Time_Point),
         #Genus = ifelse(mean(avg_abundance) >= 0.1, Genus, "< 0.1%"),
         Replicate = ifelse(Replicate == "R2A", "R2", Replicate),
         Replicate = ifelse(Replicate == "R2B", "R3", Replicate)) %>%
  filter(Replicate == "R1" | Replicate == "R2" | Replicate == "R3") %>%
  group_by(Substrate, Inocula, Time_Point, Replicate) %>%
  summarise(
    count = n_distinct(Genus), .groups = "keep"
  ) 
head(all5)


#5) Create plot!
ggplot(all5, mapping = aes(x = Time, y = count, fill = Substrate))+
  geom_col(mapping = aes(x = Time, y = count, fill = Substrate), color = "black", position = "dodge", show.legend = TRUE)+
  facet_grid(rows = vars(Inocula), cols = vars(Substrate, Replicate))+
  ylab("Number of Genera") +
  scale_fill_manual(values = genus_colors) +
  xlab(NULL)+
  theme_bw()+
  theme(axis.text.y.left = element_text(size = 12), legend.position = "bottom",
        axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5),
        axis.title.y = element_text(size = 12),
        title = element_text(size = 25))



#Make similar plot showing abundance of rare taxa only
all6 <- genusabundance %>%
  select(Genus, Substrate, Abundance, Inocula, Time, Replicate) %>%
  filter(Substrate != "Blank" & Inocula != "Blank") %>%
  filter(Abundance > 0 & Abundance <= 0.1) %>%
  mutate(Genus = as.character(Genus),
         Time = as.factor(Time),
         #Genus = ifelse(mean(avg_abundance) >= 0.1, Genus, "< 0.1%"),
         Replicate = ifelse(Replicate == "R2A", "R2", Replicate),
         Replicate = ifelse(Replicate == "R2B", "R3", Replicate)) %>%
  filter(Replicate == "R1" | Replicate == "R2" | Replicate == "R3") %>%
  group_by(Substrate, Inocula, Time, Replicate) %>%
  summarise(
    count = n_distinct(Genus), .groups = "keep"
  ) 
head(all6)


#5) Create taxa plot!
ggplot(all6, mapping = aes(x = Time, y = count, fill = Substrate))+
  geom_col(mapping = aes(x = Time, y = count, fill = Substrate), color = "black", position = "dodge", show.legend = TRUE)+
  facet_grid(rows = vars(Inocula), cols = vars(Substrate, Replicate))+
  ylab("Number of Genera") +
  scale_fill_manual(values = genus_colors) +
  ggtitle("Number of Low Abundance Genera") +
  xlab(NULL)+
  theme_bw()+
  theme(axis.text.y.left = element_text(size = 12), legend.position = "bottom",
        axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5),
        axis.title.y = element_text(size = 12),
        title = element_text(size = 25))

