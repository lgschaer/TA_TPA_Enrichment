#Packages used
library(csv)
library(tidyverse)
library(phyloseq)

#Load phyloseq objects (previously saved as rds objects)
ps1<-read_rds("/home/lgschaer/old/Plastic_Deg/TA_TPA_Enrichment_02052021/ps_unrarefied_02052021.rds")
ps1         #phyloseq object from run 1

ps2<-read_rds("/home/lgschaer/old/Plastic_Deg/TA_TPA_Enrichment_03232021/ps_unrarefied_03232021.rds")
ps2         #phyloseq object from run 2

#Merge phyloseq objects
merged_ps <- merge_phyloseq(ps1,ps2)   #merge phyloseq objects from both runs together
merged_ps                              #378 samples, 4 sample variables

head(sample_data(merged_ps)) #check that columns matched and merged correctly

#now add an "ID" column to the metadata that we can separate into metadata later
new_meta <- sample_data(merged_ps) %>%                                   #extract metadata from merged phyloseq object
  as_tibble(rownames=c("SampleID")) %>%                                  #convert to tibble format
  select(c("SampleID", "Inocula", "Substrate", "Replicate", "Time")) %>% #select the rows that we want to keep
  unite("ID", Inocula:Time, sep = "_", remove = FALSE) %>%               #unite relevant rows into an "ID" column
  column_to_rownames("SampleID")                                         #convert "SampleID" back to rownames
dim(new_meta)       #dimensions should be the same as the phyloseq object
head(new_meta)      #look at the data to make sure everything is as expected

merged_ps2 <- merged_ps                      #rename the phyloseq object before adding the new metadata
merged_ps2@sam_data <- sample_data(new_meta) #add the revised metadata to the phyloseq object
head(sample_data(merged_ps2))                #check to make sure the new metadata was added correctly
merged_ps2                                   #378 samples, 5 sample variables


combined_samples <- merge_samples(merged_ps2, "ID")   #merge samples together based on the "ID" column
combined_samples                                      #181 samples and 5 sample variables
head(sample_data(combined_samples))                   #metadata is all NA

final_meta <- sample_data(combined_samples) %>%                           #extract data from combined sample phyloseq object
  as_tibble(rownames=c("SampleID")) %>%                                   #convert to tibble format
  select(-c("ID", "Inocula", "Substrate", "Replicate", "Time")) %>%       #remove the columns with NAs
  separate(SampleID, into = c("Inocula", "Substrate", "Replicate", "Time"), remove = FALSE) %>% #separate SampleID column into new metadata columns
  column_to_rownames(var = "SampleID")                                    #convert SampleID column back to rownames
head(final_meta)          #check head and dimensions of new metadata
dim(final_meta)           #181 samples, 4 sample variables, no NAs

combined_samples2 <- combined_samples                 #rename phyloseq object
combined_samples2@sam_data <- sample_data(final_meta) #add the revised metadata to the phyloseq object
head(sample_data(combined_samples2))                  #make sure it worked
combined_samples2                                     #181 samples, 4 sample variables


#Stats on merged phyloseq object
mean(sample_sums(combined_samples2))         #mean sample sum: 32,742
median(sample_sums(combined_samples2))       #median sample sum: 13,166
max(sample_sums(combined_samples2))          #maximum sample sum: 283,465
min(sample_sums(combined_samples2))          #minimum sample sum: 0
sum(sample_sums(combined_samples2) == 0)     #36 samples have 0 reads
sum(sample_sums(combined_samples2) < 1000)   #72 samples have less than  1000 reads
sum(sample_sums(combined_samples2) < 500)    #64 samples have less than  500 reads

sampleSums <- as.data.frame(sample_sums(combined_samples2)) %>%
  rownames_to_column(var = "SampleID") %>%
  mutate(sample_sums = sample_sums(combined_samples2)) %>%
  select(SampleID, sample_sums)
head(sampleSums)
dim(sampleSums)

lowSums <- sampleSums %>%
  filter(sample_sums < 100 & sample_sums > 20)
View(lowSums)

as.csv(sampleSums, "/home/lgschaer/old/Plastic_Deg/sampleSums.csv")

plotting_info <- sample_data(combined_samples2) %>%
  as_tibble(rownames=c("SampleID")) %>%
  left_join(sampleSums, by = "SampleID") %>%
  filter(sample_sums < 50000) 
head(plotting_info)
dim(plotting_info)

ggplot(plotting_info, aes(x = sample_sums))+
  geom_histogram(binwidth = 1000)+
  facet_grid(rows = vars(Substrate))+
  scale_x_continuous(labels = scales::comma)+
  theme_classic()+
  theme(axis.text.y.left = element_text(size = 20),                  #adjust y-axis text
        axis.text.x = element_text(size = 12, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 20))                     #adjust y-axis title

plotting_info2 <- sample_data(combined_samples2) %>%
  as_tibble(rownames=c("SampleID")) %>%
  left_join(sampleSums, by = "SampleID") %>%
  filter(sample_sums < 50000) %>%
  mutate(
    Replicate = ifelse(Replicate == "R2A" | Replicate == "R2B", "R2", Replicate)
  ) %>%
  filter(Replicate=="R1"|Replicate=="R2"|Replicate=="R3")
head(plotting_info2)

ggplot(plotting_info2, aes(x = sample_sums, color = Substrate))+
  geom_histogram(binwidth = 1000)+
  facet_grid(rows = vars(Inocula), cols = vars(Time, Replicate))+
  scale_color_manual(values = sample_colors)+
  scale_x_continuous(labels = scales::comma)+
  theme_bw()+
  theme(axis.text.y.left = element_text(size = 10),                  #adjust y-axis text
        axis.text.x = element_text(size = 12, hjust = 1, vjust = 0.5, angle = 90),           #adjust x-axis label position
        axis.title.y = element_text(size = 10))                     #adjust y-axis title



### Section 3: Normalize the data
samplesover1000_all <- subset_samples(combined_samples2, sample_sums(combined_samples2) > 1000)

any(taxa_sums(samplesover1000_all) == 0)

sum(taxa_sums(samplesover1000_all) == 0)

prune_samplesover1000_all <- prune_taxa(taxa_sums(samplesover1000_all) > 0, samplesover1000_all)

any(taxa_sums(prune_samplesover1000_all) == 0)

#for reproducible data
set.seed(81)

rarefy_samplesover1000_all <- rarefy_even_depth(prune_samplesover1000_all, rngseed= 81, sample.size = min(sample_sums(prune_samplesover1000_all)))

#rename phyloseq object
pps <- rarefy_samplesover1000_all
pps

write_rds(pps, "/home/lgschaer/old/Plastic_Deg/TA_TPA_Enrichment_03232021/ps_rarefied_03252021.rds")

###Section 4: plotting

#sample colors
#sample_colors <- c("Control" = "indianred2", "Terephthalate" = "orange", "Terephthalamide" = "lightblue")
#sample_types <- c("Control", "Terephthalate", "Terephthalamide")
#sample_labels <- c("bilge" ="Bilge", "water" = "Water", "boatback" = "Boat\nBack", "boatside" = "Boat\nSide", "boatrear" = "Boat\nRear", "dock" = "Dock", "air" = "Air")

sample_colors <- c("red", "purple", "orange", "blue", "lightblue", "green", "cyan", "magenta")

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
        axis.text.x = element_text(size = 12, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
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
    x = "Inocula",                                                #compare diversity of datatype
    measures = c("Observed", "Shannon")) +                           #choose diversity measures
  geom_violin(aes(fill = Inocula), show.legend = FALSE)+          #make violin plot, set fill aes to sampletype
  geom_boxplot(width=0.1) +                                          #add boxplot, set width
  theme_classic()+                                                   #change theme to classic
  xlab(NULL)+                                                        #no label on x-axis
  theme(axis.text.y.left = element_text(size = 20),                  #adjust y-axis text
        axis.text.x = element_text(size = 12, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 20))+                     #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 25))+        #adjust headings
  scale_fill_manual(values = sample_colors)+   #set fill colors
  # scale_x_discrete(                                                  #change x-axis labels
  # breaks = sample_types)+                   
  ggtitle("Alpha Diversity") +                                       #add title
  theme(plot.title=element_text(size = 25, face = "bold", hjust = 0.5)) #change title size, face and position


#violin plot
phyloseq_object_all %>%                                                              #phyloseq object
  plot_richness(
    x = "Time",                                                #compare diversity of datatype
    measures = c("Observed", "Shannon")) +                           #choose diversity measures
  geom_violin(aes(fill = Time), show.legend = FALSE)+          #make violin plot, set fill aes to sampletype
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
pps_filtered <- subset_samples(pps, Substrate != "Blank" & Inocula != "Blank")

pps_diversity <- estimate_richness(pps_filtered, measures = "Observed") %>%
  rownames_to_column(var = "SampleID") %>%
  left_join(sdata, by = "SampleID") %>%
  group_by(Substrate, Inocula, Time) %>%
  summarise(
    avgObs = mean(Observed),
    maxObs = max(Observed),
    minObs = min(Observed)
  ) %>%
  filter(!is.na(Substrate))
head(pps_diversity)
#dim(pps_diversity)

ggplot(pps_diversity, aes(x = Time, y = avgObs))+
  geom_col(aes(fill = Inocula), show.legend = TRUE, position = "dodge", color = "black")+          #make violin plot, set fill aes to sampletype
  geom_errorbar(aes(ymin = minObs, ymax = maxObs, color = Inocula), position = "dodge") +                                          #add boxplot, set width
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

#t-SNE plot
library(tsnemicrobiota)

#filter phyloseq object
justbacteria <- combined_samples2 %>%
  subset_taxa(
    Kingdom == "Bacteria" &                   #only bacteria
      Family  != "mitochondria" &             #filter out mitochondria
      Class   != "Chloroplast"                #filter out chloroplasts
  ) %>%
  subset_samples(Substrate != "Blank" & !is.na(Time) & !is.na(Substrate)) %>%
  subset_samples(Inocula == "Compost")
justbacteria

tsne <- tsne_phyloseq(justbacteria, distance = "bray", perplexity = 25, dimensions = 2,
                      precomputed_distance = NULL, pseudocounts = 1, verbose = 1,
                      rng_seed = 81, philr_options = list(), control = list())
summary(tsne)

sample_colors

#tSNE Plot

plot_tsne_phyloseq(justbacteria, tsne, color = "Substrate", shape = "Inocula") +
  geom_point(aes(fill = Substrate), color = "black", size = 5, show.legend = TRUE) +     #sets fill color to sampletype
  scale_shape_manual(values = c(21, 22, 23, 24, 25))+
  scale_fill_manual(values = sample_colors) +
  theme_classic() +                                                      #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 20),                               #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(color = FALSE)+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command

#t-SNE plot by time

#subset phyloseq
T1ps <- subset_samples(pps, Time == "T1")
T2ps <- subset_samples(pps, Time == "T2")
T3ps <- subset_samples(pps, Time == "T3")
T4ps <- subset_samples(pps, Time == "T4")

tsneT1 <- tsne_phyloseq(T1ps, distance = "bray", perplexity = 25, dimensions = 2,
                      precomputed_distance = NULL, pseudocounts = 1, verbose = 1,
                      rng_seed = 81, philr_options = list(), control = list())

tsneT2 <- tsne_phyloseq(T2ps, distance = "bray", perplexity = 25, dimensions = 2,
                      precomputed_distance = NULL, pseudocounts = 1, verbose = 1,
                      rng_seed = 81, philr_options = list(), control = list())

tsneT3 <- tsne_phyloseq(T3ps, distance = "bray", perplexity = 25, dimensions = 2,
                      precomputed_distance = NULL, pseudocounts = 1, verbose = 1,
                      rng_seed = 81, philr_options = list(), control = list())

tsneT4 <- tsne_phyloseq(T4ps, distance = "bray", perplexity = 25, dimensions = 2,
                      precomputed_distance = NULL, pseudocounts = 1, verbose = 1,
                      rng_seed = 81, philr_options = list(), control = list())

#tSNE Plot

T1 <- plot_tsne_phyloseq(T1ps, tsneT1, color = "Substrate", shape = "Inocula") +
  geom_point(aes(fill = Substrate), color = "black", size = 5, show.legend = TRUE) +     #sets fill color to sampletype
  scale_shape_manual(values = c(21, 22, 23, 24, 25))+
  scale_fill_manual(values = sample_colors) +
  theme_classic() +                                                      #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 20),                               #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(color = FALSE)+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command

T2 <- plot_tsne_phyloseq(T2ps, tsneT2, color = "Substrate", shape = "Inocula") +
  geom_point(aes(fill = Substrate), color = "black", size = 5, show.legend = TRUE) +     #sets fill color to sampletype
  scale_shape_manual(values = c(21, 22, 23, 24, 25))+
  scale_fill_manual(values = sample_colors) +
  theme_classic() +                                                      #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 20),                               #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(color = FALSE)+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command

T3 <- plot_tsne_phyloseq(T3ps, tsneT3, color = "Substrate", shape = "Inocula") +
  geom_point(aes(fill = Substrate), color = "black", size = 5, show.legend = TRUE) +     #sets fill color to sampletype
  scale_shape_manual(values = c(21, 22, 23, 24, 25))+
  scale_fill_manual(values = sample_colors) +
  theme_classic() +                                                      #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 20),                               #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(color = FALSE)+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command

T4 <- plot_tsne_phyloseq(T4ps, tsneT4, color = "Substrate", shape = "Inocula") +
  geom_point(aes(fill = Substrate), color = "black", size = 5, show.legend = TRUE) +     #sets fill color to sampletype
  scale_shape_manual(values = c(21, 22, 23, 24, 25))+
  scale_fill_manual(values = sample_colors) +
  theme_classic() +                                                      #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 20),                               #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(color = FALSE)+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command

library(ggpubr)

labs <- c("T1","T2", "T3", "T4")

ggarrange(T1, T2, T3, T4, labels = labs,  nrow = 2, ncol = 2, common.legend = TRUE)

#Making a PCoA plot
#ordination
all_pcoa <- ordinate(
  physeq = justbacteria, 
  method = "PCoA", 
  distance = "bray"
)

#plot
plot_ordination(
  physeq = justbacteria,                                                         #phyloseq object
  ordination = all_pcoa)+                                                #ordination
  geom_point(aes(fill = Inocula, shape = Substrate), size = 3) +     #sets fill color to sampletype
  scale_shape_manual(values = c(21, 22, 23, 24))+
  scale_fill_manual(values = sample_colors) +
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
  physeq = justbacteria,                                                         #phyloseq object
  ordination = all_pcoa)+                                                #ordination
  geom_point(aes(fill = Time, shape = Substrate), size = 3) +     #sets fill color to sampletype
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

#plot color by substrate
plot_ordination(
  physeq = justbacteria,                                                         #phyloseq object
  ordination = all_pcoa)+                                                #ordination
  geom_point(aes(fill = Substrate, shape = Inocula), size = 3) +     #sets fill color to sampletype
  scale_shape_manual(values = c(21, 22, 23, 24, 25))+
  scale_fill_manual(values = sample_colors) +
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
  filter(Substrate != "Blank" & Inocula != "Blank")
head(vectors)

ggplot(vectors, aes(x = Time, y = Axis.1))+
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

ggplot(vectors, aes(x = Time, y = Axis.1, fill = Substrate))+
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
  geom_point(aes(fill = Inocula, shape = Substrate), size = 3) +     #sets fill color to sampletype
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



## Time PCOA



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
  geom_point(aes(fill = Inocula, shape = Substrate), size = 3) +     #sets fill color to sampletype
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
  geom_point(aes(fill = Inocula, shape = Substrate), size = 3) +     #sets fill color to sampletype
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
  geom_point(aes(fill = Inocula, shape = Substrate), size = 3) +     #sets fill color to sampletype
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
  geom_point(aes(fill = Inocula, shape = Substrate), size = 3) +     #sets fill color to sampletype
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
TAps <- subset_samples(pps, Substrate == "Terephthalamide")
TPAps <- subset_samples(pps, Substrate == "Terephthalate")
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
  geom_point(aes(fill = Inocula, shape = Time), size = 3) +     #sets fill color to sampletype
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
  geom_point(aes(fill = Inocula, shape = Time), size = 3) +     #sets fill color to sampletype
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
  geom_point(aes(fill = Inocula, shape = Time), size = 3) +     #sets fill color to sampletype
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
  select(Phylum, Substrate, Abundance,Inocula, Time) %>%
  group_by(Phylum, Substrate, Inocula, Time) %>%
  filter(Abundance > 0) %>%
  summarize(
    avg_abundance = sum(Abundance)
  ) %>%
  filter(Substrate != "Blank" & Inocula != "Blank")
head(all)

#4) Save phylum colors
phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
)

#5) Create taxa plot!
ggplot(all)+
  geom_col(mapping = aes(x = Substrate, y = avg_abundance, fill = Phylum), position = "fill", color = "black", show.legend = TRUE)+
  facet_grid(rows = vars(Inocula), cols = vars(Time))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = genus_colors) +
  xlab(NULL)+
  theme_minimal()+
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 10, angle = 0, vjust = 1, hjust = 0.5),
        axis.title.y = element_text(size = 20),
        title = element_text(size = 25))

abundanceplot <- phylumabundance %>%
  select(Phylum, Abundance, Inocula, Time) %>%
  group_by(Phylum, Inocula, Time) %>%
  summarize(
    avg_abundance = mean(Abundance)
  ) %>%
  filter(avg_abundance > 0.10) %>% #filter everything less than 10% relative abundance
  filter(Inocula != "Blank") 
head(abundanceplot)


ggplot(abundanceplot, mapping = aes(x = Time, y = avg_abundance, fill = Phylum), color = "black", show.legend = TRUE)+
  geom_point(position = "fill", shape = 21, size = 3, show.legend = TRUE)+
  facet_grid(rows = vars(Inocula))+
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
  select(Class, Substrate, Abundance,Inocula, Time, Replicate) %>%
  mutate( Class = as.character(Class),
          Class = ifelse(Abundance >= 0.05, Class, "< 0.05%"))%>%
  group_by(Class, Substrate,Inocula, Time) %>%
  summarize(
    avg_abundance = sum(Abundance)
  ) %>%
  filter(Substrate != "Blank" & Inocula != "Blank")
View(all2)

TA <- filter(all2, Substrate == "Terephthalamide")
TPA <- filter(all2, Substrate == "Terephthalate")
C <- filter(all2, Substrate == "Control")

#4) Save phylum colors
phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
)

#5) Create taxa plot!
ggplot(all2)+
  geom_col(mapping = aes(x = Substrate, y = avg_abundance, fill = Class), position = "fill", color = "black", show.legend = TRUE)+
  facet_grid(rows = vars(Inocula), cols = vars(Time))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = phylum_colors) +
  xlab(NULL)+
  theme_minimal()+
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5),
        axis.title.y = element_text(size = 20),
        title = element_text(size = 25))

#2) Convert phyloseq object into a data frame with relative abundance counts
familyabundance <- justbacteria %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Family) 
head(familyabundance)

#Summarize the most abundant taxa, filter to only include taxa with > 1% abundance
all3 <- familyabundance %>%
  select(Family, Substrate, Abundance, Inocula, Time) %>%
  mutate( Family = as.character(Family),
          Family = ifelse(Abundance >= 0.025, Family, "< 0.025%"))%>%
  group_by(Family, Substrate,Inocula, Time) %>%
  summarize(
    avg_abundance = sum(Abundance)
  ) %>%
  filter(Substrate != "Blank" & Inocula != "Blank")
head(all3)

#4) Save colors
genus_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
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
  "white",   "dodgerblue",   "firebrick", "yellowgreen", "magenta", "black", "lightpink",
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
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
  "white",   "dodgerblue",   "firebrick", "yellowgreen", "magenta", "black", "lightpink",
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
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
ggplot(all3B)+
  geom_col(mapping = aes(x = Time, y = avg_abundance, fill = Family), color = "black", position = "fill", show.legend = TRUE)+
  facet_grid(cols = vars(Substrate))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = genus_colors) +
  xlab(NULL)+
  theme_minimal()+
  theme(strip.text = element_text(face = "bold", size = 20))+ 
  theme(legend.text = element_text(size=15), legend.position = "right",
    axis.text.y.left = element_text(size = 15),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 1, hjust = 0.5),
        axis.title.y = element_text(size = 20),
        title = element_text(size = 25))

all3B <- filter(all3, Substrate != "Control")


#Summarize the most abundant taxa, filter to only include taxa with > 1% abundance
all4 <- genusabundance %>%
  select(Genus, Substrate, Abundance, Inocula, Time, Replicate) %>%
  filter(Substrate != "Blank" & Inocula != "Blank") %>%
  group_by(Genus, Substrate, Inocula, Time, Replicate) %>%
  summarize(
    avg_abundance = mean(Abundance)
  ) %>%
  mutate(Genus = as.character(Genus),
         Genus = ifelse(mean(avg_abundance) >= 0.025, Genus, "< 0.025%"),
         Replicate = ifelse(Replicate == "R2A", "R2", Replicate),
         Replicate = ifelse(Replicate == "R2B", "R3", Replicate)) %>%
  filter(Replicate == "R1" | Replicate == "R2" | Replicate == "R3") %>%
  group_by(Genus, Substrate, Inocula, Time, Replicate) %>%
  summarize(
    avg_abundance = sum(avg_abundance)
  )
head(all4)

TA <- filter(all4, Substrate == "Terephthalamide")
TPA <- filter(all4, Substrate == "Terephthalate")
C <- filter(all4, Substrate == "Control")


#5) Create taxa plot!
ggplot(all4)+
  geom_col(mapping = aes(x = Substrate, y = avg_abundance, fill = Genus), color = "black", position = "fill", show.legend = FALSE)+
  facet_grid(cols = vars(Time))+
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
  select(Genus, Substrate, Abundance, Inocula, Time, Replicate) %>%
  filter(Substrate != "Blank" & Inocula != "Blank") %>%
  filter(Abundance > 0) %>%
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

