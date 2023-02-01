#install DESeq2: https://bioconductor.org/packages/release/bioc/html/DESeq2.html
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("DESeq2")

#packages used
library(tidyverse)
#install.packages("matrixStats")
library(matrixStats)
library(DESeq2)
library(phyloseq)


#phyloseq object
phyloseq_object_all <- readRDS("/home/lgschaer/old/Plastic_Deg/TA_TPA_Enrichment_09062021/ps_unrarefied_09072021.rds")
head(sample_data(phyloseq_object_all))

#filter out eukaryotes and mitochondria
#will also remove inocula samples and subset to only include the final transfer
justbacteria <- phyloseq_object_all %>%
  subset_taxa(
    Kingdom == "Bacteria" &                   #only bacteria
      Family  != "Mitochondria" &             #filter out mitochondria
      Class   != "Chloroplast"                #filter out chloroplasts
  ) %>%
  subset_samples(Substrate != "Inocula") %>%
  subset_samples(Transfer == "T4")
justbacteria

#saveRDS file
saveRDS(justbacteria, file = "/home/lgschaer/old/Plastic_Deg/TA_TPA_Enrichment_09062021/tatpa_rarefied_nochloroplasts.rds")
head(sample_data(justbacteria))

#add a pseudo count of one to the OTU table
justbacteria@otu_table <- as.matrix(justbacteria@otu_table)+1

#subset phyloseq object to make sample comparison categories, Media_Carbon = DCPET_BH or TPA_BH
sample_info <- as.data.frame(sample_data(justbacteria))
unique(sample_info$Substrate)

ta_tpa <- justbacteria %>% subset_samples(Substrate == "TA" | Substrate == "TPA")
ta_tpa2 <- prune_taxa(taxa_sums(ta_tpa) > 0, ta_tpa)
ta_c <- subset_samples(justbacteria, (Substrate == "TA" | Substrate == "Control"))
ta_c2 <- prune_taxa(taxa_sums(ta_c) > 0, ta_c) 
tpa_c <- subset_samples(justbacteria, (Substrate == "TPA" | Substrate == "Control"))
tpa_c2 <- prune_taxa(taxa_sums(tpa_c) > 0, tpa_c) 

## check subsetting
unique(as.data.frame(sample_data(ta_tpa))$Substrate)
unique(as.data.frame(sample_data(ta_c))$Substrate)
unique(as.data.frame(sample_data(tpa_c))$Substrate)

## make sure there are no zero counts
any(taxa_sums(ta_tpa2)==0)
any(taxa_sums(ta_c2)==0)
any(taxa_sums(tpa_c2)==0)

# Convert phyloseq object to deseq2 format
ps.ta.tpa <- phyloseq_to_deseq2(ta_tpa2, ~ Substrate)
ps.ta.c <- phyloseq_to_deseq2(ta_c2, ~ Substrate)
ps.tpa.c <- phyloseq_to_deseq2(tpa_c2, ~ Substrate)

ps.ta.tpa$Substrate<-relevel(ps.ta.tpa$Substrate, "TA")
ps.ta.c$Substrate<-relevel(ps.ta.c$Substrate, "TA")
ps.tpa.c$Substrate<-relevel(ps.tpa.c$Substrate, "TPA")

# Run DESeq2 analysis (all taxa at once!)
dds_ps.ta.tpa <- DESeq(ps.ta.tpa)
dds_ps.ta.c <- DESeq(ps.ta.c)
dds_ps.tpa.c <- DESeq(ps.tpa.c)

# Investigate results
resultsNames(dds_ps.ta.tpa)
resultsNames(dds_ps.ta.c)
resultsNames(dds_ps.tpa.c)

# Put DESeq output into data frames
tatpa <- as.data.frame(results(dds_ps.ta.tpa, contrast=c("Substrate","TPA","TA"))) %>% mutate(Comparison="Terephthalate vs. Terephthalamide") %>% rownames_to_column(var = "taxon")
tac <- as.data.frame(results(dds_ps.ta.c, contrast=c("Substrate","Control","TA"))) %>% mutate(Comparison="Control vs. Terephthalamide") %>% rownames_to_column(var = "taxon")
tpac <- as.data.frame(results(dds_ps.tpa.c, contrast=c("Substrate","Control","TPA"))) %>% mutate(Comparison="Control vs. Terephthalate") %>% rownames_to_column(var = "taxon")
head(tpac)

# Add taxonomy to DESeq output
taxps.tatpa <- as.data.frame(tax_table(ta_tpa)) %>% rownames_to_column(var = "taxon") %>% full_join(tatpa)
taxps.tac <- as.data.frame(tax_table(ta_c)) %>% rownames_to_column(var = "taxon") %>% full_join(tac)
taxps.tpac <- as.data.frame(tax_table(tpa_c)) %>% rownames_to_column(var = "taxon") %>% full_join(tpac)
head(taxps.tpac)

##check dimensions
dim(tatpa)
dim(taxps.tatpa)

# Join everything together
enriched <- taxps.tatpa %>% 
  full_join(taxps.tac) %>%
  full_join(taxps.tpac) %>%
  filter(!is.na(padj)) %>%
  mutate(
    threshold = ifelse(padj <= 0.001 & abs(log2FoldChange) >= 2, "Enriched", "Not_Enriched"),
    Enriched_Genus = ifelse(threshold == "Enriched", as.character(Genus), "Not Enriched"),
    Enriched_Genus = ifelse(is.na(Enriched_Genus), "Unclassified Genus", Enriched_Genus)
  ) 
head(enriched)
#View(enriched)

# Are any ASVs enriched?
any(enriched$threshold == "Enriched")
sum(enriched$threshold == "Enriched")

# Save a csv of the results
write.csv(enriched,"/home/lgschaer/old/Plastic_Deg/TA_TPA_Enrichment_09062021/dsq_out/enriched_w_taxonomy_01052023.csv")

enriched_counts <- enriched %>%
  mutate(Category = ifelse(log2FoldChange > 2, "Right", "not_sig"),
         Category = ifelse(log2FoldChange < -2 & Category == "not_sig", "Left", Category)) %>%
  group_by(Comparison, Category) %>%
  summarise(Count = n()) %>%
  filter(Category != "not_sig")
enriched_counts

#any(is.na(enriched_w_tax))
#View(enriched_counts)

# Re-Order data to organize legend
sort(unique(enriched$Enriched_Genus))
length(unique(enriched$Enriched_Genus))

#enriched_w_tax$Enriched_Genus <- factor(enriched_w_tax$Enriched_Genus, 
 #                                       levels = c("Achromobacter",      "Aminobacter",        "Ancylobacter",       "Aquamicrobium",      "Bauldia",            "Bosea",
  #                                                 "Brevundimonas",      "Bryobacter",         "Candidimonas",       "Chelatococcus",     
   #                                                "Chitinophaga",       "Cryobacterium",      "Devosia",            "Hydrogenophaga",     "Hyphomonas",         "Legionella",         
    #                                               "Luteimonas",         "Mesorhizobium",      "Microbacterium",     "Millisia",           "Orrella",            "Paramesorhizobium",
     #                                              "Parapedobacter",     "Parapusillimonas",   "Parvibaculum",       "Pedomicrobium",      "Pelagibacterium",   
      #                                             "Persicitalea",       "Planktosalinus",     "Pseudaminobacter",   "Pseudolabrys",       "Pseudomonas",        "Pseudoxanthomonas", 
       #                                            "Pusillimonas",       "Rhodobacter",        "Rhodococcus",        "Shinella",           "SN8",                "Sphingobacterium",  
        #                                           "Tepidimonas",        "Thermomonas",        "Tianweitania",       "Variovorax",         "Verticiella",        "Youhaiella",   
         #                                          "Unclassified Genus", "Not Enriched"))
# Save color palette

cl <- colors(distinct = TRUE)
set.seed(15887) # to set random generator seed
colors11 <- sample(cl, length(unique(enriched$Enriched_Genus)))

#colors11 <- c(
# "magenta","lightpink","yellow","orange","cyan","blue","mediumseagreen","firebrick","lightblue",
#"magenta","lightpink","yellow","orange","cyan","blue","mediumseagreen","firebrick","lightblue",
#  "magenta","lightpink","yellow","orange","cyan","blue","mediumseagreen","firebrick","lightblue",
# "magenta","lightpink","yellow","orange","cyan","blue","mediumseagreen","firebrick","lightblue",
#"magenta","lightpink","yellow","orange","cyan","blue","mediumseagreen","gray63","white"
#)

shapes <- c(
  21, 22, 23, 24,
  21, 22, 23, 24,
  21, 22, 23, 24,
  21, 22, 23, 24,
  21, 22, 23, 24,
  21, 22, 23, 24, 
  21, 22, 23, 24, 
  21, 22, 23, 24,
  21, 22, 23, 24,
  21, 22, 23, 24,
  22, 21, 21
)
length(shapes)

#colors11 <- c(
 # "orangered",      "purple",        "green",           "cyan",          "orange",        "khaki4",             "mediumslateblue",
  #"mediumpurple1",  "darkmagenta",   "darkgreen",       "wheat2",        "yellow",        "lawngreen",          "plum",  
#  "royalblue",      "magenta",       "mediumseagreen",  "palegoldenrod", "grey47",        "chocolate4",         "darkorange3",        
 # "lightblue",      "firebrick",     "yellowgreen",     "turquoise3",    "purple4",       "blue",               "red",            
  #"lightcyan",       "coral1",       "cyan",            "goldenrod",     "yellowgreen",   "turquoise3",    "purple4",       "blue",               "red",            
#  "lightcyan",       "coral1",       "cyan",            "goldenrod",     "black",         "white"   
#) 
#length(colors11)

#colors11 <- c(
 # "palegoldenrod","palegoldenrod","palegoldenrod","palegoldenrod",
  #"orange",  "orange",  "orange",  "orange",  
#  "firebrick","firebrick","firebrick","firebrick",
 # "pink","lightpink","lightpink","lightpink",
  #"purple","purple","purple","purple",
#  "darkblue","darkblue","darkblue","darkblue",
 # "darkcyan", "darkcyan", "darkcyan", "darkcyan", 
  #"lightblue","lightblue","lightblue","lightblue",
#  "olivedrab2", "olivedrab2", "olivedrab2", "olivedrab2", 
 # "darkgreen", "darkgreen", "darkgreen", "darkgreen", 
  #"grey77","grey77","white"
#)

#volcano plot:
ggplot(data=enriched, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(fill=Enriched_Genus), shape = 21, color = "black", size=6) +
  facet_grid(cols = vars(Comparison))+
  #scale_fill_manual(values=colors11) +
  #scale_shape_manual(values=shapes) +
  labs(x = "log2 fold change", 
       y = "-log10 p-value") +
  theme_classic(base_size = 14)+
  geom_hline(yintercept = -log10(0.001), colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = 2, colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = -2, colour="#990000", linetype="dashed")+ 
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 25, face = "bold", angle = 0),
        legend.position = "bottom",
        title = element_text(size = 18))#+
  #guides(fill = guide_legend(override.aes = list(shape = shapes)))

# another way to plot
library(ggh4x)

unique(enriched$Comparison)

enriched_filt <- enriched %>%
  filter(threshold == "Enriched") %>%
  mutate(Genus2 = ifelse(is.na(Genus), "", Genus),
         Species2 = ifelse(is.na(Species), "", Species)) %>%
  unite(Taxonomy, Genus2, Species2, sep = " ") %>%
  mutate(Taxonomy = ifelse(Taxonomy == " ", "Unclassified Organism", Taxonomy),
    Where_Enriched = case_when(
    (Comparison == "Terephthalate vs. Terephthalamide" & log2FoldChange > 2) ~ "Terephthalamide",
    (Comparison == "Terephthalate vs. Terephthalamide" & log2FoldChange < -2) ~ "Terephthalate",
    (Comparison == "Control vs. Terephthalamide" & log2FoldChange > 2) ~ "Terephthalamide",
    (Comparison == "Control vs. Terephthalamide" & log2FoldChange < -2) ~ "Control",
    (Comparison == "Control vs. Terephthalate" & log2FoldChange > 2) ~ "Terephthalate",
    (Comparison == "Control vs. Terephthalate" & log2FoldChange < -2) ~ "Control"
  )) %>%
  group_by(Comparison, Enriched_Genus, Where_Enriched) #%>%
 # summarise(Count = n())
head(enriched_filt)
#View(enriched_filt)

colors <- c("olivedrab", "lightgoldenrod")

max <- 0.001
min <- round(min(enriched_filt$padj), digits = 3)
#min <- 0.5

ggplot(data=enriched_filt, aes(y=Taxonomy, x=log2FoldChange)) +
  geom_point(aes(fill=padj, shape = Where_Enriched), color = "black", size = 6) +
  facet_nested(cols = vars(Comparison), rows = vars(Class), space = "free_y", scales = "free_y")+
  scale_fill_gradientn(colors = colors, breaks = seq(min, max, by = 0.0001),#)+#,
                       limits = c(min, max), labels = as.character(seq(min, max, by = 0.0001))) +
  scale_shape_manual(values=c(21, 22, 23)) +
  xlab("log2 fold change")+ 
  ylab("Significantly Enriched Genera") +
  theme_linedraw(base_size = 14)+
  theme(axis.text.x = element_text(size = 20, angle = 0, vjust = 1, hjust = 0.5),
        axis.text.y = element_text(size = 20, angle = 0, hjust = 1, vjust = 0.5),
        axis.title.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        legend.text = element_text(size = 18),
        #legend.title = element_blank(),
        strip.text.x = element_text(size = 25, face = "bold", angle = 0),
        strip.text.y = element_text(size = 18, face = "bold", angle = 0),
        legend.position = "bottom",
        title = element_text(size = 18)) +
  guides(fill = guide_colourbar(barwidth = 40, barheight = 1, title = "Adjusted P-value"))

# same figure filtered even more (genus_list from phyloseq script)
genus_list <- c("Arthrobacter",       "Bradyrhizobium",     "Niabella",           "Eoetvoesia",        
                "Hydrogenophaga",     "Ramlibacter",        "Variovorax",         "Pseudomonas",       
                "Hoeflea",            "Ellin6055",          "Pseudoxanthomonas",  "Terrimonas",        
                "Taibaiella",         "UTBCD1",             "Rhodococcus",        "Bordetella",        
                "Noviherbaspirillum", "Micromonospora",     "Pseudorhodoplanes",  "Thermomonas",       
                "Acidovorax",         "Paenibacillus",      "Mesorhizobium",      "Luteimonas",        
                "Rhodococcus",        "Leadbetterella",     "SWB02",              "Mesorhizobium",     
                "Bordetella",         "Bosea",              "Aliihoeflea",        "Verticiella",       
                "Aminobacter",        "Pseudomonas",        "Sphingobacterium",   "Hydrogenophaga",    
                "Flavobacterium",     "Achromobacter",      "Thermomonas",        "Bacillus",          
                "Micromonospora",     "Pseudoxanthomonas")

enriched_filt2 <- enriched_filt %>%
  filter(Genus %in% genus_list)

colors <- c("olivedrab", "lightgoldenrod")

max <- 0.001
min <- round(min(enriched_filt$padj), digits = 3)
#min <- 0.5

ggplot(data=enriched_filt2, aes(y=Taxonomy, x=log2FoldChange)) +
  geom_point(aes(fill=padj, shape = Where_Enriched), color = "black", size = 6) +
  facet_nested(cols = vars(Comparison), rows = vars(Class), space = "free_y", scales = "free_y")+
  scale_fill_gradientn(colors = colors, breaks = seq(min, max, by = 0.0001),#)+#,
                       limits = c(min, max), labels = as.character(seq(min, max, by = 0.0001))) +
  scale_shape_manual(values=c(21, 22, 23)) +
  xlab("log2 fold change")+ 
  ylab("Significantly Enriched Genera") +
  theme_linedraw(base_size = 14)+
  theme(axis.text.x = element_text(size = 20, angle = 0, vjust = 1, hjust = 0.5),
        axis.text.y = element_text(size = 20, angle = 0, hjust = 1, vjust = 0.5),
        axis.title.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        legend.text = element_text(size = 18),
        #legend.title = element_blank(),
        strip.text.x = element_text(size = 25, face = "bold", angle = 0),
        strip.text.y = element_text(size = 18, face = "bold", angle = 0),
        legend.position = "bottom",
        title = element_text(size = 18)) +
  guides(fill = guide_colourbar(barwidth = 40, barheight = 1, title = "Adjusted P-value"))


#### Looking at alternate pairings of the data


#phyloseq object
phyloseq_object_all <- readRDS("/home/lgschaer/old/Plastic_Deg/TA_TPA_Enrichment_09062021/ps_unrarefied_09072021.rds")
head(sample_data(phyloseq_object_all))

new_sample_data <- as_tibble(sample_data(phyloseq_object_all)) %>%
  mutate(isTA = ifelse(Substrate == "TA", "TA", "notTA"),
         isTPA = ifelse(Substrate == "TPA", "TPA", "notTPA"),
         RowNames = SampleID) %>%
  column_to_rownames(var = "RowNames")
head(new_sample_data)

phyloseq_object_all@sam_data <- sample_data(new_sample_data)
head(sample_data(phyloseq_object_all))

#filter out eukaryotes and mitochondria
#will also remove inocula samples and subset to only include the final transfer
justbacteria <- phyloseq_object_all %>%
  subset_taxa(
    Kingdom == "Bacteria" &                   #only bacteria
      Family  != "Mitochondria" &             #filter out mitochondria
      Class   != "Chloroplast"                #filter out chloroplasts
  ) %>%
  subset_samples(Substrate != "Inocula") #%>%
  #subset_samples(Transfer == "T4")
justbacteria

#saveRDS file
saveRDS(justbacteria, file = "/home/lgschaer/old/Plastic_Deg/TA_TPA_Enrichment_09062021/ps_nomito_nochloro_noinocula_T4.rds")
head(sample_data(justbacteria))

#add a pseudo count of one to the OTU table
justbacteria@otu_table <- as.matrix(justbacteria@otu_table)+1

#subset phyloseq object to make sample comparison categories, Media_Carbon = DCPET_BH or TPA_BH
sample_info <- as.data.frame(sample_data(justbacteria))
unique(sample_info$isTA)
unique(sample_info$isTPA)

## make sure there are no zero counts
any(taxa_sums(justbacteria)==0)

# Convert phyloseq object to deseq2 format
ps.isTA <- phyloseq_to_deseq2(justbacteria, ~ isTA)
ps.isTPA <- phyloseq_to_deseq2(justbacteria, ~ isTPA)

# Run DESeq2 analysis (all taxa at once!)
dds_ps.isTA <- DESeq(ps.isTA)
dds_ps.isTPA <- DESeq(ps.isTPA)

# Investigate results
resultsNames(dds_ps.isTA)
resultsNames(dds_ps.isTPA)

# Put DESeq output into data frames
isTA <- as.data.frame(results(dds_ps.isTA, contrast=c("isTA","notTA","TA"))) %>% mutate(Comparison="Terephthalamide Enriched") %>% rownames_to_column(var = "taxon")
isTPA <- as.data.frame(results(dds_ps.isTPA, contrast=c("isTPA","notTPA","TPA"))) %>% mutate(Comparison="Terephthalate Enriched") %>% rownames_to_column(var = "taxon")
head(isTPA)

# Add taxonomy to DESeq output
taxps.isTA <- as.data.frame(tax_table(justbacteria)) %>% rownames_to_column(var = "taxon") %>% full_join(isTA)
taxps.isTPA <- as.data.frame(tax_table(justbacteria)) %>% rownames_to_column(var = "taxon") %>% full_join(isTPA)
head(taxps.isTPA)

##check dimensions
dim(isTPA)
dim(taxps.isTPA)

# Join everything together
enriched <- taxps.isTA %>% 
  full_join(taxps.isTPA) %>%
  filter(!is.na(padj)) %>%
  mutate(
    threshold = ifelse(padj <= 0.001 & abs(log2FoldChange) >= 2, "Enriched", "Not_Enriched"),
    Enriched_Genus = ifelse(threshold == "Enriched", as.character(Genus), "Not Enriched"),
    Enriched_Genus = ifelse(is.na(Enriched_Genus), "Unclassified Genus", Enriched_Genus)
  ) 
head(enriched)
#View(enriched)

# Are any ASVs enriched?
any(enriched$threshold == "Enriched")
sum(enriched$threshold == "Enriched")

# Save a csv of the results
write.csv(enriched,"/home/lgschaer/old/Plastic_Deg/TA_TPA_Enrichment_09062021/dsq_out/enriched_w_taxonomy_01052023.csv")

enriched_counts <- enriched %>%
  mutate(Category = ifelse(log2FoldChange > 2, "Right", "not_sig"),
         Category = ifelse(log2FoldChange < -2, "Left", Category)) %>%
  group_by(Comparison, Category) %>%
  summarise(Count = n()) %>%
  filter(Category != "not_sig")
enriched_counts

#any(is.na(enriched_w_tax))
#View(enriched_counts)

# Re-Order data to organize legend
sort(unique(enriched$Enriched_Genus))
length(unique(enriched$Enriched_Genus))

#enriched_w_tax$Enriched_Genus <- factor(enriched_w_tax$Enriched_Genus, 
#                                       levels = c("Achromobacter",      "Aminobacter",        "Ancylobacter",       "Aquamicrobium",      "Bauldia",            "Bosea",
#                                                 "Brevundimonas",      "Bryobacter",         "Candidimonas",       "Chelatococcus",     
#                                                "Chitinophaga",       "Cryobacterium",      "Devosia",            "Hydrogenophaga",     "Hyphomonas",         "Legionella",         
#                                               "Luteimonas",         "Mesorhizobium",      "Microbacterium",     "Millisia",           "Orrella",            "Paramesorhizobium",
#                                              "Parapedobacter",     "Parapusillimonas",   "Parvibaculum",       "Pedomicrobium",      "Pelagibacterium",   
#                                             "Persicitalea",       "Planktosalinus",     "Pseudaminobacter",   "Pseudolabrys",       "Pseudomonas",        "Pseudoxanthomonas", 
#                                            "Pusillimonas",       "Rhodobacter",        "Rhodococcus",        "Shinella",           "SN8",                "Sphingobacterium",  
#                                           "Tepidimonas",        "Thermomonas",        "Tianweitania",       "Variovorax",         "Verticiella",        "Youhaiella",   
#                                          "Unclassified Genus", "Not Enriched"))
# Save color palette

cl <- colors(distinct = TRUE)
set.seed(15887) # to set random generator seed
colors11 <- sample(cl, length(unique(enriched$Enriched_Genus)))

#colors11 <- c(
# "magenta","lightpink","yellow","orange","cyan","blue","mediumseagreen","firebrick","lightblue",
#"magenta","lightpink","yellow","orange","cyan","blue","mediumseagreen","firebrick","lightblue",
#  "magenta","lightpink","yellow","orange","cyan","blue","mediumseagreen","firebrick","lightblue",
# "magenta","lightpink","yellow","orange","cyan","blue","mediumseagreen","firebrick","lightblue",
#"magenta","lightpink","yellow","orange","cyan","blue","mediumseagreen","gray63","white"
#)

shapes <- c(
  21, 22, 23, 24,
  21, 22, 23, 24,
  21, 22, 23, 24,
  21, 22, 23, 24,
  21, 22, 23, 24,
  21, 22, 23, 24, 
  21, 22, 23, 24, 
  21, 22, 23, 24,
  21, 22, 23, 24,
  21, 22, 23, 24,
  22, 21, 21
)
length(shapes)

#colors11 <- c(
# "orangered",      "purple",        "green",           "cyan",          "orange",        "khaki4",             "mediumslateblue",
#"mediumpurple1",  "darkmagenta",   "darkgreen",       "wheat2",        "yellow",        "lawngreen",          "plum",  
#  "royalblue",      "magenta",       "mediumseagreen",  "palegoldenrod", "grey47",        "chocolate4",         "darkorange3",        
# "lightblue",      "firebrick",     "yellowgreen",     "turquoise3",    "purple4",       "blue",               "red",            
#"lightcyan",       "coral1",       "cyan",            "goldenrod",     "yellowgreen",   "turquoise3",    "purple4",       "blue",               "red",            
#  "lightcyan",       "coral1",       "cyan",            "goldenrod",     "black",         "white"   
#) 
#length(colors11)

#colors11 <- c(
# "palegoldenrod","palegoldenrod","palegoldenrod","palegoldenrod",
#"orange",  "orange",  "orange",  "orange",  
#  "firebrick","firebrick","firebrick","firebrick",
# "pink","lightpink","lightpink","lightpink",
#"purple","purple","purple","purple",
#  "darkblue","darkblue","darkblue","darkblue",
# "darkcyan", "darkcyan", "darkcyan", "darkcyan", 
#"lightblue","lightblue","lightblue","lightblue",
#  "olivedrab2", "olivedrab2", "olivedrab2", "olivedrab2", 
# "darkgreen", "darkgreen", "darkgreen", "darkgreen", 
#"grey77","grey77","white"
#)

#volcano plot:
ggplot(data=enriched, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(fill=Enriched_Genus), shape = 21, color = "black", size=6) +
  facet_grid(cols = vars(Comparison))+
  #scale_fill_manual(values=colors11) +
  #scale_shape_manual(values=shapes) +
  labs(x = "log2 fold change", 
       y = "-log10 p-value") +
  theme_classic(base_size = 14)+
  geom_hline(yintercept = -log10(0.001), colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = 2, colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = -2, colour="#990000", linetype="dashed")+ 
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 25, face = "bold", angle = 0),
        legend.position = "bottom",
        title = element_text(size = 18))#+
#guides(fill = guide_legend(override.aes = list(shape = shapes)))

# another way to plot
library(ggh4x)

unique(enriched$Comparison)

enriched_summary <- enriched %>%
  filter(threshold == "Enriched") %>%
  mutate(Genus2 = ifelse(is.na(Genus), "", Genus),
         Species2 = ifelse(is.na(Species), "", Species)) %>%
  unite(Taxonomy, Genus2, Species2, sep = " ") %>%
  mutate(Taxonomy = ifelse(Taxonomy == " ", "Unclassified Organism", Taxonomy),
         Where_Enriched = case_when(
           (Comparison == "Terephthalate Enriched" & log2FoldChange > 2) ~ "Terephthalamide+Controls",
           (Comparison == "Terephthalate Enriched" & log2FoldChange < -2) ~ "Terephthalate",
           (Comparison == "Terephthalamide Enriched" & log2FoldChange > 2) ~ "Terephthalate+Controls",
           (Comparison == "Terephthalamide Enriched" & log2FoldChange < -2) ~ "Terephthalamide")) %>%
  select(Where_Enriched, Phylum, Class, Order, Family, Genus, Species, Taxonomy, log2FoldChange, padj)
head(enriched_summary)
#View(enriched_summary)

write_csv(enriched_summary, "/home/lgschaer/old/Plastic_Deg/TA_TPA_Enrichment_09062021/deseq_enriched_w_tax_01062023.csv")

sum(enriched_summary$Where_Enriched=="Terephthalate")
sum(enriched_summary$Where_Enriched=="Terephthalamide")

dim(enriched_summary)


enriched_filt <- enriched %>%
  filter(threshold == "Enriched") %>%
  mutate(Genus2 = ifelse(is.na(Genus), "", Genus),
         Species2 = ifelse(is.na(Species), "", Species)) %>%
  unite(Taxonomy, Genus2, Species2, sep = " ") %>%
  mutate(Taxonomy = ifelse(Taxonomy == " ", "Unclassified Organism", Taxonomy),
         Where_Enriched = case_when(
           (Comparison == "Terephthalate Enriched" & log2FoldChange > 2) ~ "All Other Samples",
           (Comparison == "Terephthalate Enriched" & log2FoldChange < -2) ~ "Terephthalate",
           (Comparison == "Terephthalamide Enriched" & log2FoldChange > 2) ~ "All Other Samples",
           (Comparison == "Terephthalamide Enriched" & log2FoldChange < -2) ~ "Terephthalamide"),
         Where_Enriched = factor(Where_Enriched, levels = c("Terephthalamide", "Terephthalate", "All Other Samples")))
head(enriched_filt)
#View(enriched_filt)

ggplot(data=enriched_filt, aes(y=Taxonomy, x=log2FoldChange)) +
  facet_nested(cols = vars(Comparison), rows = vars(Class), space = "free_y", scales = "free_y")+
  geom_col(aes(fill = Where_Enriched, y=Taxonomy, x=log2FoldChange), color = "black")+
  scale_fill_manual(values = c("maroon", "olivedrab4", "white"))+
  xlab("log2 fold change")+ 
  ylab("Significantly Enriched Taxa") +
  theme_linedraw(base_size = 14)+
  theme(axis.text.x = element_text(size = 30, angle = 0, vjust = 1, hjust = 0.5),
        axis.text.y = element_text(size = 30, angle = 0, hjust = 1, vjust = 0.5),
        axis.title.y = element_text(size = 30),
        axis.title.x = element_text(size = 25),
        legend.text = element_text(size = 28),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 40, face = "bold", angle = 0),
        strip.text.y = element_text(size = 40, face = "bold", angle = 0),
        legend.position = "bottom",
        title = element_text(size = 18))


# same figure filtered even more (genus_list from phyloseq script)
genus_list <- c("Arthrobacter",       "Bradyrhizobium",     "Niabella",           "Eoetvoesia",        
                "Hydrogenophaga",     "Ramlibacter",        "Variovorax",         "Pseudomonas",       
                "Hoeflea",            "Ellin6055",          "Pseudoxanthomonas",  "Terrimonas",        
                "Taibaiella",         "UTBCD1",             "Rhodococcus",        "Bordetella",        
                "Noviherbaspirillum", "Micromonospora",     "Pseudorhodoplanes",  "Thermomonas",       
                "Acidovorax",         "Paenibacillus",      "Mesorhizobium",      "Luteimonas",        
                "Rhodococcus",        "Leadbetterella",     "SWB02",              "Mesorhizobium",     
                "Bordetella",         "Bosea",              "Aliihoeflea",        "Verticiella",       
                "Aminobacter",        "Pseudomonas",        "Sphingobacterium",   "Hydrogenophaga",    
                "Flavobacterium",     "Achromobacter",      "Thermomonas",        "Bacillus",          
                "Micromonospora",     "Pseudoxanthomonas")

enriched_filt2 <- enriched_filt %>%
  filter(Genus %in% genus_list)
head(enriched_filt2)

ggplot(data=enriched_filt2, aes(y=Taxonomy, x=log2FoldChange)) +
  facet_nested(cols = vars(Comparison), rows = vars(Class), space = "free_y", scales = "free_y")+
  geom_col(aes(fill = Where_Enriched, y=Taxonomy, x=log2FoldChange), color = "black")+
  scale_fill_manual(values = c("maroon", "olivedrab4", "white"))+
  xlab("log2 fold change")+ 
  ylab("Significantly Enriched Taxa") +
  theme_linedraw(base_size = 14)+
  theme(axis.text.x = element_text(size = 30, angle = 0, vjust = 1, hjust = 0.5),
        axis.text.y = element_text(size = 30, angle = 0, hjust = 1, vjust = 0.5),
        axis.title.y = element_text(size = 30),
        axis.title.x = element_text(size = 25),
        legend.text = element_text(size = 28),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 40, face = "bold", angle = 0),
        strip.text.y = element_text(size = 40, face = "bold", angle = 0),
        legend.position = "bottom",
        title = element_text(size = 18))

