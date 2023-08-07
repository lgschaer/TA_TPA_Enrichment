#load packages
library(csv)
library(tidyverse)
library(ggh4x)

#load hplc data
hplc <- read_csv("/home/lgschaer/old/Plastic_Deg/TA_TPA_Enrichment_09062021/HPLC_TPA_TA_Experiment.csv")
head(hplc)


hplc2 <- hplc %>%
  pivot_longer(cols = c("TA_g_L", "TPA_g_L", "TPAMA_g_L"), names_to = "Compound", values_to = "Concentration_g_L") %>%
  separate(SampleID, into = c(NA, NA, "Replicate", NA), sep = " ", remove = FALSE) %>%
  separate(Compound, into = c("Compound", NA), sep = "_") %>%
  select(SampleID, Substrate, Environment, Location, Transfer, Replicate, Compound, Concentration_g_L)
head(hplc2)
#View(hplc2)

tpa <- hplc2 %>% filter(Substrate == "Terephthalate")
unique(tpa$Substrate)

ta <- hplc2 %>% filter(Substrate == "Terephthalamide") %>%
  mutate(Concentration_g_L = Concentration_g_L*14)
unique(ta$Substrate)

colors <- c("maroon4", "olivedrab4", "lightblue")

ggplot(ta, aes(x = Transfer, y = Concentration_g_L))+
  facet_nested(cols = vars(Location, Replicate), rows = vars(Substrate), scales = "free_x", space = "free")+
  geom_col(aes(fill = Compound), color = "black")+
  scale_fill_manual(values = colors)+
  theme_bw(base_size = 14)+
  theme(axis.text.y.left = element_text(size = 10),
        axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "bottom",
        title = element_text(size = 18))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))


ggplot(tpa, aes(x = SampleID, y = Concentration_g_L))+
  facet_nested(cols = vars(Location, Transfer), rows = vars(Substrate), scales = "free_x", space = "free")+
  geom_col(aes(fill = Compound), color = "black")+
  scale_fill_manual(values = colors)+
  theme_bw(base_size = 14)+
  theme(axis.text.y.left = element_text(size = 10),
        axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "bottom",
        title = element_text(size = 18))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))



### Terephthalamide as a percent
head(ta)
View(ta)

ta2 <- ta %>%
  group_by(Substrate, Location, Transfer, Compound) %>%
  mutate(n = 1,
         Percent = sum(Concentration_g_L)/1) %>%
  summarise(
    n = sum(n),
    sum = sum(Percent),
    Percent = sum/n
  )
head(ta2)  
#View(ta2)

ggplot(ta2, aes(x = Transfer, y = Percent))+
  facet_nested(cols = vars(Location), rows = vars(Substrate), scales = "free_x", space = "free")+
  geom_col(aes(fill = Compound), color = "black", position = "fill")+
  scale_fill_manual(values = colors)+
  theme_bw(base_size = 14)+
  theme(axis.text.x = element_text(size = 18, angle = 0, vjust = 0.5, hjust = 0.5, color = "black"),
        axis.title.y = element_text(size = 20),
        strip.text = element_text(size = 22, face = "bold", color = "black"),
        strip.background = element_rect(fill = "white"),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        legend.position = "bottom",
        title = element_text(size = 18))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))
  
  
##########################
#load packages
library(tidyverse)
library(csv)

#load data terephthalamide_HPLC.csv
ta_hplc <- as.csv("/home/lgschaer/old/Plastic_Deg/TA_TPA_Enrichment_09062021/terephthalamide_HPLC.csv", row.names = 1, header = TRUE, sep = ",", check.names = TRUE, stringsAsFactors = TRUE)
head(ta_hplc)

unique(ta_hplc$Environment)

#plot degradation
colors <- c("maroon")

ta_hplc2 <- ta_hplc %>%
  filter(Type == "HPLC") %>%
  group_by(Substrate, Environment) %>%
  summarise(minTA = min(Terephthalamide_g_L),
            maxTA = max(Terephthalamide_g_L),
            Terephthalamide_g_L = mean(Terephthalamide_g_L)) %>%
  mutate(
    Environment = case_when(
      Environment == "Starting" ~ "Starting\nConc.",
      Environment == "Compost" ~ "Compost",
      Environment == "Hydrocarbon Seep" ~ "Hydrocarbon\nSeep",
      Environment == "Lake Sediment" ~ "Lake\nSediment",
      Environment == "Shoreline Sediment" ~ "Shoreline\nSediment",
      Environment == "Stream Sediment" ~ "Stream\nSediment"),
    Environment = factor(Environment, levels = c("Starting\nConc.", "Compost", "Hydrocarbon\nSeep", "Lake\nSediment", "Shoreline\nSediment", "Stream\nSediment")))
head(ta_hplc2)
#View(ta_hplc2)

ggplot(ta_hplc2, aes(x = Environment, y = Terephthalamide_g_L))+
  #facet_grid(cols = vars(Environment), scales = "free_x", space = "free_x")+
  geom_col(aes(fill = Substrate), color = "black", size = 1)+
  geom_errorbar(aes(ymax = maxTA, ymin = minTA), width = 0.5, color = "black", size = 1)+
  scale_fill_manual(values = colors)+
  scale_color_manual(values = colors)+
  xlab(NULL)+
  ylab("Terephthalamide (g/L)")+ 
  theme_linedraw()+
  theme(legend.text = element_text(size = 22),
        legend.title = element_text(face = "bold", size = 22),
        legend.position = "bottom",
        axis.text.y = element_text(size = 22),                  
        axis.text.x = element_text(size = 22),           
        axis.title.x = element_text(size = 30), 
        axis.title.y = element_text(size = 30),                     
        strip.text = element_text(face = "bold", size = 30, color = "black"),
        strip.background = element_rect(fill="white"))+
  guides(fill = guide_legend(nrow = 1))
