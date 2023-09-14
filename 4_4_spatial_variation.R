### Spatial variation 

# Load required libraries
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)

# Define subregions of Europe
subregions <- c("Northern Europe", "Southern Europe", "Eastern Europe", "Western Europe")

# Get world data with a resolution of 1:50m
world_data <- ne_countries(scale = "medium", returnclass = "sf")
#world_data <- world_data[world_data$region_wb == "Europe & Central Asia",]

# Filter data to include only European countries
#europe_data <- world_data[world_data$subregion %in% subregions,]

# Filter to only countries in our data
countries_ecdc <- read.csv("data/countries.csv")[,-1]
# GB NOT UK in iso_a2
countries_ecdc[which(countries_ecdc == "UK")] <- "GB"
# GR NOT EL for greece in iso_a2
countries_ecdc[which(countries_ecdc == "EL")] <- "GR"
europe_data <- world_data[world_data$iso_a2 %in% countries_ecdc,]

# Create a ggplot map
europe_map <- ggplot() +
  geom_sf(data = world_data , fill = "gray96") +
  theme_minimal() +
  geom_sf(data = europe_data , aes(fill = name_long)) +
 # scale_fill_manual(values = c("Northern Europe" = "red",
   #                            "Southern Europe" = "blue",
  #                             "Eastern Europe" = "green",
  #                             "Western Europe" = "orange")) +
 labs(title = "Countries in analysis",
       fill = "Country") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") + 
  coord_sf(xlim = c(-25, 37), ylim = c(33, 73), expand = FALSE)

# Display the map
print(europe_map)
ggsave("plots/countries.tiff")

## Colour by subregion
# Create a ggplot map
europe_map <- ggplot() +
  geom_sf(data = world_data , fill = "gray96") +
  theme_minimal() +
  geom_sf(data = europe_data , aes(fill = subregion)) +
  scale_fill_manual(values = c("Northern Europe" = "red",
                             "Southern Europe" = "blue",
                              "Eastern Europe" = "green",
                              "Western Europe" = "orange")) +
  labs(title = "Countries in analysis",
       fill = "Subregion") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") + 
  coord_sf(xlim = c(-25, 37), ylim = c(33, 73), expand = FALSE)

# Display the map
print(europe_map)
ggsave("plots/countries_subregion.tiff")

### 
europe_data[,c("subregion","iso_a2")]

### Mean age
## From https://worldpopulationreview.com/country-rankings/median-age
mean_age <- read.csv("data/mean_age_by_country.csv") %>% filter(cca2 %in% countries_ecdc) %>% rename(iso_a2 = cca2)
mean_age_grp <- mean_age %>% group_by(subregion) %>% summarise(av_age = median(MedianAge, na.rm = TRUE), av_age_female = median(MedianAgeFemale, na.rm = TRUE), av_age_male = median(MedianAgeMale, na.rm = TRUE))


new_eu_data_raw <- left_join(europe_data, mean_age)
new_eu_data <- left_join(europe_data, mean_age_grp)

# Create a ggplot map
europe_map <- ggplot() +
  geom_sf(data = world_data , fill = "gray96") +
  theme_minimal() +
  geom_sf(data = new_eu_data , aes(fill = av_age)) +
 scale_fill_continuous(lim = c(40,50), type = "viridis") + 
  labs(title = "Median age", fill = "Median age") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") + 
  coord_sf(xlim = c(-25, 37), ylim = c(33, 73), expand = FALSE)

# Display the map
print(europe_map)
ggsave("plots/countries_subregion_meanage.tiff")

### By gender
# Create a ggplot map
europe_mapm <- ggplot() +
  geom_sf(data = world_data , fill = "gray96") +
  theme_minimal() +
  geom_sf(data = new_eu_data , aes(fill = av_age_male)) +
  scale_fill_continuous(lim = c(40,50), type = "viridis") + 
  labs(title = "Men", fill = "Median age") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") + 
  coord_sf(xlim = c(-25, 37), ylim = c(33, 73), expand = FALSE)

# Create a ggplot map
europe_mapf <- ggplot() +
  geom_sf(data = world_data , fill = "gray96") +
  theme_minimal() +
  geom_sf(data = new_eu_data , aes(fill = av_age_female)) +
  scale_fill_continuous(lim = c(40,50), type = "viridis") + 
  labs(title = "Women", fill = "Median age") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") + 
  coord_sf(xlim = c(-25, 37), ylim = c(33, 73), expand = FALSE)

# Display the map
europe_mapf + europe_mapm
ggsave("plots/countries_subregion_meanage_bygender.tiff")

### Explore a bit
summary(mean_age$MedianAge)
summary(mean_age$MedianAgeFemale)
summary(mean_age$MedianAgeMale)

### By gender
# Create a ggplot map
europe_mapm <- ggplot() +
  geom_sf(data = world_data , fill = "gray96") +
  theme_minimal() +
  geom_sf(data = new_eu_data , aes(fill = av_age_male)) +
  scale_fill_continuous(lim = c(40,50), type = "viridis") + 
  labs(title = "Men", fill = "Median age") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") + 
  coord_sf(xlim = c(-25, 37), ylim = c(33, 73), expand = FALSE)

# Create a ggplot map
europe_mapf <- ggplot() +
  geom_sf(data = world_data , fill = "gray96") +
  theme_minimal() +
  geom_sf(data = new_eu_data , aes(fill = av_age_female)) +
  scale_fill_continuous(lim = c(40,50), type = "viridis") + 
  labs(title = "Women", fill = "Median age") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") + 
  coord_sf(xlim = c(-25, 37), ylim = c(33, 73), expand = FALSE)

# Display the map
europe_mapf + europe_mapm
ggsave("plots/countries_subregion_meanage_bygender.tiff")

### By gender
# Create a ggplot map
europe_mapm <- ggplot() +
  geom_sf(data = world_data , fill = "gray96") +
  theme_minimal() +
  geom_sf(data = new_eu_data_raw , aes(fill = MedianAgeMale)) +
  scale_fill_continuous(lim = c(30,50), type = "viridis") + 
  labs(title = "Men", fill = "Median age") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") + 
  coord_sf(xlim = c(-25, 37), ylim = c(33, 73), expand = FALSE)

# Create a ggplot map
europe_mapf <- ggplot() +
  geom_sf(data = world_data , fill = "gray96") +
  theme_minimal() +
  geom_sf(data = new_eu_data_raw , aes(fill = MedianAgeFemale)) +
  scale_fill_continuous(lim = c(30,50), type = "viridis") + 
  labs(title = "Women", fill = "Median age") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") + 
  coord_sf(xlim = c(-25, 37), ylim = c(33, 73), expand = FALSE)
europe_mapf + europe_mapm
ggsave("plots/countries_medianage_bygender.tiff")
