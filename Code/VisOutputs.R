# Try something
library(tidyverse)
library(sf)

# First let's see if queries for SERC

# Load up point data
pw_points <- read_csv("Raw_data/All_Data_Porewater_Biomass_Richness_XYZPrecision_Zstar.csv")

# Filter
serc_pw <- pw_points %>% 
  filter(Site %in% c("SERC", "LUMC")) %>% 
  group_by(Site, Subsite) %>% 
  # Avg lat lon
  summarise(lat = mean(Latitude_dd),
            lon = mean(Longitude_dd))

ll_crs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
aea_crs <- "+proj=aea +ellps=WGS84 +lat_1=29.5 +lat_2=45.5 +lon_0=-96 +x_0=0 +y_0=0"

serc_pw_sf <- st_as_sf(serc_pw,
         coords = c("lon", "lat"),
         crs = ll_crs)

st_write(serc_pw_sf, "SERC_and_lumcon_veg_sites.shp")

# Project
serc_pw_aea <- st_transform(serc_pw_sf, 
                            aea_crs)

# Load up grid data
pheno_points1 <- read_csv("Raw_data/smithsonianenvironmentalresearchcenter_evi_series.csv")
pheno_points2 <- read_csv("Raw_data/louisianauniversitiesmarineconsortium_evi_series.csv")

pheno_points <- bind_rows(pheno_points1,
                          pheno_points2)

# Get unique pixels
pheno_xy <- pheno_points %>% 
  dplyr::select(latitude, longitude) %>% 
  distinct_all() %>% 
  mutate(index = 1:n())

pheno_points <- pheno_points %>% 
  left_join(pheno_xy)

# Project
pheno_xy_sf <- st_as_sf(pheno_xy,
           coords = c("longitude", "latitude"),
           crs = ll_crs)

pheno_xy_aea <- st_transform(pheno_xy_sf, 
                             aea_crs)

# Spatial join pixels
serc_pw_aea_join <- serc_pw_aea %>% 
  st_join(pheno_xy_aea, join = st_nearest_feature)

serc_pw_aea_join_ll <- st_transform(serc_pw_aea_join, ll_crs) 

# Filter to pixels
pheno_points_subset <- pheno_points %>% 
  filter(index %in% serc_pw_aea_join$index) %>% 
  left_join(serc_pw_aea_join_ll)



# Graph
ggplot(pheno_points_subset, aes(x = img_doy, y = evi)) +
  geom_point() +
  facet_grid(Site~Subsite~img_year)


# Check on model output
pheno_files_list <- c(list.files(path = "Output/", 
                                 pattern = "_PhenoBeckLong", recursive = TRUE, full.names = TRUE))

for (i in 1:length(pheno_files_list)) {
  if (i == 1) {
    model_values <- read_csv(pheno_files_list[i])
  } else {
    model_values <- bind_rows(model_values, read_csv(pheno_files_list[i]))
  }
}

model_values_wide <- model_values %>% 
  rename(year = img_year) %>% 
  dplyr::select(site_id, latitude, longitude, year, params, value) %>% 
  spread(key = params, value = value)

model_values_wide_se <- model_values %>% 
  rename(year = img_year) %>% 
  dplyr::select(site_id, latitude, longitude, year, params, se) %>%
  mutate(params = paste0(params, "_se")) %>% 
  spread(key = params, value = se)
  
model_values_wide <- model_values_wide %>% 
  left_join(model_values_wide_se)

model_values_mat1 <- model_values_wide %>% filter(site_id == "SERC")

pca_analysis_SERC <- prcomp(model_values_mat1 %>%  select(eos:sos))
plot(pca_analysis_SERC)
(pca_analysis_SERC)
biplot(pca_analysis_SERC)


normalize <- function(x, ...) {
  return(as.integer(265 * (x - min(x, ...)) /(max(x, ...) - min(x, ...))))
}


model_values_mat1[,c("PC1", "PC2", "PC3")] <- pca_analysis_SERC$x[,1:3]

model_values_mat1 <- model_values_mat1 %>% 
  mutate(R = normalize(PC1),
         G = normalize(PC2), 
         B = normalize(PC3))


model_values_mat2 <- model_values_wide %>% filter(site_id == "LUMCON")

pca_analysis_LUMCON <- prcomp(model_values_mat2 %>%  select(eos:sos))
plot(pca_analysis_LUMCON)
(pca_analysis_LUMCON)
biplot(pca_analysis_LUMCON)

model_values_mat2[,c("PC1", "PC2", "PC3")] <- pca_analysis_LUMCON$x[,1:3]

model_values_mat2 <- model_values_mat2 %>% 
  mutate(R = normalize(PC1),
         G = normalize(PC2), 
         B = normalize(PC3))

ouptut_w_pcas <- model_values_mat1 %>% 
  bind_rows(model_values_mat2)

write_csv(ouptut_w_pcas, "Output/fits_beck_wide_analysis_ready.csv")

model_values_xy <- model_values_wide  %>% 
  dplyr::select(latitude, longitude) %>% 
  distinct_all() %>% 
  mutate(index_2 = 1:n())

model_values_wide <- model_values_wide %>% 
  left_join(model_values_xy)

# Project
model_values_xy_sf <- st_as_sf(model_values_xy,
                        coords = c("longitude", "latitude"),
                        crs = ll_crs)

model_values_aea <- st_transform(model_values_xy_sf, 
                             aea_crs)

serc_pw_aea_join2 <- serc_pw_aea %>% 
  st_join(model_values_aea, join = st_nearest_feature)

model_values_filter <- model_values_wide %>% 
  filter(index_2 %in% serc_pw_aea_join2$index_2) %>% 
  left_join(serc_pw_aea_join2)

# View(model_values_filter)

pred_evi <- expand_grid(model_values_filter, doy = 1:365) %>% 
  mutate(evi = (mn + (mx - mn) * ((1/(1+exp(-rsp*(doy-sos))))  + (1/(1+exp(rau*(doy-eos))))))) %>% 
  rename(img_year=year)

ggplot(pheno_points_subset, aes(x = img_doy, y = evi)) +
  geom_point() +
  geom_line(data = pred_evi, color = "red", aes(x = doy), lty = 2) +
  facet_grid(Site~Subsite~img_year) +
  theme_minimal() +
  xlab("Day of Year") +
  ylab("Enhanced Vegetation Index")

ggsave("EVI_timer_series_picture.jpg", height = 6.5, width = 6.5)
