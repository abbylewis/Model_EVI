# More vis outputs
library(tidyverse)
library(sf)
library(rjags)
library(ggmcmc)

biomass_points <- read_csv("Raw_data/NASA_blue_methane_plant_plot_summary.csv")

biomass_xy <- biomass_points %>% 
  filter(complete.cases(Latitude_dd, Longitude_dd)) %>% 
  group_by(Site, Subsite) %>% 
  # Avg lat lon
  summarise(lat = mean(Latitude_dd),
            lon = mean(Longitude_dd))

ll_crs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
aea_crs <- "+proj=aea +ellps=WGS84 +lat_1=29.5 +lat_2=45.5 +lon_0=-96 +x_0=0 +y_0=0"

biomass_sf <- st_as_sf(biomass_xy,
                       coords = c("lon", "lat"),
                       crs = ll_crs)

# Project
biomass_sf_aea <- st_transform(biomass_sf, 
                            aea_crs)

uvvr_tab <- data.frame(Site = c("LUMCON", "LUMCON", "LUMCON", 
                                "SERC", "SERC","SERC"),
                       Subsite = c("LUM1","LUM2","LUM3",
                                   "Mud","Hog","Kirk"),
                       UVVR = c(0.9576, 0.0307, 0.0147,
                                0, 0.0661, 0.1003)) %>% 
  mutate(fraction_veg = 1/(UVVR+1))

model_values_wide <- read_csv("Output/fits_beck_wide_analysis_ready.csv")


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

biomass_sf_aea_join <- biomass_sf_aea %>% 
  st_join(model_values_aea, join = st_nearest_feature)

model_values_filter <- model_values_wide %>% 
  filter(index_2 %in% biomass_sf_aea_join$index_2) %>% 
  left_join(biomass_sf_aea_join) %>% 
  select(-geometry) %>% 
  rename(Year = year) %>% 
  select(Site, Subsite, Year, eos:sos)

biomass_points_plot <- biomass_points %>% 
  left_join(uvvr_tab) %>% 
  mutate(Biomass_g_per_m2_se = Biomass_g_per_m2_se*fraction_veg,
         Biomass_g_per_m2 = Biomass_g_per_m2*fraction_veg) %>% 
  left_join(model_values_filter) %>% 
  mutate(doy = yday(ymd(paste(Year, Month, Day, sep="-"))),
    evi = (mn + (mx - mn) * ((1/(1+exp(-rsp*(doy-sos))))  + (1/(1+exp(rau*(doy-eos))))))) 
  
ggplot(biomass_points_plot, aes(x = evi, y = Biomass_g_per_m2)) +
  geom_point(aes(color = Site, 
                 shape = as.character(Year))
             ) + 
  geom_segment(aes(color = Site, 
                   y = Biomass_g_per_m2-Biomass_g_per_m2_se, 
                   yend = Biomass_g_per_m2+Biomass_g_per_m2_se)) +
  geom_smooth(method = "lm", color = "black") +
  theme_minimal() +
  xlab("Enhanced Vegetation Index") +
  ylab(expression("Biomass (g m"^"-2"~")")) +
  theme(legend.title = element_blank()) +
  scale_y_log10()

biomass_jaggified <- biomass_points_plot %>% 
  filter(Biomass_g_per_m2 > 0) %>% 
  mutate(Plot = as.numeric(factor(paste(Site, Subsite, Plot, Year, Month, Day))),
         Sampling_Event = as.numeric(factor(paste(Site, Subsite, Plot, Year, Month, Day))),
         Site = as.numeric(factor(Site)),
         ln_Biomass_g_per_m2 = log(Biomass_g_per_m2),
         tau_obs = 1/(Biomass_g_per_m2_se/Biomass_g_per_m2)^2) %>% 
  select(Site, Sampling_Event, Plot, evi, ln_Biomass_g_per_m2, tau_obs, Area_m2) %>% 
  arrange(Site, Sampling_Event, Plot)

biomass_subsite_site <- biomass_jaggified %>% 
  select(Site, Sampling_Event, evi) %>% 
  distinct_all()

jags_data <- list(N_sampling_events = max(biomass_jaggified$Sampling_Event),
                  Site_short_vect = biomass_subsite_site$Site,
                  evi = biomass_subsite_site$evi,
                  
                  N = nrow(biomass_jaggified),
                  Site_long_vect = biomass_jaggified$Site,
                  Sampling_Event_long_vect = biomass_jaggified$Sampling_Event,
                  Plot = biomass_jaggified$Plot,
                  
                  ln_Biomass_g_per_m2 = biomass_jaggified$ln_Biomass_g_per_m2,
                  tau_obs = biomass_jaggified$tau_obs
                  )

# Jags model 
jags_model <- "model{

  # Set priors

  # Uninformed priors for relationship between EVI and biomass
  #evi_slope_mean ~ dnorm(0, 0.01)
  #evi_intercept_mean ~ dnorm(0, 0.01)
  
  #evi_slope_tau ~ dgamma(0.001, 0.001)
  #evi_intercept_tau ~ dgamma(0.001, 0.001)
  
  #sub_pixel_alpha ~ dgamma(0.001,0.001)
  #sub_pixel_beta ~ dgamma(0.001,0.001)
  
  #model_alpha ~ dgamma(0.001,0.001)
  #model_beta ~ dgamma(0.001,0.001)
  
  tau_model ~ dgamma(0.001, 0.001)
  evi_slope ~ dnorm(0, 0.01)
  evi_intercept ~ dnorm(0, 0.01)
  
  # Uninformed for sub-pixel heterogeneity
  for (k in 1:2) {
  
    tau_sub_pixel[k] ~ dgamma(0.001, 0.001)
    #tau_model[k] ~ dgamma(0.001, 0.001)
    
    #evi_slope[k] ~ dnorm(0, 0.01)
    #evi_intercept[k] ~ dnorm(0, 0.01)
    
  }
  
  # For each site-subplot combo
  for (j in 1:N_sampling_events) {
  
    # Sub pixel biomass is a linear function of EVI
    modeled_bmass[j] <- evi[j] * evi_slope + evi_intercept
    
    z_bmass[j] ~ dnorm(modeled_bmass[j], tau_model)
    
  }
  
  # For each observation
  for (i in 1:N) {

    # A true plot biomass is a function of pixel biomass, sub-pixel random effect, and plot level precision
    z_plot_biomass[i] ~ dnorm(z_bmass[Sampling_Event_long_vect[i]], tau_sub_pixel[Site_long_vect[i]])
  
    # Observed biomass is a function of true biomass and observation precision
    ln_Biomass_g_per_m2[i] ~ dnorm(z_plot_biomass[i], tau_obs[i])
  
  }
  

}"

# Burn in 
j.model   <- jags.model(file = textConnection(jags_model),
                        data = jags_data,
                        n.chains = 4)

# Track slope, and intercept
var.out   <- coda.samples(model = j.model,
                          variable.names = c("evi_slope",
                                             "evi_intercept",
                                             "tau_model",
                                             "tau_sub_pixel"),
                          n.iter = 2000)

library(ggmcmc)

tidy_jags <- ggs(var.out)

ggs_traceplot(tidy_jags)

tidy_sum <- tidy_jags %>% 
  group_by(Parameter) %>% 
  filter(Iteration > 1000) %>% 
  summarise(mean = mean(value),
            se = sd(value),
            upper_CI = quantile(value, 0.975),
            lower_CI= quantile(value, 0.025))
(tidy_sum)
write_csv(tidy_sum, "Output/EVI_biomass_model_summary.csv")

create_nice_line <- tidy_jags %>% 
  filter(Iteration > 1000) %>% 
  separate(Parameter, into =c("Parameter", "Site"), sep = "\\[") %>% 
  mutate(Site = str_remove_all(Site, "\\]"),
         Site = recode(Site, "1"="LUMCON", "2"="SERC")) %>% 
  filter(Parameter %in% c("evi_intercept", "evi_slope")) %>% 
  spread(key = Parameter, value = value) %>% 
  mutate(index = 1:n())

join_these <- expand.grid(index = 1:nrow(create_nice_line), evi = seq(min(biomass_jaggified$evi), max(biomass_jaggified$evi), by = 0.01))

join_these <- join_these %>% left_join(create_nice_line)

bmass_pred <- join_these %>% 
  mutate(Biomass_g_per_m2 = exp(evi_intercept + evi * evi_slope)) %>% 
  group_by(evi) %>% 
  summarise(median_Biomass_g_per_m2 = median(Biomass_g_per_m2),
            upperCI = quantile(Biomass_g_per_m2, 0.975),
            lowerCI = quantile(Biomass_g_per_m2, 0.025))

ggplot(biomass_points_plot, aes(x = evi, y = Biomass_g_per_m2)) +
  geom_segment(aes(color = Site, 
                   y = Biomass_g_per_m2-Biomass_g_per_m2_se, 
                   yend = Biomass_g_per_m2+Biomass_g_per_m2_se)) +
  geom_line(data = bmass_pred, color = "black", aes(x = evi, y = median_Biomass_g_per_m2)) +
  geom_point(aes(color = Site, 
                 shape = as.character(Year))
  ) + 
  geom_ribbon(data = bmass_pred, aes(x = evi, ymin = lowerCI, ymax=upperCI), color = "black", fill = "grey", alpha = 0.3,
              inherit.aes = FALSE) +
  theme_minimal() +
  # facet_wrap(.~Site) +
  xlab("Predicted Enhanced Vegetation Index") +
  ylab(expression("Biomass (g m"^"-2"~")")) +
  theme(legend.title = element_blank())

ggsave("Biomass_prediction_figure.jpg", width = 5.5, height = 4)
            