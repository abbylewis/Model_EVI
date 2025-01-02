# Jim's try at a function

# Internal functions
library(tidyverse)
library(purrr)
library(furrr)
source("Code/FitDoubleLogBeck.R")
options(future.rng.onMisuse = "ignore")

# Input data
evi_serc <- read_csv("Raw_data/smithsonianenvironmentalresearchcenter_evi_series.csv") %>% 
  mutate(site_id = "SERC") %>% 
  select(-img_month, -img_day) %>% 
  group_by(site_id, img_year, latitude, longitude) %>% 
  mutate(n = n()) %>% 
  filter(n > 10) %>% 
  ungroup() %>% 
  arrange(site_id, img_year, latitude, longitude, img_doy)

evi_lumcon <- read_csv("Raw_data/louisianauniversitiesmarineconsortium_evi_series.csv") %>% 
  mutate(site_id = "LUMCON") %>% 
  select(-img_month, -img_day) %>% 
  group_by(site_id, img_year, latitude, longitude) %>% 
  mutate(n = n()) %>% 
  filter(n > 10) %>% 
  ungroup() %>% 
  arrange(site_id, img_year, latitude, longitude, img_doy)

all_evi <- evi_serc %>% 
  bind_rows(evi_lumcon)

# Iterate by site and by year
iterator_df <- all_evi %>%
  select(site_id, img_year) %>% 
  distinct_all()

FitDoubleLogBeck_tidy <- function(id_i=1, df) {
  df <- df[id_i, ]
  t <- df$data[[1]]$img_doy
  evi <- df$data[[1]]$evi
    
  output_df <- data.frame( 
    params = c("mn",  "mx",  "sos", "rsp", "eos", "rau"),
                           value = rep(NA, 6),
                           se = rep(NA, 6))
   
  try({
    output_pheno <- FitDoubleLogBeck(x=evi, t=t, hessian=T, ninit = 100)
    output_df <- data.frame(params = names(output_pheno$params),
                            value = output_pheno$params,
                            se = output_pheno$stdError)
  })
  
  output_df <- expand_grid(df %>% select(-data), output_df)
  
  return(output_df)
}

for (i in 1:nrow(iterator_df)) {

  subset_df <- all_evi %>% 
    filter(img_year == iterator_df$img_year[i],
           site_id == iterator_df$site_id[i])
  
  out_csv_name <- paste0("Output/", iterator_df$site_id[i], "_",  iterator_df$img_year[i], "_PhenoBeckLong.csv")
  
  print(paste0("Running ", iterator_df$site_id[i], " ", iterator_df$img_year[i], "..."))
  
  if (! file.exists(out_csv_name)) {
    temp_df <- subset_df %>% 
      group_by(site_id, img_year, latitude, longitude) %>% 
      nest()
    
    plan(multisession,  workers = parallel::detectCores()-1)  
    pheno_out <- future_map(1:nrow(temp_df), ~FitDoubleLogBeck_tidy(id_i = .x, df=temp_df), .progress = T) %>% 
      bind_rows()
    
    write_csv(pheno_out, out_csv_name)
    
  }
  
}
 
