
FAOSTAT_pest_data <- read.csv("data/FAOSTAT_pest_data_en_3-4-2025.csv") # From https://www.fao.org/faostat/en/#data/RP
FAOSTAT_area_data <- read.csv("data/FAOSTAT_area_data_en_3-5-2025.csv") # From https://www.fao.org/faostat/en/#data/RL

landuse <- FAOSTAT_area_data %>% 
  filter(Element %in% c("Area", "Value of agricultural production (Int. $) per Area"), 
         Item %in% c("Land area", "Cropland", "Permanent crops", "Temporary crops", "Agricultural land")) %>%
  filter(!(Element == "Area" & Item == "Agricultural land")) %>%
  mutate(Item = ifelse(Element == "Value of agricultural production (Int. $) per Area", "agri_value", Item)) %>%
  select(Area, Year, Item, Value) %>% 
  pivot_wider(id_cols = c(Area, Year), values_from = Value, names_from = Item) %>%
  mutate(active_crops = `Temporary crops` + `Permanent crops`)

FAOSTAT_added_data <- FAOSTAT_pest_data %>%
  select(Area, Element, Item, Year, Value) %>%
  filter(!is.na(Value), 
         Year == 2019) %>%
  pivot_wider(id_cols = c(Area, Item, Year), names_from = Element, values_from = Value) %>%
  left_join(., landuse, by = c("Area", "Year")) %>%
  mutate(Country = case_when(Area == "Belgium" ~ "BE",
                             Area == "France" ~ "FR",
                             Area == "Germany" ~ "DE",
                             Area == "Netherlands (Kingdom of the)" ~ "NL",
                             Area == "Portugal" ~ "PT",
                             Area == "Romania" ~ "RO",
                             Area == "Switzerland" ~ "CH",
                             Area == "United Kingdom of Great Britain and Northern Ireland" ~ "UK"
  )) %>%
  mutate(Country = factor(Country, levels = c("PT", "FR", "UK", "BE", "NL", "CH", "DE", "RO"))) %>%
  mutate(`Use per area of cropland` = ifelse(is.na(`Use per area of cropland`), `Agricultural Use` / Cropland, `Use per area of cropland`),
         `Use per value of agricultural production` = ifelse(is.na(`Use per value of agricultural production`), `Agricultural Use` / agri_value, `Use per value of agricultural production`),
         use_per_active_cropland = `Agricultural Use` / active_crops,
         use_per_land_area = `Agricultural Use` / `Land area`) %>%
  left_join(., cropland_fraction, by = "Country") %>%
  mutate(Cropland_in_2km_radius = (pi*2000^2 / 10000) * cropland_fraction_2k_radius) %>%
  mutate(est_use_in_2k_radius = `Use per area of cropland` * Cropland_in_2km_radius) %>%
  select(Country, Item, Year, `Agricultural Use`, `Use per area of cropland`, 
         `Use per capita`, `Use per value of agricultural production`, 
         use_per_active_cropland, use_per_land_area, est_use_in_2k_radius, cropland_fraction_2k_radius, Cropland_in_2km_radius)
