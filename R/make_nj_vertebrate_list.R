library(dplyr)
source("R/align_taxonomy.R")

# list of NJ vertebrates compiled from:
# https://dep.nj.gov/njfw/fishing/freshwater/freshwater-fish-of-new-jersey/
# https://dep.nj.gov/njfw/education/online-field-guide-for-reptiles-and-amphibians/
# https://www.nj.gov/dep/fgw/chkbirds.htm
# https://www.nj.gov/dep/fgw/chkmamls.htm

nj <- read.csv("data/nj_vertebrates.csv") %>%
  select(-status)

# map_of_life_link	threat_status	expected	recorded	dataset_title	dataset_ids
nj_tax <- align_taxonomy_to_others(nj)

write.csv(nj_tax, "data/nj_vertebrates_taxonomy.csv", row.names = F)
