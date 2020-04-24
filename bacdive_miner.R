# R setup

library(tidyverse)
library(BacDiveR)

# Read bacdive

# Start by pulling all entries for a species
taxon_1 <- "Bacillus halodurans"
bac_data <- bd_retrieve_taxon(name = taxon_1)


# Let's try a specific id
bug <- BacDiveR:::bd_retrieve_data(searchTerm="133323", searchType = "bacdive_id")
# Let's try pulling the fatty acid profile from an entry
bug$`133323`$morphology_physiology$FA_profile

# Now we can try this on the large scale
# This takes 20 - 40 minutes
# We're using an advanced search query where:
#   BacDive entries that have fatty acids with number 1 in their FA profiles
#   have FA profiles that exist.

fa_bugs <- bd_retrieve_by_search(
  queryURL=
    paste("https://bacdive.dsmz.de/advsearch?site=advsearch",
          "searchparams%5B1149%5D%5Bcontenttype%5D=text",
          "searchparams%5B1149%5D%5Btypecontent%5D=contains",
          "searchparams%5B1149%5D%5Bsearchterm%5D=1",
          "advsearch=search",
          sep="&")
)

# cache it to avoid redundant downloading
fa_bugs %>% write_rds(file.path("bacdive_fatty_acids_all_data.rds"))

# Load from Cache
fa_bugs <- read_rds(file.path("bacdive_fatty_acids_all_data.rds"))

# We've got our data structure, now let's parse it for what we care about

strains <- 
  tibble(
    ID = names(fa_bugs),
    data = map(fa_bugs, ~.x$taxonomy_name$strains)
    ) %>% unnest(data)

fas <- 
  tibble(
  ID = names(fa_bugs),
  data = map(fa_bugs, ~.x$morphology_physiology$FA_profile) 
  ) %>% unnest(data)

# download the isolation sources from https://bacdive.dsmz.de/isolation-sources
source <- read_csv(
  file.path("export_bacdive_iso_table.csv"),
  col_types = cols(
    ID = col_character(),
    Species = col_character(),
    `Culture collection number` = col_character(),
    `Isolation source` = col_character(),
    Country = col_character(),
    Continent = col_character(),
    `Category 1` = col_character(),
    `Category 2` = col_character(),
    `Category 3` = col_character()
  )
)

unique_source <- select(source, ID, `Category 1`) %>% unique()

# Count the number of isolates from each source
all_count <- unique_source %>% count(`Category 1`)
environmental_count <- 
  source %>% filter(`Category 1` == "#Environmental") %>% 
  select(`ID`, `Category 2`) %>% unique() %>% count(`Category 2`)


