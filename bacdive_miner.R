# R setup

library(tidyverse) # Main library used in this project
library(magrittr) # Useful for pipe %>% and T %T>%
library(BacDiveR) # Allows us to access bacdive RESTful API and download entries
library(ggrepel) # Useful for plotting labels on ggplots

# Read bacdive

# Start by pulling all entries for a species
taxon_1 <- "Bacillus halodurans"
bac_data <- bd_retrieve_taxon(name = taxon_1)


# Let's try a specific id
bug <- BacDiveR:::bd_retrieve_data(searchTerm="133323", searchType = "bacdive_id")
# Let's try pulling the fatty acid profile from an entry
bug$`133323`$morphology_physiology$FA_profile

# Assuming this works:
# Now we can try this on the large scale
# This takes 20 - 40 minutes
# We're using an advanced search query where:
#   BacDive entries that have fatty acids with number 1 in their FA profiles
#   have FA profiles that exist.


##########
# ONLY EXECUTE THIS CODE IF YOU NEED A FRESH COPY OF THE FATTY ACID DATA
fa_bugs <- bd_retrieve_by_search(
  queryURL=
    paste("https://bacdive.dsmz.de/advsearch?site=advsearch",
          "searchparams%5B1149%5D%5Bcontenttype%5D=text",
          "searchparams%5B1149%5D%5Btypecontent%5D=contains",
          "searchparams%5B1149%5D%5Bsearchterm%5D=1",
          "advsearch=search",
          sep="&")
)
# cache the rds to avoid redundant downloading
fa_bugs %>% write_rds(file.path("bacdive_fatty_acids_all_data.rds"))

#########


# Load from Cache
fa_bugs <- read_rds(file.path("/Users/Tristan/bacdive_mineR/bacdive_fatty_acids_all_data.rds"))

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

# We can download the isolation sources from https://bacdive.dsmz.de/isolation-sources
source <- read_csv(
  file.path("/Users/Tristan/bacdive_mineR/export_bacdive_iso_table.csv"),
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



# Read in top 20 taxa from aerobiology study
tax_table <- read.csv("top_20_tax.csv")

# Create a tibble of the strains that are in the tax_table
matching_strains <- strains %>%
  filter(genus %in% tax_table$Genus)

matching_genus <- matching_strains %>% 
  select(ID, genus)

# Create a tibble of the fatty acids of the strains in the tax_table
matching_fas <- fas %>%
  filter(ID %in% matching_strains$ID) %>% # check if the ID is in the matching strains tibble
  mutate(FA_STDV = NULL) %>% # get rid of the FA_STDV column
  group_by(ID, FA_name) %>%
  mutate(FA_mn_ECL = mean(FA_ECL), FA_mn_val = mean(FA_percent)) %>% # average the different FA measurements
  left_join(matching_genus, by = "ID") # add the genus strings to the tibble
matching_fas$ID = as.numeric(as.character(matching_fas$ID))
matching_fas %<>% arrange(genus) #arrange in ascending order
# This tibble is useful for plotting individual fatty acid profiles.
# Let's take Alishewanella fetalis (ID = 440) as an example:
test_440_fa <- matching_fas %>%
  filter(ID == "440")

test_440_str <- matching_strains %>% 
  filter(ID == "440")

test_440_plt <- ggplot(test_440_fa, aes(x = FA_mn_ECL, y = FA_mn_val)) +
  geom_bar(stat="identity") +
  geom_text_repel(label=test_440_fa$FA_name, nudge_y = 0.5, size = 2, color = "red") +
  ggtitle(paste(test_440_str$species, "Fatty Acid Profile")) +
  ylab("Percent composition") +
  xlab("Elution Time (minutes)") +
  theme_bw()
test_440_plt # It works!
ggsave(plot=test_440_plt, filename=as.character(paste(test_440_str$species, "Fatty Acid Profile", ".png")))




# In order to see FA profiles at the genus level, we need to compress
# the FA profiles of bacteria within the same genus:
require(dplyr) # just to be safe
matching_fas_genus <- matching_fas %>%
  group_by(genus, FA_name) %>% # we want FA data to be described for each genus
  drop_na() %>%
  summarize(FA_std_val = sd(FA_mn_val),
            FA_std_ECL = sd(FA_mn_ECL),
            FA_mn_val = mean(FA_mn_val),
            FA_mn_ECL = mean(FA_mn_ECL) # take sd of % composition
            ) %>%
  drop_na()
  
for (gns in unique(matching_fas_genus$genus)) {
  print(gns)
  df <- matching_fas_genus %>% filter(genus == gns)
  plt <- ggplot(df, aes(x = FA_mn_ECL, y = FA_mn_val)) +
          geom_bar(stat = "identity") +
          geom_text_repel(label=df$FA_name, nudge_y = 0.5, size = 2, color = "red", segment.size = 0.1) +
          ggtitle(paste(gns, "fatty acid profile")) +
          ylab("Percent composition") +
          xlab("Elution Time (minutes)") +
          theme_bw()
  print(plt)
  ggsave(plot = plt, filename=as.character(paste(gns, "fatty_acid_profile.png")))
}
# Done!