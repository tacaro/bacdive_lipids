# bacdive_project
This tool allows us to pull bacterial strain information from the [BacDive Database](https://bacdive.dsmz.de/), extract their fatty acid profiles, and select fatty acid profiles based on bacterial taxa we care about.

I demonstrate this function by using the 20 most abundant bacterial taxa identified in an [aircraft sampling campaign](https://github.com/tacaro/aeroDADA2).


## Repository Contents
- README.md: This file
- bacdive_miner.R: R script that parses fatty acid profiles for certain bacterial taxa
- top_20_taxa.csv: A comma-separated file that lists the most abundant bacterial taxa from the GLDS-170 dataset. This file is structured as a taxa table that results from the [DADA2 pipeline](https://benjjneb.github.io/dada2/).
- export_bacdive_iso_table.csv: Downloaded from BacDive, a list of isolation sources.
- bacdive_fatty_acids_all_data.rds: A large cache. Contains every BacDive entry with a valid fatty acid profile.
- .PNG images: Results from bacdive_miner.R

## Method
Goal: Given a list of bacterial taxa, access the BacDive database and retrieve fatty acid profiles.
0. Set working directory to the location of your clone (line 7).

1. **Create cache of BacDive Fatty Acid Data**
**Option 1: Using the pre-packaged cache**
The easiest option is to use the cache that I have provided. This contains every BacDive entry containing a valid Fatty Acid profile. To import this cache, execute on line 46:
`fa_bugs <- read_rds(file.path("bacdive_fatty_acids_all_data.rds"))`

**Option 2: Creating your own BacDive Cache**
**NOTE**: *To request data from BacDive using their RESTful API, an account needs to be created and [BacDiveR](https://github.com/TIBHannover/BacDiveR) library must be installed. To avoid this, follow option 1. If following option 1, do not worry that library(BacDiveR) returns an error.*
Begin by reading in a single taxon, as a sanity check (lines 16-18).
Then, use `bd_retrieve_by_search()` followed by `write_rds()` to write the entries as a .rds file.

2. Now that we have our BacDive fatty acid database, we need to extract some important information from it. First, let's isolate the strain names in a new tibble called `strains`:
```strains <-
  tibble(
    ID = names(fa_bugs),
    data = map(fa_bugs, ~.x$taxonomy_name$strains)
    ) %>% unnest(data)
```
And the same with the fatty acid data:
```fas <-
  tibble(
  ID = names(fa_bugs),
  data = map(fa_bugs, ~.x$morphology_physiology$FA_profile)
  ) %>% unnest(data)
```
*Optional:* We can import the bacdive isolation sources, if we are curious about correlating fatty acid data to where the bacteria was isolated. This is outside the scope of this project. See lines 62 - 84.

3. Import the taxa-of-interest. In our case, the top 20 taxa identified in the aforementioned aerosol sampling study. `tax_table <- read.csv("top_20_tax.csv")` Note that the amplicon sequencing done in this experiment only resolves bacterial taxonomy to the genus level.

4. Identify which of our taxa-of-interest are included in the BacDive database:
```
matching_strains <- strains %>%
  filter(genus %in% tax_table$Genus)

matching_genus <- matching_strains %>%
  select(ID, genus)
```

5. Next, we need to create a tibble containing the fatty acid profiles of the taxa-of-interest (lines 99-106). We begin with the large fatty acid database we defined as `fas`
First, we use `filter()` using BacDive IDs such that we include only taxa-of-interest.
Second, we use `mutate()` to remove the FA_STDV column as this column is empty.
Third, grouping by BacDive ID and the fatty acid, we average the elution times and relative abundance of the fatty acids for each bacteria across different MS runs.
Lastly, we use a `left_join()` to add genus names to the tibble.

6. We then collapse multiple species of the same genus into a single fatty acid profile by averaging the relative fatty acid compositions. We do this with the `summarize()` function.

7. Finally, we use a for loop to plot each fatty acid profile for our taxa-of-interest. This results in the plots contained in this repository. For example:
![Acinetobacter Fatty Acid Profile](/Acinetobacter\ fatty_acid_profile.png)
