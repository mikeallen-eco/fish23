#' Get NCBI taxonomy from MOL species
#'
#' This function takes MOL species and tries to match them to NCBI taxonomy and also generates a MOL to GBIF key
#'
#' @param mol_df input df of species names from mol, where species column is `scientific_name`
#' @param mol_synonym_df df of synonyms from mol, set = NULL if none
#' @param save_gbif_key_path path to save mol to gbif key, include .csv. example: 'data/test.csv'
#' @param  verbose logical. logs updates in console
#'
#' @return a dataframe of mol taxonomy to matched NCBI taxonomy
#'
#' @examples align_mol_taxonomy_to_others(mol_df = mol_mammals, mol_synonym_df = NULL, save_gbif_key_path = "data/gbif/mol_gbif_mammals.csv")
#'
#' @export
align_taxonomy_to_others =
  function(mol_df,
           mol_synonym_df = NULL, 
           save_gbif_key_path = NULL, 
           verbose = TRUE)
  {
    library(dplyr)
    library(taxize)
    library(rgbif)
    library(stringr)
    library(glue)
    
    
    # set sleep timers
    short_sleep = 2
    long_sleep = 10
    
    
    # NCBI ----
    # get NCBI names for input taxa where available
    if(verbose == TRUE) message("# Starting Taxize match to NCBI #")
    
    results_ncbi_taxize = list()
    for(i in 1:nrow(mol_df)){
      
      if(verbose == TRUE) message(glue("NCBI processing {i} of {nrow(mol_df)} : {mol_df$scientific_name[i]}"))
      
      results_ncbi_taxize[[i]] <-
        tryCatch(
          {taxize::tax_name(sci = mol_df$scientific_name[i],
                            get = c("species", "genus", "family", "order", "class"),
                            db = "ncbi",
                            ask = FALSE, # NA returned for multiple
                            messages = FALSE) %>% 
              tibble() %>% 
              mutate(mol_name = mol_df$scientific_name[i], .before = 'db') %>% 
              select(-query)
          },
          error = function(e) {
            # Error handling code
            tryCatch(
              {
                # wait and try again
                if(verbose == TRUE) message(glue("Name {mol_df$scientific_name[i]} didnt work the first try, waiting a little and trying again."))
                Sys.sleep(long_sleep)
                taxize::tax_name(sci = mol_df$scientific_name[i],
                                 get = c("species", "genus", "family", "order", "class"),
                                 db = "ncbi",
                                 ask = FALSE, # NA returned for multiple
                                 messages = FALSE) %>% 
                  tibble() %>% 
                  mutate(mol_name = mol_df$scientific_name[i], .before = 'db') %>% 
                  select(-query)
              },
              error =  function(e2) {
                if(verbose == TRUE) message(glue("Name {mol_df$scientific_name[i]} didnt work the second try, not trying again, but noting this in 'notes'column"))
                tibble(notes = "error when looking up", mol_name = mol_df$scientific_name[i])
              } # end error second try
            ) # end second try
          }, # end error first try
          finally = {
            # wait for n seconds
            Sys.sleep(short_sleep)
          }
        ) # end first try
    } # end for loop
    
    
    # to one df
    results_ncbi1 = bind_rows(results_ncbi_taxize) %>% 
      mutate(match_ncbi = if_else(!is.na(species), 1, 0)) %>% 
      mutate(match_type = if_else(match_ncbi == 1, "direct MOL match", NA_character_))
    
    
    
    # GBIF ----
    # check which names are in GBIF (this also finds close matches that are spelling errors)
    if(verbose == TRUE) message("# Starting match to GBIF #")
    
    results_gbif = list()
    for(i in 1:nrow(mol_df)){
      if(verbose == TRUE) message(glue("GBIF processing {i} of {nrow(mol_df)} : {mol_df$scientific_name[i]}"))
      # need to do this next step in 2 parts as rgbif doesn't return order for reptiles or class for fish
      results_gbif_step1 <- rgbif::name_backbone(name = mol_df$scientific_name[i])
      if(is.null(results_gbif_step1$class)) {results_gbif_step1$class <- NA}
      if(is.null(results_gbif_step1$order)) {results_gbif_step1$order <- NA}
      
      results_gbif[[i]] <- results_gbif_step1 %>%
        mutate(mol_name = mol_df$scientific_name[i])
    } # end for loop
    
    gbif_df <- bind_rows(results_gbif) %>% 
      select(mol_name, gbif_rank = rank, matchType, kingdom, 
             phylum, class, order, family, genus, species)
    
    
    # NCBI part 2 ----
    # test GBIF synonyms for MOL names with still no valid NCBI name
    
    # which MOL spp that aren't yet matched & gbif match has diff spp name?
    mol_gbif_diff = results_ncbi1 %>% 
      filter(match_ncbi == 0) %>% 
      select(mol_name) %>% 
      left_join(gbif_df %>% select(mol_name, species), by = "mol_name") %>% 
      filter(mol_name != species)
    
    if(nrow(mol_gbif_diff) > 0){
      
      if(verbose == TRUE) message("# Starting match to NCBI Part 2 #")
      results_ncbi_taxize2 = list()
      for(i in 1:nrow(mol_gbif_diff)){
        if(verbose == TRUE) message(glue("NCBI 2nd processing {i} of {nrow(mol_gbif_diff)} : {mol_gbif_diff$species[i]}"))
        
        results_ncbi_taxize2[[i]] <- 
          tryCatch(
            {taxize::tax_name(sci = mol_gbif_diff$species[i],
                              get = c("species", "genus", "family", "order", "class"),
                              db = "ncbi",
                              ask = FALSE, # NA returned for multiple
                              messages = FALSE) %>% 
                tibble() %>% 
                mutate(gbif_name = mol_gbif_diff$species[i], .before = 'db') %>% 
                select(-query)
            },
            error = function(e) {
              # Error handling code
              tryCatch(
                {
                  # wait and try again
                  if(verbose == TRUE) message(glue("Name {mol_gbif_diff$species[i]} didnt work the first try, waiting a little and trying again."))
                  Sys.sleep(long_sleep)
                  taxize::tax_name(sci = mol_gbif_diff$species[i],
                                   get = c("species", "genus", "family", "order", "class"),
                                   db = "ncbi",
                                   ask = FALSE, # NA returned for multiple
                                   messages = FALSE) %>% 
                    tibble() %>% 
                    mutate(gbif_name = mol_gbif_diff$species[i], .before = 'db') %>% 
                    select(-query)
                },
                error =  function(e2) {
                  if(verbose == TRUE) message(glue("Name {mol_gbif_diff$species[i]} didnt work the second try, not trying again, but noting this in 'notes'column"))
                  tibble(notes = "error when looking up", gbif_name = mol_gbif_diff$mol_name[i])
                } # end error second try
              ) # end second try
            }, # end error first try
            finally = {
              # wait for n seconds
              Sys.sleep(short_sleep)
            }
          ) # end first try
      } # end for loop
      
      
      
      # to one df
      results_ncbi2 = bind_rows(results_ncbi_taxize2) %>% 
        left_join(mol_gbif_diff, by = c("gbif_name" = "species")) %>% 
        mutate(match_ncbi = if_else(!is.na(species), 1, 0)) %>% 
        mutate(match_type = if_else(match_ncbi == 1, "direct GBIF match", NA_character_))
      
    }else{ # end if
      # make empty df with colnames of results_ncbi1
      results_ncbi2 = data.frame(matrix(ncol = 2, nrow = 0))
      colnames(results_ncbi2) <- c("mol_name", "gbif_name")
      results_ncbi2$mol_name = as.character(results_ncbi2$mol_name)
      results_ncbi2$gbif_name = as.character(results_ncbi2$gbif_name)
    } 
    
    ## combine with original
    results_ncbi3 = 
      results_ncbi1 %>% 
      filter(!mol_name %in% (results_ncbi2 %>% select(-gbif_name) %>% pull(mol_name))  ) %>% 
      bind_rows(results_ncbi2 %>% select(-gbif_name)) %>%
      left_join(gbif_df %>% select(mol_name, gbif_name = species, gbif_genus = genus,
                                   gbif_family = family, gbif_order = order, 
                                   gbif_class = class)) %>%
      mutate(class = case_when(is.na(class) ~ gbif_class,
                               TRUE ~ class),
             order = case_when(is.na(order) ~ gbif_order,
                               TRUE ~ order),
             family = case_when(is.na(family) ~ gbif_family,
                                TRUE ~ family),
             genus = case_when(is.na(genus) ~ gbif_genus,
                               TRUE ~ genus)) %>%
      select(mol_name, ncbi_name = species, gbif_name, genus, family, order, class,
             match_ncbi, match_type) %>%
      # this fixes an issue with an unexpected many-to-many relationship 
      # but probably better to fix it before the join above
      distinct() 
    
    # write.csv(notes, paste0("data/", results_ncbi3$class[1], "_align_tax_error_notes.csv"), row.names = F)
    
    return(results_ncbi3)
    
  }
# # Synonyms ----
# 
# # Get all synonyms of unmatched species
# # Check they haven't already been searched
# unmatched_ncbi1 = results_ncbi3 %>% 
#   filter(match_ncbi == 0) %>% 
#   pull(mol_name)
# 
# 
# ## GBIF Synonyms ----
# ## grab syn names of unmatched species
# gbif_syns_list = list()
# for(i in 1:length(unmatched_ncbi1)){
#   
#   if(verbose == TRUE) message(glue("GBIF getting synonyms {i} of {length(unmatched_ncbi1)} : {unmatched_ncbi1[i]}"))
#   
#   gbif_syns_list[[i]] <- 
#   tryCatch(
#     {rgbif::name_lookup(unmatched_ncbi1[i])$data %>% 
#         filter(taxonomicStatus == "SYNONYM" & species != unmatched_ncbi1[i]) %>%
#         distinct(species) %>%
#         mutate(mol_name = unmatched_ncbi1[i])
#     },
#     error = function(e) {
#       # Error handling code
#       tryCatch(
#         {
#           # wait and try again
#           if(verbose == TRUE) message(glue("Name {unmatched_ncbi1[i]} didnt work the first try, waiting a little and trying again."))
#           rgbif::name_lookup(unmatched_ncbi1[i])$data %>% 
#             filter(taxonomicStatus == "SYNONYM" & species != unmatched_ncbi1[i]) %>%
#             distinct(species) %>%
#             mutate(mol_name = unmatched_ncbi1[i])
#         },
#         error =  function(e2) {
#           if(verbose == TRUE) message(glue("Name {unmatched_ncbi1[i]} didnt work the second try, not trying again, but noting this in 'notes'column"))
#           tibble(notes = "error when looking up", mol_name = unmatched_ncbi1[i])
#         } # end error second try
#       ) # end second try
#     }, # end error first try
#     finally = {
#       # wait for n seconds
#       Sys.sleep(short_sleep)
#     }
#   ) # end first try
# }# end for
# 
# gbif_syns = bind_rows(gbif_syns_list) %>% 
#   rename(syn_species = species) %>% 
#   mutate(syn_source = "gbif")
# 
# 
# ## MOL Syns ----
# 
# ## if syns are provided
# if(!is.null(mol_synonym_df)){
#   
#   # get syns of unmatched spp
#   mol_syns = mol_synonym_df %>%
#     distinct(mol_syn, .keep_all = T) %>%
#     rename(syn_species = mol_syn) %>% 
#     mutate(syn_source = "mol") %>% 
#     # filter only unmatched ones
#     filter(mol_name != syn_species) %>% 
#     filter(mol_name %in%
#              (results_ncbi3 %>% filter(match_ncbi == 0) %>% pull(mol_name))
#            )
#   
#   # combine with gbif
#   all_syns = bind_rows(gbif_syns, mol_syns)
# 
# }else{
#   all_syns = gbif_syns
# }
# 
# 
# ## Search syns ----
# syn_list = list()
# for(i in 1:nrow(all_syns)){
#   
#   if(verbose == TRUE) message(glue("NCBI 3rd processing {i} of {nrow(all_syns)} : {all_syns$syn_species[i]}"))
#   
#   syn_list[[i]] <- 
#     tryCatch(
#       {taxize::tax_name(sci = all_syns$syn_species  [i],
#                         get = c("species", "genus", "family", "order", "class"),
#                         db = "ncbi",
#                         messages = FALSE) %>% 
#           tibble() %>% 
#           mutate(mol_name = all_syns$mol_name[i], .before = 'db') %>% 
#           select(-query)
#       },
#       error = function(e) {
#         # Error handling code
#         tryCatch(
#           {
#             # wait and try again
#             if(verbose == TRUE) message(glue("Name {all_syns$species[i]} didnt work the first try, waiting a little and trying again."))
#             Sys.sleep(long_sleep)
#             taxize::tax_name(sci = all_syns$syn_species[i],
#                              get = c("species", "genus", "family", "order", "class"),
#                              db = "ncbi",
#                              messages = FALSE) %>% 
#               tibble() %>% 
#               mutate(mol_name = all_syns$mol_name[i], .before = 'db') %>% 
#               select(-query)
#           },
#           error =  function(e2) {
#             if(verbose == TRUE) message(glue("Name {all_syns$syn_species[i]} didnt work the second try, not trying again, but noting this in 'notes'column"))
#             tibble(notes = "error when looking up", mol_name = all_syns$mol_name[i])
#           } # end error second try
#         ) # end second try
#       }, # end error first try
#       finally = {
#         # wait for n seconds
#         Sys.sleep(short_sleep)
#       }
#     ) # end first try    
# }
# 
# 
# syn_results = bind_rows(syn_list) %>% 
#   mutate(match_ncbi = if_else(!is.na(species), 1, 0)) %>% 
#   mutate(match_type = if_else(match_ncbi == 1, "syn match", NA_character_))
# 
# # combine
# results_ncbi4 = 
#   results_ncbi3 %>% 
#   filter(!mol_name %in% (syn_results %>% filter(match_ncbi ==1) %>% pull(mol_name))  ) %>% 
#   bind_rows(syn_results %>% filter(match_ncbi ==1))
#
# return(results_ncbi4)
#   
# }
