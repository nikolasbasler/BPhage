pick_ambiguous_taxa <- function(vcontact_output, taxlevel_to_pick) {
  # Pick the subfamily/genus that corresponds to the family/subfamily 
  # prediction. This is necessary because with the vcontact3 beta version used
  # here, sometimes several equally valid taxa are listed on one taxonomic 
  # level (here only for subfamily and genus). These are then mentioned in a 
  # double pipe-separated list (||, read "or of"), but the name of every is 
  # truncated, apparently to save characters. E.g.:
  # novel_subfamily_0_of_novel_family_1_of_novel_order_116_of_Caudoviricetes||novel_family_8_of_novel_order_86_of_Caudoviricetes
  # Unfortunately, the first item in this list does't necessarily correspond to
  # the higher taxonomic level, so in the example above, the family prediction
  # might be "novel_family_8_of..." So it's a bit tricky to pick and complete 
  # the subfamily/genus name (to "novel_subfamily_0_of_novel_family_8_of..." in 
  # this case).
  
  # So run this function twice, once with taxlevel_to_pick="Subfamily", once 
  # with taxlevel_to_pick="Genus"
  # This function works for the bphage dataset and the contigs for the
  # additional datasets, but there is no guarantee that it would universally
  # work. In fact, I know it won't work if the subfamily isn't 
  # "novel_subfamily_xx" but one with an actual name (of which there isn't any
  # case in this dataset).
  # Hopefully the vcontact3 devs will one day make it so that in case there are
  # several possible taxa on one level, the taxon that corresponds to the higher
  # level is mentioned first in this double pipe-separated list.
  
  if (taxlevel_to_pick == "Subfamily") {
    higher_tax <- "family..prediction."
    lower_tax <- "subfamily..prediction."
    start_string <- "novel_subfamily"
  }
  if (taxlevel_to_pick == "Genus") {
    higher_tax <- "Subfamily"
    lower_tax <- "genus..prediction."
    start_string <- "novel_genus"
  }
  
  # In the Genus pass of this function, prepend the subfamily to the fields in 
  # the genus..prediction. column where it's missing. E.g. in cases like this:
  # novel_genus_5_of_novel_subfamily_0_of_novel_family_4_of_novel_order_170_of_Caudoviricetes||novel_family_17_of_novel_order_172_of_Caudoviricetes
  if (taxlevel_to_pick == "Genus") {
    vcontact_output <- vcontact_output %>%
      mutate(genus..prediction. = str_replace_all(
        genus..prediction.,
        "\\|\\|novel_family",
        paste0("||", str_extract(Subfamily, "^([^_]+_){3}[^_]+"), "_novel_family")))
    }
  
  # Make a new column where each row contains a list with the doublepipe-
  # separated fields of the "subfamily..prediction." column (catch empty 
  # characters and set to NA). Then pick the item from the list that contains 
  # the substring in "family..prediction.".
  vcontact_output %>%
    mutate(new_column = str_split(.data[[lower_tax]], "\\|\\|")) %>% 
    mutate(
      new_column = map2_chr(new_column, .data[[higher_tax]], ~{
        if (.y != "" && !is.na(.y)) {
          match <- .x[str_detect(.x, .y)]
          if (length(match) > 0) match else NA_character_
          } else {
            NA_character_
            }
        })
      ) %>%
    # If the new family entry doesn't start with it already, take the
    # "novel_subfamily_x" part from the beginning of "subfamily..prediction". 
    # and prepend it.
    mutate(new_column = if_else(
      str_starts(new_column, start_string),
      new_column,
      paste(
        str_extract(.data[[lower_tax]], "^([^_]+_){3}[^_]+"),
        new_column,
        sep = "_")
      )) %>%
    # Copy all the original "subfamily..prediction." entries that don't need any 
    # of this treatment 
    mutate(
      new_column = if_else(
        str_count(.data[[lower_tax]], "\\|\\|") < 1,
        .data[[lower_tax]],
        new_column
      )
    ) %>%
    rename(!!sym(taxlevel_to_pick) := new_column)
}
