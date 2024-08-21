pick_ambiguous_taxa <- function(vcontact_output, taxlevel_to_correct) {
  # Pick the subfamily that corresponds to the family prediction (or 
  # equvalently pick the genus that corresponds to the subfamily prediction).
  # So run this function twice, once with taxlevel_to_correct="Subfamily", once
  # with taxlevel_to_correct="Genus"
  
  if (taxlevel_to_correct == "Subfamily") {
    higher_tax <- "family..prediction."
    lower_tax <- "subfamily..prediction."
    start_string <- "novel_subfamily"
  }
  if (taxlevel_to_correct == "Genus") {
    higher_tax <- "Subfamily"
    lower_tax <- "genus..prediction."
    start_string <- "novel_genus"
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
    rename(!!sym(taxlevel_to_correct) := new_column) 
}
