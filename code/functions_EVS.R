### Title:    Functions to understand and clean EVS data
### Author:   Edoardo Costantini
### Created:  2020-06-22

# Function to tabulate the answers of a variable by country
tab_country <- function(dt, country, var_name){
  # dt = int.df
  # country = int.df$country
  # var_name = "v173"
  out.put <- lapply(
    unique(country), function(x){
      cs_var <- (dt[country == x, ][, var_name])
      cs_tab <- summary(factor(cs_var, levels = val_labels(cs_var)))
      cs_out <- as.data.frame(cs_tab)
      colnames(cs_out) <- names(which(val_labels(country) == x))
      return(cs_out)
    }
  )
  fun.out <- list(t(do.call(cbind, out.put)))
  names(fun.out) <- paste0(var_name, ": ", 
                           var_label(dt[, var_name]))
  return( fun.out )
}

# Example Use
# tab_country(int.df, int.df$country, "v173")

# Clean data
# Put inportant information first, convert to formats you like

clean_up <- function(dt, id, country, age, inc, ord, edu){
  dt.character <- sapply(dt[, id], as.character)
  # dt.country   <- sapply(dt[, country], factor,
  #                        levels = val_labels(country),
  #                        labels = names(val_labels(country)))
  dt.country <-  factor(dt[, country], 
                        levels = val_labels(dt[, country]), 
                        labels = names(val_labels(dt[, country])))
  dt.numeric   <- sapply(dt[, c(age, inc)], as.numeric)
  dt.factors   <- dt[, c(ord, edu)]
  
  data.frame(id = dt.character, country = dt.country, dt.numeric, dt.factors)
}
