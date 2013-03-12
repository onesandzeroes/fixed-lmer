# Overview ##################
# TODO ##############
# GetVarInfo fails when the contrasts have just been set as, for example
# running the example analysis of the lme4::cake dataset fails because the
# contrasts are listed as:
# $recipe: "contr.treatment"
# Think this can be diagnosed by checking the type/class of each element of
# contrast_info, as it will be character for these kinds of default contrasts,
# and numeric for manually specified contrasts
# Should be easy enough to fix, but need to make sure that original contrast
# names get preserved
# Seems to be failing on cake$temperature's polynomial contrasts, as these
# don't provide labels for each level in quite the same way

# Functions #######
GetVarInfo <- function(mod, effect_name) {
  # Split the effect name into its individual components
  # NB: This will fail if the variable names contain colons!
  all_vars <- unlist(strsplit(effect_name, ":"))
  # Create a mapping from the effect name (which may have a contrast name/number
  # tacked on the end) to the variable name
  effect_map <- list()
  for (varname in all_vars) {
    effect_map[[varname]] <- varname
  }
  # Variables that use contrasts set using the contrasts() attribute have
  # slightly modified names (since the name/number of the contrast is added at
  # the end), so need to grab and match them
  contrast_info <- attr(mod@X, "contrasts")
  contrast_vars <- names(contrast_info)
  for (var in contrast_vars) {
    contrast_info[[var]] = contrasts(mod@frame[[var]])
  }
  contrast_values <- list()
  for (c_var in contrast_vars) {
    var_search <- grep(c_var, all_vars)
    # grep returns a 0-length integer vector if the term is not found, so:
    if (length(var_search) > 0) {
      full_effect_name <- all_vars[var_search[1]]
      # Name/number of the contrast will be the part of the variable name
      # that's not in c_var
      contrast_name <- unlist(strsplit(full_effect_name, c_var))[2]
      # Check to see if the contrast name can be inted, in which case the
      # contrasts (probably) didn't have names set and should be accessed
      # by number
      as_int <- as.integer(contrast_name)
      if ( ! is.na(as_int[1])) {
        contrast_name <- as_int
      }
      contrast_values[[c_var]] <- contrast_info[[c_var]][, contrast_name]
      effect_map[[full_effect_name]] <- c_var
    }
  }
  var_values <- list()
  print("Contrast vars:")
  print(contrast_vars)
  print("Effect map:")
  print(effect_map)
  for (varname in effect_map) {
    if (varname %in% contrast_vars) {
      labels <- names(contrast_values[[varname]])
      if (! is.null(labels)) {
        # Keep the level labels for labelled contrast variables
        var_values[[paste0(varname, "_labels")]] <- labels
      } else {
        var_values[[varname]] <- contrast_values[[varname]]
      }
    } else {
      var_values[[varname]] <- unique(mod@frame[[varname]])
    }
  }
  res_list <- list(
    var_values=var_values,
    all_vars=all_vars,
    contrast_info=contrast_info,
    contrast_vars=contrast_vars,
    contrast_values=contrast_values,
    effect_map=effect_map
  )
  return(res_list)
}

# Takes a row from the fixef_df created by FixedInteraction, and returns
# a string (or maybe call/expression when I figure out how to turn these things
# into functions) containing the calculation for that particular coefficient
CreateCoefficient <- function(df_row) {
  # Columns 1-4 are estimate, s.e., t, and effect name, 5 onwards
  # are the effects
  effects <- colnames(df_row)[5:ncol(df_row)]
  est <- as.character(df_row$Estimate)
  effect_logical <- as.logical(df_row[1, 5:ncol(df_row)])
  effects_to_keep <- c(est, effects[effect_logical])
  coef_string <- paste0(effects_to_keep, collapse=" * ")
  df_row$coef_string <- coef_string
  return(df_row)
}

FixedInteraction <- function(mod, effect_name) {
  cat("Warning: this function is still in progress and untested: \n")
  cat("Correct results not guaranteed")
  require(plyr)
  var_info <- GetVarInfo(mod, effect_name)
  var_values <- var_info[["var_values"]]
  effect_map <- var_info[["effect_map"]]
  contrast_values <- var_info[["contrast_values"]]
  print("Contrast values:")
  print(contrast_values)
  print("Variable values:")
  print(var_values)
  # By this point, we have all variables that enter into the interaction separated out,
  # as well as the values of each variable's contrast coefficients
  level_grid <- expand.grid(var_values)
  label_cols <- colnames(level_grid)[grep("_labels", colnames(level_grid))]
  for (label_col in label_cols) {
    base_name <- strsplit(label_col, "_labels")[[1]][1]
    level_grid[[base_name]] <- contrast_values[[base_name]][
      match(
        level_grid[[label_col]],
        names(contrast_values[[base_name]])
      )
      ]
  }
  print("Level grid: ")
  print(head(level_grid))
  fixef_mat <- summary(mod)@coefs
  # Get the table of fixed effects, estimates, sds and t.values- currently
  # the effect names are stored as the row names of this dataframe
  fixef_df <- data.frame(fixef_mat)
  fixef_df$effect <- rownames(fixef_df)
  model_intercept <- as.character(fixef_df$Estimate[fixef_df$effect == "(Intercept)"][1])
  rownames(fixef_df) <- NULL
  # For contrasts, we need to grab the effects specific
  # to that contrasts, and therefore need the full effect name,
  # found in the names() of effect_map
  for (ind in seq_along(effect_map)) {
    effect_name <- names(effect_map)[ind]
    var_name <- effect_map[[effect_name]]
    fixef_df[[var_name]] <- FALSE
    fixef_df[[var_name]][grep(effect_name, fixef_df$effect)] <- TRUE
  }
  # Need to keep only effects made up of combinations of the interaction variables,
  # if any other variables are present, that effect is not relevant to the
  # current interaction
  row_effects <- strsplit(fixef_df$effect, ":")
  contains_only_ivars <- sapply(
    row_effects,
    function(x) {
      all(x %in% names(effect_map))
    })
  fixef_df <- fixef_df[contains_only_ivars, ]
  # Create coefficient strings, still not sure how to turn these into
  # functions, but it seems like the obvious next step
  fixef_df <- ddply(
    fixef_df,
    .(effect),
    CreateCoefficient
  )
  final_formula_string <- paste0(
    c(model_intercept, fixef_df$coef_string),
    collapse=" + "
  )
  print(final_formula_string)
  EffectFunction <- function(level_grid_df) {
    parsed <- parse(text=final_formula_string)
    effect_val <- eval(parsed, envir=level_grid_df)
    return(effect_val)
  }
  level_grid$effect_val <- EffectFunction(level_grid)
  res_list <- list(
                   results=level_grid,
                   predict_func=EffectFunction
                   )
  return(res_list)
}



