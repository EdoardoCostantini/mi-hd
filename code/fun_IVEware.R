# Title:    Imputation using IVEware software
# Project:  Imputing High Dimensional Data
# Author:   Edoardo Costantini
# Created:  2023-03-20
# Modified: 2023-03-20

impute_IVEware <- function(Z, perform = TRUE, parms = parms){
  
  ## Input: 
  # @Z: dataset w/ missing values
  # @parms: the initialization object parms
  
  ## Example inputs
  # Z = Xy_mis # or
  # parms = parms
  
  ## output: 
  # - a list of chains imputed datasets at iteration iters
  # - per variable list of imputed values
  # - imputation run time
  
  ## body:
  if(perform == TRUE){
    
    tryCatch({
      
    start.time <- Sys.time()

    # Save the current data set as a temp data for processing
    save(Z, file = "Z.rda")

    # Define instruction name
    instr <- "fun_IVEware_instructions"

    # Assemble .set script
    cat(
        "title Multiple imputation; \n",
        "datain Z; \n",
        "dataout Z_imputed; \n",
        "default continuous; \n",
        "minrsqd .01; \n",
        paste0("iterations ", parms$iters, "; \n"),
        paste0("multiples ", parms$mice_ndt, "; \n"),
        "run;",
        file = paste0(instr, ".set")
    )

    # Define the location of scrlib app
    srclib <<- "/Library/Srclib/R"

    # Initialize srclib
    source(
      file.path(srclib,
        "init.R",
        fsep = .Platform$file.sep
      )
    )

    # Perform imputation
    impute(name = instr)

    # Define end time
    end.time <- Sys.time()

    # Store the outputs
    lapply(1:parms$mice_ndt, function(m) {
      putdata(
        name = instr,
        dataout = paste0("imp_", m),
        mult = m
      )
    })

    # Load and save each data in internal objects local to this environment
    imp_list <- lapply(1:parms$mice_ndt, function(m) {
      # Load a dataset
      load(file = paste0("imp_", m, ".rda"))

      # Return the dataset "imp" so that it can be stored separately
      imp
    })

    # List files in the directory
    files <- list.files()

    # Identify the ones that are not R scripts or .set (IVEware instructions)
    files.to.remove <- files[!grepl(pattern = "*.R|*.set", x = files)]

    # Remove files
    lapply(files.to.remove, unlink)

    # Store results
    imp_IVEware_dats <- imp_list
    imp_IVEware_imps <- NULL
    imp_IVEware_time <- difftime(end.time, start.time, units = "mins")
    
    return(list(dats = imp_IVEware_dats,
                imps = imp_IVEware_imps,
                time = imp_IVEware_time))
    
    
    }, error = function(e){
      err <- paste0("Original Error: ", e)
      print(err)
      return(list(dats = NULL,
                  imps = NULL,
                  time = NULL)
      )
    }
    )

  } else {
    return(list(dats = NULL,
                imps = NULL,
                time = NULL))
  }
}
