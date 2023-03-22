# Title:    Imputation using IVEware software
# Project:  Imputing High Dimensional Data
# Author:   Edoardo Costantini
# Created:  2023-03-20
# Modified: 2023-03-22

impute_IVEware <- function(Z, minR2 = 0.01, rep_status = 1, perform = TRUE, parms = parms){
  
  ## Input: 
  # @Z: dataset w/ missing values
  # @parms: the initialization object parms
  
  ## Example inputs
  # Z = Xy_mis # or
  # minR2 = 0.01
  # parms = parms
  # rep_status = 1
  
  ## output: 
  # - a list of chains imputed datasets at iteration iters
  # - per variable list of imputed values
  # - imputation run time
  
  ## body:
  if(perform == TRUE){
    
    tryCatch({

    # Define a temp folder name
    tempfold <- paste0("IVEware-temp-", rep_status)

    # Create temp folder for processing
    dir.create(tempfold)

    # Move to the tempfold
    setwd(tempfold)

    # Identify variables under imputations in Z
    target.Z <- names(which(colSums(is.na(Z)) != 0))

    # Store a new version of Z
    Z_IVEware <- Z

    # Store current NA option
    default.na.action <- options("na.action")$na.action

    # Make sure options for na is set to pass
    options(na.action = "na.pass")

    # Make Z a design matrix as a data.frame
    Z_IVEware <- data.frame(model.matrix(~., data = Z_IVEware)[, -1])

    # Sanitize the names of the factors for easy processing in IVEware (required for evs)
    colnames(Z_IVEware) <- paste0("v", 1:ncol(Z_IVEware))

    # Identify variables under imputations in Z_IVEware
    target.ZIVEware <- names(which(colSums(is.na(Z_IVEware)) != 0))

    # Save the current data set as a temp data for processing
    save(Z_IVEware, file = "Z_IVEware.rda")

    # Define instruction name
    instr <- paste0("fun_IVEware_instructions")

    # Assemble .set script
    cat(
        "title Multiple imputation; \n",
        "datain Z_IVEware; \n",
        "dataout Z_IVEware_imputed; \n",
        "default continuous; \n",
        paste0("minrsqd ", minR2, "; \n"),
        paste0("iterations ", parms$iters, "; \n"),
        paste0("maxpred ", ncol(Z)-1, "; \n"),
        paste0("multiples ", parms$mice_ndt, "; \n"),
        "print all; \n",
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

    # Define a log file to store all the console output
    IVEware_log <- file(paste0(instr, ".txt"))

    # Define start time
    start.time <- Sys.time()

    # Sink the log
    # sink(IVEware_log, append = TRUE, type = "output")

    # Perform imputation
    impute(name = instr)

    # Close the full log
    # closeAllConnections()

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

      # Make a copy of z
      Z_temp <- Z

      # Return the dataset "imp" so that it can be stored separately
      Z_temp[, target.Z] <- imp[, target.ZIVEware]

      # Return
      Z_temp
    })

    # Move back to previous wd location
    setwd("../")

    # Remove temp folder
    unlink(tempfold, recursive = TRUE)

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
