# Title:    Example imputation with IVEware 
# Author:   Edoardo Costantini
# Created:  2023-02-23
# Modified: 2023-03-28

# Simple example ---------------------------------------------------------------

# Load IVEwareExampleData example data for imputation
load(file = "../data/IVEwareExampleData.rda")

# Define the location of scrlib app
srclib <<- "/Library/Srclib/R"

# Initialize srclib
source(
    file.path(srclib,
        "init.R",
        fsep = .Platform$file.sep
    )
)

# Number of variables with missing values ("imp" in  the output)
sum(colSums(is.na(IVEwareExampleData)) != 0)

# Store the file name containing the imputation instructions
instr <- "IVEwareExampleInstr"

# Perform the imputations
impute(name = instr)

# Define the number of multiple imputations used in the instr script
M <- 5

# Store the outputs
lapply(1:M, function(m){
    putdata(
        name = "IVEwareExampleInstr", 
        dataout = paste0("imp_", m), 
        mult = m)
})

# Load and save each data in internal objects local to this environment
imp_list <- lapply(1:M, function(m) {
    # Load a dataset
    load(file = paste0("imp_", m, ".rda"))

    # Return the dataset "imp" so that it can be stored separately
    imp
})

# Define some rows to check imputations
var_to_check <- "CHOLESTH"
NA_index <- is.na(IVEwareExampleData[, var_to_check])

# Check values are different
cbind(
    OG = IVEwareExampleData[NA_index, var_to_check],
    imp1 = imp_list[[1]][NA_index, var_to_check],
    imp2 = imp_list[[2]][NA_index, var_to_check],
    imp3 = imp_list[[3]][NA_index, var_to_check],
    imp4 = imp_list[[4]][NA_index, var_to_check],
    imp5 = imp_list[[5]][NA_index, var_to_check]
)

# List files in the directory
files <- list.files()

# Identify the ones that are not R scripts or .set (IVEware instructions)
files.to.remove <- files[!grepl(pattern = "*.R|*.set", x = files)]

# Remove files
lapply(files.to.remove, unlink)
