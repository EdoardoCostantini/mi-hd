### Title:    Generating lavan model text
### Author:   Edoardo Costantini
### Created:  2020-06-24

## To get a model, create the object parms$z_m_id
## and run this script

# Save location
model_loc_path <- "../txt/"
model_name <- "lavaan_model_sat.txt"

# Define varibale names
var_names <- paste0("z", parms$z_m_id)

# Text objects ------------------------------------------------------------

# Means
head_means <- "# Intercepts\n"
all_means <- paste0(var_names, " ~ ", "1")
all_means <- paste(all_means, collapse = "\n")

# Variances
head_vars <- "# Variances \n"
all_vars <- paste0(var_names, " ~~ ", var_names)
all_vars <- paste(all_vars, collapse = "\n")

# Coivariances
head_covs <- "# Covariances \n"
all_covs <- combn(var_names, 2)
all_covs <- apply(all_covs, 2, paste0, collapse = " ~~ ")
all_covs <- paste(all_covs, collapse = "\n")

# Compile model file ------------------------------------------------------

# Start and intercepts
cat(paste(head_means,
          all_means),
    file = paste0(model_loc_path, model_name),
    sep = "\n")

# Add variances
cat(paste(head_vars, all_vars),
    file = paste0(model_loc_path, model_name),
    sep = "\n",
    append = TRUE)

# Add covariances and end
cat(paste(head_covs, 
          all_covs),
    file = paste0(model_loc_path, model_name),
    sep = "\n",
    append = TRUE)

# # Read model --------------------------------------------------------------
# 
# sat_model <- read.table(paste0(model_loc_path, model_name),
#                         as.is = TRUE)$V1
# class(sat_model)
# 
# # Example
# sta_model <- '
#   # Interecepts
#   z1 ~ 1
#   z2 ~ 1
#   z3 ~ 1
#   z4 ~ 1
#   z5 ~ 1
#   z11 ~ 1
#   z12 ~ 1
#   z13 ~ 1
#   z14 ~ 1
#   z15 ~ 1
#   
#   # Variances
#   z1  ~~ z1 
#   z2  ~~ z2 
#   z3  ~~ z3 
#   z4  ~~ z4 
#   z5  ~~ z5 
#   z11 ~~ z11
#   z12 ~~ z12
#   z13 ~~ z13
#   z14 ~~ z14
#   z15 ~~ z15
#   
#   # Covariances
#   z1  ~~ z1
#   z1  ~~ z2
#   z1  ~~ z3
#   z1  ~~ z4
#   z1  ~~ z5
#   z1  ~~ z11
#   z1  ~~ z12
#   z1  ~~ z13
#   z1  ~~ z14
#   z1  ~~ z15
#   z2  ~~ z3
#   z3  ~~ z3 
#   z4  ~~ z4 
#   z5  ~~ z5 
#   z11 ~~ z11
#   z12 ~~ z12
#   z13 ~~ z13
#   z14 ~~ z14
#   z15 ~~ z15
#   '
# 
# class(sta_model)