# Title:    Initialization for Blade Computing run
# Author:   Edoardo Costantini
# Created:  2020-05-16
# Modified:  2022-09-15

# Packges to install:

  # 1. Install all packages you can
  source("./init_general.R")
  install.packages(pack_list)
  
  # 2. Install PcAux
  install.packages("devtools")
  devtools::install_github("PcAux-Package/PcAux/source/PcAux")
  
  # 3. Install blasso package on a Windows computer:
  # Install w/ point and click interface. Use the gz.tar file in the package 
  # folder in the impute_HD folder or
  install.packages("../input/blasso_0.3.tar.gz", repos = NULL, type = "source")