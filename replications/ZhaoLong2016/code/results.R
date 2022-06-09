### Title:    Replication Zhao Long 2016 - Compare methods
### Author:   Anonymized for peer review
### Created:  2020-02-10
### Modified: 2020-02-14

rm(list = ls())

source("./functions.R")
source("./init.R")

out_old <- readRDS("../output/pooled_ZL2016-mc-20200309_1936.rds") # DO NOT DELETE
  # old version, so it does not have the dt_rep stored inside. If you want to get
  # the resutls you need to write 500 for the store_sum repetitions
  # old valid (pre correct blasso) keep to comapre results with fixed blasso
  # and previous incorrect blasso
out <- readRDS("../output/pooled_ZL2016-mc-20200507_1521.rds")

# Results from simulation -------------------------------------------------
# Old
res_old <- lapply(1:4, function(x){ 
  extract_results(cond_name = names(out_old[[1]])[x], 
                  output = out_old, 
                  dt_rep = 500)
})
names(res_old) <- names(out[[1]]) 
res_old

# New
res <- lapply(1:4, function(x){
  extract_results(cond_name = names(out[[1]])[x], 
                  output = out, 
                  dt_rep = out[[1]]$cond_200_4$parms$dt_rep)
})
names(res) <- names(out[[1]]) 
res

# Results from paper ------------------------------------------------------
results_paper <- list(cond_200_4 = data.frame(bias_p = c(.074, .017, -.023, 
                                                       -.008, -.356, -.206),
                                              ci_p = c(.894, .930, .946,
                                                     .938, .378, .630)),
                      cond_1000_4 = data.frame(bias_p = c(.066, .030, -.030, 
                                                        -.008, -.356, -.206),
                                               ci_p = c(.904, .918, .934,
                                                      .938, .378, .630)),
                      cond_200_20 = data.frame(bias_p = c(.053, -.012, -.037, 
                                                        -.046, -.264, -.224),
                                               ci_p = c(.886, .918, .942,
                                                      .940, .530, .412)),
                      cond_1000_20 = data.frame(bias_p = c(.032, -.003, -.056, 
                                                         -.046, -.264, -.224),
                                                ci_p = c(.914, .934, .956,
                                                       .940, .530, .412))
)

# Compare paper and simulation study --------------------------------------

compare_b1 <- vector("list", length(out[[1]]))
  names(compare_b1) <- names(out[[1]]) 
for (i in 1:length(out[[1]])) {
  compare_b1[[i]] <- cbind(res_old[[i]], res[[i]], results_paper[[i]])
}
compare_b1
