### Title:    Replication Zhao Long 2016 - Compare methods
### Author:   Edoardo Costantini
### Created:  2020-02-10
### Modified: 2020-02-14

rm(list = ls())

source("./functions.R")
source("./init.R")

out <- readRDS("../output/pooled_ZL2016-mc-20200309_1936.rds") # last valid DO NOT DELETE
# old version, so it does not have the dt_rep stored inside. If you want to get
# the resutls you need to write 500 for the store_sum repetitions
out <- readRDS("../output/pooled_ZL2016-mc-20200319_1032.rds") # last valid
out <- readRDS("../output/pooled_ZL2016-mc-20200504_1703.rds") # last valid

# P = 200 condition

(cond_name <- names(out[[1]])[1])

# Average Bias ------------------------------------------------------------

store_sum <- vector("list", out[[1]]$cond_200_4$parms$dt_rep)

for (i in 1:parms$dt_rep) {
  store_sum[[i]] <- out[[i]][[cond_name]]$cond_bias
}

bias_out <- round(Reduce("+", store_sum)/parms$dt_rep, 3)

bias_b1 <- as.data.frame(t(bias_out))[2]

# Average Coverage --------------------------------------------------------

store_sum <- vector("list", out[[1]]$cond_200_4$parms$dt_rep)

for (i in 1:parms$dt_rep) {
  store_sum[[i]] <- out[[i]][[cond_name]]$cond_CIco
}

CI_out <- Reduce("+", store_sum)/parms$dt_rep

rownames(CI_out) <- rownames(bias_out)
CI_b1 <- as.data.frame(t(CI_out))[2]

resu <- cbind(bias_b1, CI_b1)
colnames(resu) <- c("bias", "ci")
resu

# Compare to paper
paper <- data.frame(bias = resu$bias,
                    bias_p = c(.074, .017, -.023, -.023, -.008, -.356, -.206),
                    ci = resu$ci,
                    ci_p = c(.894, .930, .946, .946,  .938, .378, .630))
rownames(paper) <- rownames(resu)
paper
