### Title:    Replication Zhao Long 2016 - Compare methods
### Author:   Edoardo Costantini
### Created:  2020-02-10
### Modified: 2020-02-14

source("./functions.R")
source("./init.R")

out <- readRDS("../output/pooled_ZL2016-mc-20200309_1936.rds")

# P = 200 condition

(cond_name <- names(out[[1]])[1])

# Average Bias ------------------------------------------------------------

store_sum <- vector("list", parms$dt_rep)

for (i in 1:parms$dt_rep) {
  store_sum[[i]] <- out[[i]][[cond_name]]$cond_bias
}

bias_out <- round(Reduce("+", store_sum)/parms$dt_rep, 3)

bias_b1 <- as.data.frame(t(bias_out))[2]

# Average Coverage --------------------------------------------------------

store_sum <- vector("list", parms$dt_rep)

for (i in 1:parms$dt_rep) {
  store_sum[[i]] <- out[[i]][[cond_name]]$cond_CIco
}

CI_out <- Reduce("+", store_sum)/parms$dt_rep

rownames(CI_out) <- rownames(bias_out)
CI_b1 <- as.data.frame(t(CI_out))[2]

resu <- cbind(bias_b1, CI_b1)
colnames(resu) <- c("bias", "ci")
resu