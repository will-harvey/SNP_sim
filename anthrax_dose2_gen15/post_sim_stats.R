### Population stats

## Calculate some summary statistics for evolved populations
## In each simulation, in each generation

## parameters relevant to  simulation
res.dir <- "sim_dose2_gen20"
save_as <- "stats_dose2_gen20.Rdata"

n_sim <- 100
dose_size <- 2
gen_to_sample <- 20
gen.max <- max(gen_to_sample) - 1 # this is used to determine which files to read in

#########################################
results_files <- list.files(path = res.dir)
stat <- c()

for (i in 1:n_sim) {
  dat <- read.csv(paste0(res.dir, "/sim", i, "_dose", dose_size, ".csv"))
  
  sim <- data.frame(sim = i,
                    dose = dat[1,1],
                    gen = 1:gen_to_sample,
                    pop.size = NA,
                    genomes = NA,
                    max_snp = NA,
                    prop_0 = NA,
                    prop_1 = NA,
                    prop_2 = NA,
                    prop_3 = NA,
                    prop_4 = NA,
                    prop_5 = NA,
                    prop_6 = NA,
                    prop_7 = NA,
                    prop_8 = NA,
                    prop_9 = NA,
                    prop_10 = NA)
  
  for (j in 1:gen_to_sample) {
    # subset to get genomes in gen only
    temp <- dat[dat[,j] > 0,]
    pop.size <- sum(temp[,j])
    sim$pop.size[sim$gen == j] <- pop.size
    # number of rows is now the number of unique genomes
    sim$genomes[sim$gen == j] <- length(unique(temp$var))
    sim$max_snp[sim$gen == j] <- max(temp$SNP)
    
    # proportions of pop with each no. of snps
    sim$prop_0[sim$gen == j] <- sum(temp[temp$SNP == 0, j]) / pop.size
    sim$prop_1[sim$gen == j] <- sum(temp[temp$SNP == 1, j]) / pop.size
    sim$prop_2[sim$gen == j] <- sum(temp[temp$SNP == 2, j]) / pop.size
    sim$prop_3[sim$gen == j] <- sum(temp[temp$SNP == 3, j]) / pop.size
    sim$prop_4[sim$gen == j] <- sum(temp[temp$SNP == 4, j]) / pop.size
    sim$prop_5[sim$gen == j] <- sum(temp[temp$SNP == 5, j]) / pop.size
    sim$prop_6[sim$gen == j] <- sum(temp[temp$SNP == 6, j]) / pop.size
    sim$prop_7[sim$gen == j] <- sum(temp[temp$SNP == 7, j]) / pop.size
    sim$prop_8[sim$gen == j] <- sum(temp[temp$SNP == 8, j]) / pop.size
    sim$prop_9[sim$gen == j] <- sum(temp[temp$SNP == 9, j]) / pop.size
    sim$prop_10[sim$gen == j] <- sum(temp[temp$SNP == 10, j]) / pop.size
    
    sim$prop_0[sim$gen == j] <- round(sim$prop_0[sim$gen == j], 3)
    sim$prop_1[sim$gen == j] <- round(sim$prop_1[sim$gen == j], 3)
    sim$prop_2[sim$gen == j] <- round(sim$prop_2[sim$gen == j], 3)
    sim$prop_3[sim$gen == j] <- round(sim$prop_3[sim$gen == j], 3)
    sim$prop_4[sim$gen == j] <- round(sim$prop_4[sim$gen == j], 3)
    sim$prop_5[sim$gen == j] <- round(sim$prop_5[sim$gen == j], 3)
    sim$prop_6[sim$gen == j] <- round(sim$prop_6[sim$gen == j], 3)
    sim$prop_7[sim$gen == j] <- round(sim$prop_7[sim$gen == j], 3)
    sim$prop_8[sim$gen == j] <- round(sim$prop_8[sim$gen == j], 3)
    sim$prop_9[sim$gen == j] <- round(sim$prop_9[sim$gen == j], 3)
    sim$prop_10[sim$gen == j] <- round(sim$prop_10[sim$gen == j], 3)
  }
  stat <- rbind(stat, sim)
}

save(stat, file = save_as)

