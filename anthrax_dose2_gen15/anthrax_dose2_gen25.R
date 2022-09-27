## Mutation analysis
## From a homogeneous infectious does, simulate emergence of SNPs in expanding population.

## Genome size and mutation rate replicate anthrax (Bacillus anthracis)

## parameters governing simulation
genome_size <- 5200000
mutation_rate <- 0.00000000083
generations <- 20
dose_size <- 2
n_sim <- 100

#########################################
res.dir <- paste0("sim_dose", dose_size, "_gen", generations)
dir.create(res.dir)
genome_mut_rate <- genome_size*mutation_rate

## Record number of genomes with each SNP count vs initial dose
for (sim in 1:n_sim) {
  
  # Next analysis - each new genotype becomes a new row
  dose_var <- data.frame(var = "var", SNP = 0)
  results <- as.data.frame(matrix(0, 1, generations))
  names(results) <- paste0("gen_", 1:generations)
  results[1,1] <- dose_size
  results <- cbind(results, dose_var)
  
  
  for (c in 1:(generations-1)) {
    
    # use a for loop to split rows with large numbers of genotypes
    for (r in 1:nrow(results)) {
      if (results[r, c] > 1000000000) {
        # half the number of genotypes in a specific row, then make a new row
        # OK to divide by 2 as all numbers above threshold will be even numbers (because produced by a doubling process)
        results[r, c] <- results[r, c] / 2
        new.row <- results[r,]
        # for new row, remove genotypes in prev generations - they are counted in
        # original row
        new.row[1:(c-1)] <- 0
        results <- rbind(results, new.row)
      }
    }
    
    # now a second loop does replication and division
    # (copes with previous for loop potentially adding rows in the same generation)
    for (r in 1:nrow(results)) {
      
      if (results[r, c] != 0) {
        # use rpois to determine distribution of SNPs
        SNPs <- rpois(results[r, c], genome_mut_rate)
        
        ################################
        # Deal with no SNP genomes first
        SNPs_0 <- sum(SNPs == 0)
        if (SNPs_0 > 0) {
          offspring_non_mut <- SNPs_0 * 2
          results[r, (c+1)] <- offspring_non_mut
        }
        
        ################################
        # Deal with 1 SNP genomes (2 of each mutant genome in next gen)
        SNPs_1 <- sum(SNPs == 1)
        if (SNPs_1 > 0) {
          new.rows <- as.data.frame(matrix(0, SNPs_1, generations))
          # These genomes first appear in gen + 1 at frequency of 2
          new.rows[, (c+1)] <- 2
          # sample positions to be SNPs
          new.rows$var <- paste(results$var[r],
                                round(runif(n = SNPs_1, 1, genome_size), 0),
                                sep = "_")
          new.rows$SNP <- results$SNP[r] + 1 # N SNPs increases by 1
          
          names(new.rows) <- names(results)
          results <- rbind(results, new.rows)
        }
        
        ################################
        # Deal with 2 SNP genomes (2 of each mutant genome in next gen)
        SNPs_2 <- sum(SNPs == 2)
        if (SNPs_2 > 0) {
          new.rows <- as.data.frame(matrix(0, SNPs_2, generations))
          # These genomes first appear in gen + 1 at frequency of 2
          new.rows[, (c+1)] <- 2 
          # sample positions to be SNPs
          new.rows$var <- paste(results$var[r],
                                round(runif(n = SNPs_2, 1, genome_size), 0),
                                round(runif(n = SNPs_2, 1, genome_size), 0),
                                sep = "_")
          new.rows$SNP <- results$SNP[r] + 2 # N SNPs increases by 2
          
          names(new.rows) <- names(results)
          results <- rbind(results, new.rows)
          
        }
        
        ################################
        # Deal with 3 SNP genomes (2 of each mutant genome in next gen)
        SNPs_3 <- sum(SNPs == 3)
        if (SNPs_3 > 0) {
          new.rows <- as.data.frame(matrix(0, SNPs_3, generations))
          # These genomes first appear in gen + 1 at frequency of 2
          new.rows[, (c+1)] <- 2 
          # sample positions to be SNPs
          new.rows$var <- paste(results$var[r],
                                round(runif(n = SNPs_3, 1, genome_size), 0),
                                round(runif(n = SNPs_3, 1, genome_size), 0),
                                round(runif(n = SNPs_3, 1, genome_size), 0),
                                sep = "_")
          new.rows$SNP <- results$SNP[r] + 3 # N SNPs increases by 2
          
          names(new.rows) <- names(results)
          results <- rbind(results, new.rows)
        }
        
        
        ################################
        # Deal with 4 SNP genomes (2 of each mutant genome in next gen)
        SNPs_4 <- sum(SNPs == 4)
        if (SNPs_4 > 0) {
          new.rows <- as.data.frame(matrix(0, SNPs_4, generations))
          # These genomes first appear in gen + 1 at frequency of 2
          new.rows[, (c+1)] <- 2 
          # sample positions to be SNPs
          new.rows$var <- paste(results$var[r],
                                round(runif(n = SNPs_4, 1, genome_size), 0),
                                round(runif(n = SNPs_4, 1, genome_size), 0),
                                round(runif(n = SNPs_4, 1, genome_size), 0),
                                sep = "_")
          new.rows$SNP <- results$SNP[r] + 4 # N SNPs increases by 2
          
          names(new.rows) <- names(results)
          results <- rbind(results, new.rows)
        }
        
      }
    }
    
    # In every generation generation, save results
    write.csv(results, paste0(res.dir, "/sim", sim, "_dose", dose_size, ".csv"),
              quote = F, row.names = F)
  }
  
}

