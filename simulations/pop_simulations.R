getMultiplicativeFitness <- function(selection_coefficient, number_of_mutations){
  (1-selection_coefficient)^number_of_mutations
}

getEpistaticFitness <- function(selection_coefficient, number_of_mutations, alpha){
  exp(alpha*(1-1/getMultiplicativeFitness(selection_coefficient, number_of_mutations)))
} ###this function is not used, but I have it in case we want to try a different fitness landscape

getMutationProbability <- function(mutation_rate, number_of_mutations){
  dpois(number_of_mutations, mutation_rate)
}


### The following function is no longer because it is both too slow and requires too much memory
getTransitionMatrix <- function(selection_coef, mutation_rate, number_of_cpgs){ 
  fitness_vector <- sapply(0:number_of_cpgs, function(xx) getMultiplicativeFitness(selection_coef, xx))
  mutation_probability_vector <- sapply(0:number_of_cpgs, function(xx) getMutationProbability(mutation_rate, xx))
  
  m_matrix <- matrix(NA, nrow = number_of_cpgs+1, ncol = number_of_cpgs+1)
  for (i in 2:(number_of_cpgs+1)){
    m_matrix[, i] <- c(rep(0, i-1), mutation_probability_vector[1:(length(mutation_probability_vector)-(i-1))])
    m_matrix 
  }
  m_matrix[, 1] <- mutation_probability_vector
  transition_matrix <- matrix(NA, nrow = number_of_cpgs+1, ncol = number_of_cpgs+1)
  for (i in 1:(number_of_cpgs+1)){ #this loop gives  the multiplication of the M matrix by the W matrix because W is diagonal
    transition_matrix[, i] <- fitness_vector[i]*m_matrix[, i]
    transition_matrix 
  }
  transition_matrix
}

library(Rcpp)
## Function for computing the matrix product M * WP where WP is a vector formed by W*P
## and M is a lower triangular vector with a special structure, which means we can store the entries as vector
cppFunction('NumericVector transmult(NumericVector m, NumericVector wp) {
  int n = m.length();
  NumericVector out(n);
  for (int i=0; i<n; i++) {
    double total=0;
    for (int j=0; j<=i; j++) {
      total += wp[j]*m[i-j];
    }
    out[i] = total;
  }
  return out;
}')


library(fftwtools)
###the following function is what we actually use because it is much faster once the number of sites gets somewhat large 
computeConvolution <- function(x, y){
  n1 <- n_sites
  x <- c(rep.int(0, n1), x)
  n <- length(y <- c(y, rep.int(0, n_sites)))
  
  res <- fftw_r2c(x, HermConj=1) * Conj(fftw_r2c(y, HermConj=1))
  res <- fftw_c2c(res, inverse = 1)
  Re(res)/n
  #Re(fftw_c2c(fftw_r2c(x, HermConj=1) * Conj(fftw_r2c(y, HermConj=1)), inverse = 1))/n
}


runSimulationWithDrift <- function(n_vector_a, n_vector_b, n_of_possible_mutations, 
                                   fitness_vector, mutation_probability_vector_a, 
                                   mutation_probability_vector_b) {
  n_a <- n_vector_a[1] # The initial number of indiviuals with epiallele A
  n_b <- n_vector_b[1]
  total_number_of_individuals <- n_a + n_b
  
  #each entry in n_vector_a and n_vector_b corresponds to the number of individuals with a given number of mutations 
  #of type A and B, respectively 
  n_vector <- c(n_vector_a, n_vector_b) 
  p_a <- n_a/total_number_of_individuals
  p_b <- 1 - p_a
  
  while (p_b > 0 & p_b < 1){
    cpg_vector_a <- n_vector[1:n_of_possible_mutations]/n_a
    cpg_vector_b <- n_vector[(n_of_possible_mutations+1):(2*n_of_possible_mutations)]/n_b
    
    w_bar_a <- sum(cpg_vector_a*fitness_vector)
    w_bar_b <- sum(cpg_vector_b*fitness_vector)
    
    #cpg_vector_a_new <- (1/w_bar_a)*(transition_matrix_a %*% cpg_vector_a)
    #cpg_vector_b_new <- (1/w_bar_b)*(transition_matrix_b %*% cpg_vector_b)
    cpg_vector_a_times_w <- fitness_vector*cpg_vector_a
    cpg_vector_b_times_w <- fitness_vector*cpg_vector_b
    cpg_vector_a_new <- (1/w_bar_a)*computeConvolution(mutation_probability_vector_a, rev(cpg_vector_a_times_w))[1:n_of_possible_mutations]
    cpg_vector_b_new <- (1/w_bar_b)*computeConvolution(mutation_probability_vector_b, rev(cpg_vector_b_times_w))[1:n_of_possible_mutations]
    
    cpg_vector_a_new[abs(cpg_vector_a_new) < .Machine$double.eps] <- 0
    cpg_vector_b_new[abs(cpg_vector_b_new) < .Machine$double.eps] <- 0
    
    fitnesses_a <- cpg_vector_a_new * fitness_vector
    fitnesses_b <- cpg_vector_b_new * fitness_vector
    
    w_bar_a <- sum(fitnesses_a)
    w_bar_b <- sum(fitnesses_b)
    w_bar <- p_a*w_bar_a + p_b*w_bar_b
    
    n_vector_new <- rmultinom(1, total_number_of_individuals, (1/w_bar)*c(p_a*fitnesses_a, p_b*fitnesses_b))
    n_b <- sum(n_vector_new[(n_of_possible_mutations+1):(2*n_of_possible_mutations)])
    n_a <- total_number_of_individuals - n_b
    p_b <- n_b/total_number_of_individuals
    p_a <- 1 - p_b
    
    
    n_vector <- n_vector_new
    n_vector
    n_a
    n_b
    p_a
    p_b
  }
  p_b
}


## Magic constants
scaling_factor <- 1 #this can be used to adjust population size and mutation rate and make the simulation faster
#we use the scaling factor because the mut rates were estimated from human effective pop size 
#and the s_hets correspond to human effective pop size too

s_hets <- scaling_factor*c(0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01) #this is for CpGs
#for H3K36me3, because we are talking about coding regions: s_hets <- scaling_factor*c(0.001, 0.002, 0.003, 0.004, 0.005, 0.01, 0.05)
nPromoters <- c(1, 50, 100, 200)
nCpGs <- c(80) #this is roughly the average across all proximal promoters (for H3K36me3 it is 135)                                      
possible_values <- expand.grid(s_hets = s_hets, nPromoters = nPromoters, nCpGs = nCpGs)

ii <- Sys.getenv("SGE_TASK_ID")
if(!is.na(ii))
  q(save = "no")
RNGkind("L'Ecuyer-CMRG")
set.seed(ii)

this_nCpGs <- possible_values[ii, "nCpGs"]
this_s <- possible_values[ii, "s_hets"]
this_nPromoters <- possible_values[ii, "nPromoters"]

n_sites <- this_nPromoters * this_nCpGs

## Magic constants which are fixed
## Mutation rates
U_a <- scaling_factor*7.971193e-08 #for H3K36me3: U_a <- scaling_factor*1.4e-08 ##see Methods of the paper for estimation
U_b <- scaling_factor*2.474832e-07 #for H3K36me3: U_b <- scaling_factor*1.14e-08
nIndividuals <- 10000/scaling_factor

## Starting distribution                                        
n_vector_a <- c(0.5 * nIndividuals, rep(0, n_sites)) #in this and the following line the 0.5 can be modified depending on the desired starting frequencies
n_vector_b <- c(0.5 * nIndividuals, rep(0, n_sites))
fitness_vector <- sapply(0:n_sites, function(xx) getMultiplicativeFitness(this_s, xx))

mutation_probability_vector_a <- sapply(0:n_sites, function(xx) getMutationProbability(U_a*n_sites, xx))
mutation_probability_vector_b <- sapply(0:n_sites, function(xx) getMutationProbability(U_b*n_sites, xx))

eq_freq <- replicate(150, runSimulationWithDrift(n_vector_a, n_vector_b, n_sites+1, 
                                                 fitness_vector, mutation_probability_vector_a, 
                                                 mutation_probability_vector_b))
fixation_table <- table(eq_freq)

## Saving
file_name <- sprintf("initial1_%g_%g_%g_10000_individuals.rda", this_s, this_nCpGs, this_nPromoters) 
save(fixation_table, file = file_name)


####
##trials <- 10
###system.time({
#  eq_freq <- foreach(icount(150), .combine=c) %dopar% {
#   runSimulationWithDrift(n_vector_a, n_vector_b, n_sites+1, 
#                           fitness_vector, mutation_probability_vector_a, 
#                           mutation_probability_vector_b)
#  }
#})







