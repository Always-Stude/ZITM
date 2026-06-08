df <- data.frame(
  kon = runif(n = 10000, min = 0.1, max = 3),    
  koff = runif(n = 10000, min = 0.1, max = 3),   
  ksyn = runif(n = 10000, min = 10, max = 30), 
  bf = NA,                                     
  bs = NA,
  size = NA,
  prob = NA
)

rownames(df)=paste0("gene", 1:10000)

df$bf=df$kon;df$bs=df$ksyn/df$koff
df$size=df$bf
df$prob=1/(1+df$bs)

#TM
N=10
result_matrix <- matrix(NA, nrow = 10000, ncol = N)
rownames(result_matrix) <- rownames(df)  
colnames(result_matrix) <- paste0("cell", 1:N)

for (i in 1:10000) {
  kon_i <- df$kon[i]
  koff_i <- df$koff[i]
  ksyn_i <- df$ksyn[i]
  
  for (j in 1:N) {
    c <- rbeta(1, shape1 = kon_i, shape2 = koff_i)
    lambda <- c * ksyn_i
    result_matrix[i, j] <- rpois(1, lambda = lambda)
  }
}

#NB
N=10
result_matrix <- matrix(NA, nrow = 10000, ncol =N)
rownames(result_matrix) <- rownames(df)  
colnames(result_matrix) <- paste0("cell", 1:N)

for (i in 1:10000) {
  size_i <- df$size[i]
  prob_i <- df$prob[i]
  
  for (j in 1:10) {
    
    result_matrix[i, j] <- rnbinom(n = 1, size = size_i, prob = prob_i)
  }
}


#dropout
result_matrix=t(apply(as.matrix(result_matrix), 1, function(row) {
  
  if (all(row == 0)) {
    return(rep(0, length(row)))
  }
  
  else {
    non_zero = row[row != 0]   
    u = log(mean(non_zero))    
    p0 = exp(-0.1 * u^2)       
    bernoulli_success = rbinom(length(row), size = 1, prob = 1-p0)
    ss = row * bernoulli_success
  }
}))
