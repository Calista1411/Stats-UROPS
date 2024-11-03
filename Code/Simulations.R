#Finding stationary distribution for a simulated first order Markov chain
x=matrix(runif(16),nrow=4)
rowsum=apply(x,1,sum)
a=diag(1/rowsum,4,4)
P=a%*%x
rowsum=apply(P,1,sum) #transition matrix P (rows sum up to 1)
Pt=t(P) 
eigen(Pt) #take the vector corresponding to eigenvalue of 1 - this is pi transpose
pi=1/sum(-eigen(Pt)$vectors[,1])*(-eigen(Pt)$vectors[,1]) #rescale so that it sums up to 1

#Limiting/stationary distribution of Markov chain
matrix=P
for (i in 1:1000) {
  matrix=matrix%*%P
}

#Finding stationary distribution for a simulated second order Markov chain 
x=matrix(runif(64),nrow=16)
rowsum=apply(x,1,sum)
a=diag(1/rowsum,16,16)
P=a%*%x
values=as.vector(t(P))
values1=c(values[1:4],rep(0,12),rep(0,4),values[5:8],rep(0,8),rep(0,8),values[9:12],rep(0,4),rep(0,12),values[13:16],
          values[17:20],rep(0,12),rep(0,4),values[21:24],rep(0,8),rep(0,8),values[25:28],rep(0,4),rep(0,12),values[29:32],
          values[33:36],rep(0,12),rep(0,4),values[37:40],rep(0,8),rep(0,8),values[41:44],rep(0,4),rep(0,12),values[45:48],
          values[49:52],rep(0,12),rep(0,4),values[53:56],rep(0,8),rep(0,8),values[57:60],rep(0,4),rep(0,12),values[61:64])
values1=matrix(values1,byrow=TRUE,nrow=16,ncol=16)
rowsum=apply(values1,1,sum) #transition matrix P (rows sum up to 1)
Pt=t(values1) 
eigen(Pt) #take the vector corresponding to eigenvalue of 1 - this is pi transpose
pi=1/sum(-eigen(Pt)$vectors[,1])*(-eigen(Pt)$vectors[,1]) #rescale so that it sums up to 1
pi=matrix(pi,nrow=4,byrow=TRUE)
rowSums(pi)
colSums(pi)

#Simulation for test for time homogeneity

set.seed(1) #for reproducibility
results=numeric(100)
simulation=function(){
  sample=numeric(150000)
  sample[1]=sample(1:4,1)
  for (i in 2:10000) {
    current_state <- sample[i - 1]
    next_state <- sample(1:4, 1, prob = P[current_state, ])
    sample[i] <- next_state}
  counts <- matrix(0, nrow = 4, ncol = 4)
  for (i in 1:(length(sample) - 1)) {
    first_num <- sample[i]
    second_num <- sample[i + 1]
    counts[first_num, second_num] <- counts[first_num, second_num] + 1
  }
  nij=numeric(16)
  index=1
  for (i in 1:nrow(counts)) {
    for (j in 1:ncol(counts)) {
      nij[index]=counts[i,j]
      index=index+1
    }
  }
  nij=rep(nij,2)
  ni.=rep(rep(rowSums(counts),rep(4,4)),2)
  #take S=2 
  sample_1=sample[1:5000]
  counts_1 <- matrix(0, nrow = 4, ncol = 4)
  for (i in 1:length(sample_1)-1) {
    first_num <- sample_1[i]
    second_num <- sample_1[i + 1]
    counts_1[first_num, second_num] <- counts_1[first_num, second_num] + 1
  }
  nij_1=numeric(16)
  index=1
  for (i in 1:nrow(counts_1)) {
    for (j in 1:ncol(counts_1)) {
      nij_1[index]=counts_1[i,j]
      index=index+1
    }
  }
  ni._1=rep(rowSums(counts_1),rep(4,4))
  sample_2=sample[5000:10000]
  counts_2 <- matrix(0, nrow = 4, ncol = 4)
  for (i in 1:length(sample_2)-1) {
    first_num <- sample_2[i]
    second_num <- sample_2[i + 1]
    counts_2[first_num, second_num] <- counts_2[first_num, second_num] + 1
  }
  nij_2=numeric(16)
  index=1
  for (i in 1:nrow(counts_2)) {
    for (j in 1:ncol(counts_2)) {
      nij_2[index]=counts_2[i,j]
      index=index+1
    }
  }
  ni._2=rep(rowSums(counts_2),rep(4,4))
  nij_s=c(nij_1,nij_2)
  ni._s=c(ni._1,ni._2)
  test=2*sum((nij_s)*(log(nij_s/ni._s)-log(nij/ni.)))
  return(test)
}
for (i in 1:length(results)) {
  results[i] <- simulation()
}
ks.test(results,pchisq,df=12) 
  
#Simulation for LR/Pearson's Chi-Squared Test for independence (Order 0) vs Order 1
#Asymptotic chi-square dist. holds for all values of n 
set.seed(1)
x=matrix(runif(4),nrow=1) 
rowsum=apply(x,1,sum)
P=x/rowsum #probability distribution for an order 0 Markov chain 
results=numeric(1000) #generate 1000 test statistics
simulation=function(){
  sample=sample(1:4, size=31, replace = TRUE,prob=P) #size: length of each Markov chain generated
  all_combinations <- generate_combinations(4, 2)
  window_size <- 2
  combinations <- embed(sample, window_size)[, window_size:1]
  words <- apply(combinations, 1, paste, collapse = "")
  two_letter_count <- table(factor(words, levels = all_combinations))
  nab <- as.vector(two_letter_count)
  two_words_matrix <- matrix(nab, nrow = 4, ncol = 4,byrow=TRUE)
  na.=rep(rowSums(two_words_matrix),rep(4,4))
  counts=numeric(4)
  for (i in 2:(length(sample))) {
    num <- sample[i]
    counts[num] <- counts[num] + 1
  }
  n.=rep(sum(counts),16)
  nb=rep(counts,4)
  expected=(na.*nb)/(n.)
  valid = !is.na(expected) & expected > 0
  test=sum(((nab[valid]-expected[valid])^2)/expected[valid])
  return(test)
}
for (i in 1:length(results)){
  results[i]<-simulation()
}
ks.test(results,pchisq,df=9) 
par(mfrow=c(1,1))
#QQ-Plot
qqplot(qchisq(ppoints(1000), df = 9), results,
       main = "Q-Q Plot for Test of Independence Against First Order MC",
       xlab = "Theoretical Quantiles (Chi-Squared with D.F. of 9)",
       ylab = "Sample Quantiles",
       pch = 19, col = "blue")
abline(0, 1, col = "red", lwd = 2)

#Simulation for LR/Pearson's Chi-Squared Test for Order 1 vs Order 2
set.seed(1)
x=matrix(runif(16),nrow=4)
rowsum=apply(x,1,sum)
a=diag(1/rowsum,4,4)
P=a%*%x #transition matrix for a simulated first-order MC
results=numeric(1000)
generate_first_order_chain <- function(P, n) {
  chain <- numeric(n)
  chain[1] <- sample(1:4, 1) 
  for (i in 2:n) {
    current_state <- chain[i - 1]
    next_state <- sample(1:4, 1, prob = P[current_state, ])
    chain[i] <- next_state
  }
  return(chain)
}
simulation=function(){
  sample=generate_first_order_chain(P,339)  #n: length of each MC
  all_combinations <- generate_combinations(4, 3)
  window_size <- 3
  combinations <- embed(sample, window_size)[, window_size:1]
  words <- apply(combinations, 1, paste, collapse = "")
  three_letter_count <- table(factor(words, levels = all_combinations))
  nabc <- as.vector(three_letter_count)
  three_words_matrix <- matrix(nabc, nrow = 4^2, ncol = 4,byrow=TRUE)
  nab.=rep(rowSums(three_words_matrix),rep(4,4^2))
  
  window_size2 <- 2
  combinations2 <- embed(sample[-1], window_size2)[, window_size2:1]
  words2 <- apply(combinations2, 1, paste, collapse = "")
  two_letter = table(words2)
  if (length(two_letter) < 16) {
    two_letter <- c(two_letter, rep(0, 16 - length(two_letter))) # Ensure it has 16 counts
  }
  two_words_matrix=matrix(two_letter,nrow=4,ncol=4,byrow=TRUE)
  nb.=rep(rep(rowSums(two_words_matrix),rep(4,4)),4)
  nbc=rep(as.vector(two_letter),4)
  expected=(nab.*nbc)/(nb.)
  valid = !is.na(expected) & expected > 0
  test=sum(((nabc[valid]-expected[valid])^2)/expected[valid])
  return(test)
}
for (i in 1:length(results)) {
  results[i] <- simulation()
}
ks.test(results,pchisq,df=36)
results1vs2<-results #used as empirical distribution
#QQ-Plot
qqplot(qchisq(ppoints(1000), df = 36), results,
       main = "Q-Q Plot for Test of First Order MC Against Second Order MC",
       xlab = "Theoretical Quantiles (Chi-Squared with D.F. of 36)",
       ylab = "Sample Quantiles",
       pch = 19, col = "blue")
abline(0, 1, col = "red", lwd = 2)

#Simulation for LR/Pearson's Chi-Squared Test for Order 2 vs Order 3
set.seed(1)
x=matrix(runif((4^2)*4),ncol=4)
rowsum=apply(x,1,sum)
a=diag(1/rowsum,4^2,4^2)
P=a%*%x #transition matrix for 2nd-Order MC
results=numeric(1000)
generate_second_order_chain <- function(P, n) {
  chain <- numeric(n)
  chain[1:2] <- sample(1:4, 2, replace = TRUE)
  for (i in 3:n) {
    previous_state1 <- chain[i - 2]
    previous_state2 <- chain[i - 1]
    current_index <- (previous_state1 - 1) * 4 + previous_state2
    next_state <- sample(1:4, 1, prob = P[current_index, ])
    chain[i] <- next_state
  }
  return(chain)
}
simulation <- function() {
  sample=generate_second_order_chain(P,4500)  
  all_combinations <- generate_combinations(4, 4)
  window_size <- 4
  combinations <- embed(sample, window_size)[, window_size:1]
  words <- apply(combinations, 1, paste, collapse = "")
  four_letter_count <- table(factor(words, levels = all_combinations))
  nabcd <- as.vector(four_letter_count)
  if (length(nabcd) < 64) {
    nabcd <- c(nabcd, rep(0, 64 - length(nabcd)))
  }
  four_words_matrix <- matrix(nabcd, nrow = 4^3, ncol = 4,byrow=TRUE)
  nabc.=rep(rowSums(four_words_matrix),rep(4,4^3))
  
  window_size2 <- 3
  combinations2 <- embed(sample[-1], window_size2)[, window_size2:1]
  words2 <- apply(combinations2, 1, paste, collapse = "")
  three_letter = table(words2)
  if (length(three_letter) < 64) {
    three_letter <- c(three_letter, rep(0, 64 - length(three_letter)))
  }
  three_words_matrix=matrix(three_letter,nrow=16,ncol=4,byrow=TRUE)
  nbc.=rep(rep(rowSums(three_words_matrix),rep(4,16)),4)
  nbcd=rep(as.vector(three_letter),4)
  expected=(nabc.*nbcd)/(nbc.)
  valid=expected>0 & !is.na(expected)
  test=sum(((nabcd[valid]-expected[valid])^2)/expected[valid])
  return(test)
}
for (i in 1:length(results)) {
  results[i] <- simulation()
}
ks.test(results,pchisq,df=144)
results2vs3<-results #used as empirical distribution
#QQ-Plot
qqplot(qchisq(ppoints(1000), df = 144), results,
       main = "Q-Q Plot for Test of Second Order MC Against Third Order MC",
       xlab = "Theoretical Quantiles (Chi-Squared with D.F. of 144)",
       ylab = "Sample Quantiles",
       pch = 19, col = "blue")
abline(0, 1, col = "red", lwd = 2)

#Simulation for LR/Pearson's Chi-Squared Test for Order 3 vs Order 4
set.seed(1)
x=matrix(runif((4^3)*4),ncol=4)
rowsum=apply(x,1,sum)
a=diag(1/rowsum,4^3,4^3)
P=a%*%x #transition matrix for 3rd-Order MC
results=numeric(1000)
generate_third_order_chain <- function(P, n) {
  initial_sequence <- sample(1:4, 3, replace = TRUE)
  chain <- initial_sequence
  current_state <- initial_sequence[1] + (initial_sequence[2] - 1) * 4 + (initial_sequence[3] - 1) * 16
  for (i in 4:n) { #
    next_nucleotide_prob <- P[current_state, ]
    next_nucleotide <- sample(1:4, 1, prob = next_nucleotide_prob)  
    chain <- c(chain, next_nucleotide)  
    current_state <- (current_state %% (4^2)) * 4 + next_nucleotide  
  }
  return(chain)
}
simulation <- function() {
  sample=generate_third_order_chain(P,6000) 
  all_combinations <- generate_combinations(4, 5)
  window_size <- 5
  combinations <- embed(sample, window_size)[, window_size:1]
  words <- apply(combinations, 1, paste, collapse = "")
  five_letter_count <- table(factor(words, levels = all_combinations))
  nabcde <- as.vector(five_letter_count)
  if (length(nabcde) < 1024) { 
    nabcde <- c(nabcde, rep(0, 1024 - length(nabcde)))
  }
  five_words_matrix <- matrix(nabcde, nrow = 4^4, ncol = 4,byrow=TRUE)
  nabcd.=rep(rowSums(five_words_matrix),rep(4,4^4))
  
  all_combinations <- generate_combinations(4, 4)
  window_size2 <- 4
  combinations2 <- embed(sample[-1], window_size2)[, window_size2:1]
  words2 <- apply(combinations2, 1, paste, collapse = "")
  four_letter_count <- table(factor(words2, levels = all_combinations))
  nbcde <- as.vector(four_letter_count)
  if (length(nbcde) < 256) { 
    nbcde <- c(nbcde, rep(0, 256 - length(nbcde)))
  }
  four_words_matrix=matrix(four_letter_count,nrow=4^3,ncol=4,byrow=TRUE)
  nbcd.=rep(rep(rowSums(four_words_matrix),rep(4,4^3)),4)
  nbcde=rep(as.vector(four_letter_count),4)
  expected=(nabcd.*nbcde)/(nbcd.)
  valid_indices <- !is.na(expected) & expected > 0
  nabcde_valid <- nabcde[valid_indices]
  expected_valid <- expected[valid_indices]
  test=sum(((nabcde_valid-expected_valid)^2)/expected_valid)
  return(test)
}
for (i in 1:length(results)) {
  results[i] <- simulation()
}
ks.test(results,pchisq,df=576)
results3vs4<-results #used as empirical distribution
#QQ-Plot
qqplot(qchisq(ppoints(1000), df = 576), results,
       main = "Q-Q Plot for Test of Third Order MC Against Fourth Order MC",
       xlab = "Theoretical Quantiles (Chi-Squared with D.F. of 576)",
       ylab = "Sample Quantiles",
       pch = 19, col = "blue")
abline(0, 1, col = "red", lwd = 2)

#Simulation for LR/Pearson's Chi-Squared Test for Order 4 vs Order 5
#The max length of gene sequences which is 32,759 still does not asymptotically converge to a chi-squared dist. of the required D.F.
set.seed(1)
x=matrix(runif((4^4)*4),ncol=4)
rowsum=apply(x,1,sum)
a=diag(1/rowsum,4^4,4^4)
P=a%*%x #transition matrix for 4th-Order MC
results=numeric(1000)
generate_fourth_order_chain <- function(P, n) {
  initial_sequence <- sample(1:4, 4, replace = TRUE)
  chain <- initial_sequence
  current_state <- initial_sequence[1] + (initial_sequence[2] - 1) * 4 +
    (initial_sequence[3] - 1) * 16 + (initial_sequence[4] - 1) * 64
  for (i in 5:n) {
    next_nucleotide_prob <- P[current_state, ]
    next_nucleotide <- sample(1:4, 1, prob = next_nucleotide_prob)
    chain <- c(chain, next_nucleotide)
    current_state <- (current_state %% (4^3)) * 4 + next_nucleotide
  }
  return(chain)
}
simulation <- function() {
  sample <-generate_fourth_order_chain(P,32759)
  all_combinations <- generate_combinations(4, 6)
  window_size <- 6
  combinations <- embed(sample, window_size)[, window_size:1]
  words <- apply(combinations, 1, paste, collapse = "")
  six_letter_count <- table(factor(words, levels = all_combinations))
  nabcdef <- as.vector(six_letter_count)
  if (length(nabcdef) < 4096) {
    nabcdef <- c(nabcdef, rep(0, 4096 - length(nabcdef)))
  }
  six_words_matrix <- matrix(nabcdef, nrow = 4^5, ncol = 4,byrow=TRUE)
  nabcde.=rep(rowSums(six_words_matrix),rep(4,4^5))
  
  all_combinations <- generate_combinations(4, 5)
  window_size2 <- 5
  combinations2 <- embed(sample[-1], window_size2)[, window_size2:1]
  words2 <- apply(combinations2, 1, paste, collapse = "")
  five_letter_count <- table(factor(words2, levels = all_combinations))
  nbcdef <- as.vector(five_letter_count)
  if (length(nbcdef) < 1024) { 
    nbcdef <- c(nbcdef, rep(0, 1024 - length(nbcdef)))
  }
  five_words_matrix=matrix(five_letter_count,nrow=4^4,ncol=4,byrow=TRUE)
  nbcde.=rep(rep(rowSums(five_words_matrix),rep(4,4^4)),4)
  nbcdef=rep(as.vector(five_letter_count),4)
  expected=(nabcde.*nbcdef)/(nbcde.)
  valid_indices <- expected > 0 
  nabcdef_valid <- nabcdef[valid_indices]
  expected_valid <- expected[valid_indices]
  test=sum(((nabcdef_valid-expected_valid)^2)/expected_valid)
  return(test)
}
for (i in 1:length(results)) {
  results[i] <- simulation()
}
ks.test(results,pchisq,df=2304)
results4vs5_32759=results #used as empirical dist. for gene sequences
#results4vs5_24562=results #used as empirical dist. for coding DNA sequences
#QQ-Plot
qqplot(qchisq(ppoints(1000), df = 2304), results,
       main = "Q-Q Plot for Test of Fourth Order MC Against Fifth Order MC",
       xlab = "Theoretical Quantiles (Chi-Squared with D.F. of 2304)",
       ylab = "Sample Quantiles",
       pch = 19, col = "blue")
abline(0, 1, col = "red", lwd = 2)
