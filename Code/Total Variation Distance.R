#TVD between Order 0 and Order 1
var_dist_indepvs1 <- function(seq){
  sample <- unlist(strsplit(seq, ""))
  sample=gsub("A",1,sample)
  sample=gsub("C",2,sample)
  sample=gsub("G",3,sample)
  sample=gsub("T",4,sample)
  sample=as.numeric(sample)
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
  indep=numeric(4)
  for (i in 1:length(sample)) {
    num=sample[i]
    indep[num]=indep[num]+1
  }
  indep=indep/length(sample)
  T1=matrix(nab/na.,byrow=T,nrow=4)
  T1[is.na(T1)] <- 0
  abs_dist=abs(sweep(T1,2,indep,FUN="-")) #calculate the absolute differences between p_ab and p_b for all a,b=1,2,3 and 4 
  var_dist=1/2*1/4*sum(abs_dist) #every row of the transition matrix T1 is a probability distribution, so there are 4 probability distributions
  return(var_dist)
}
genome_genes_sample$var_dist_indepvs1 <- sapply(genome_genes_sample$Sequence, var_dist_indepvs1) #for the sample of gene sequences

#TVD between Order 1 and Order 2
var_dist_1vs2 <- function(seq){
  sample <- unlist(strsplit(seq, ""))
  sample=gsub("A",1,sample)
  sample=gsub("C",2,sample)
  sample=gsub("G",3,sample)
  sample=gsub("T",4,sample)
  sample=as.numeric(sample)
  all_combinations <- generate_combinations(4, 2)
  window_size <- 2
  combinations <- embed(sample, window_size)[, window_size:1]
  words <- apply(combinations, 1, paste, collapse = "")
  two_letter_count <- table(factor(words, levels = all_combinations))
  nab <- as.vector(two_letter_count)
  two_words_matrix <- matrix(nab, nrow = 4, ncol = 4,byrow=TRUE)
  na.=rep(rowSums(two_words_matrix),rep(4,4))
  T1=matrix(nab/na.,byrow=T,nrow=4)
  T1[is.na(T1)] <- 0 #handles cases where there are no counts of letters starting with a (leading to a division by 0, producing NA values)
  all_combinations <- generate_combinations(4, 3)
  window_size <- 3
  combinations <- embed(sample, window_size)[, window_size:1]
  words <- apply(combinations, 1, paste, collapse = "")
  three_letter_count <- table(factor(words, levels = all_combinations))
  nabc <- as.vector(three_letter_count)
  if (length(nabc) < 64) {
    nabc <- c(nabc, rep(0, 64 - length(nabc)))
  }
  three_words_matrix <- matrix(nabc, nrow = 4^2, ncol = 4,byrow=TRUE)
  nab.=rep(rowSums(three_words_matrix),rep(4,4^2))
  T2=matrix(nabc/nab.,byrow=T,nrow=16)
  T2[is.na(T2)] <- 0
  rows_list2 <- lapply(1:4, function(i) seq(i, 16, by=4)) 
  abs_dist_list2 <- lapply(1:4, function(i) abs(sweep(T2[rows_list2[[i]], ], 2, T1[[i]]))) 
  #the above code calculates the absolute differences between (p_1bc, p_2bc, p_3bc, p_4bc) and p_bc for all b,c=1,2,3 and 4 
  abs_dist_combined2 <- do.call(rbind, abs_dist_list2)
  var_dist <- (1/2)*(1/16)*sum(abs_dist_combined2,na.rm=TRUE) #T2 has 16 rows, hence 16 probability distributions
  return(var_dist)
} 
genome_genes_sample$var_dist_1vs2 <- sapply(genome_genes_sample$Sequence, var_dist_1vs2)

#TVD between Order 2 and Order 3
var_dist_2vs3 <- function(seq){
  sample <- unlist(strsplit(seq, ""))
  sample=gsub("A",1,sample)
  sample=gsub("C",2,sample)
  sample=gsub("G",3,sample)
  sample=gsub("T",4,sample)
  sample=as.numeric(sample)
  all_combinations <- generate_combinations(4, 3)
  window_size <- 3
  combinations <- embed(sample, window_size)[, window_size:1]
  words <- apply(combinations, 1, paste, collapse = "")
  three_letter_count <- table(factor(words, levels = all_combinations))
  nabc <- as.vector(three_letter_count)
  three_words_matrix <- matrix(nabc, nrow = 4^2, ncol = 4,byrow=TRUE)
  nab.=rep(rowSums(three_words_matrix),rep(4,4^2))
  T2=matrix(nabc/nab.,byrow=T,nrow=16)
  T2[is.na(T2)] <- 0
  all_combinations <- generate_combinations(4, 4)
  window_size <- 4
  combinations <- embed(sample, window_size)[, window_size:1]
  words <- apply(combinations, 1, paste, collapse = "")
  four_letter_count <- table(factor(words, levels = all_combinations))
  nabcd <- as.vector(four_letter_count)
  if (length(nabcd) < 256) {
    nabcd <- c(nabcd, rep(0, 256 - length(nabcd)))
  }
  four_words_matrix <- matrix(nabcd, nrow = 4^3, ncol = 4,byrow=TRUE)
  nabc.=rep(rowSums(four_words_matrix),rep(4,4^3))
  T3=matrix(nabcd/nabc.,byrow=T,nrow=64)
  T3[is.na(T3)] <- 0
  rows_list3 <- lapply(1:16, function(i) seq(i, 64, by=16))
  abs_dist_list3 <- lapply(1:16, function(i) abs(sweep(T3[rows_list3[[i]], ], 2, T2[[i]])))
  abs_dist_combined3 <- do.call(rbind, abs_dist_list3)
  var_dist <- (1/2)*(1/64)*sum(abs_dist_combined3,na.rm=TRUE)
  return(var_dist)
} 
genome_genes_sample$var_dist_2vs3 <- sapply(genome_genes_sample$Sequence, var_dist_2vs3)

#TVD between Order 3 and Order 4
var_dist_3vs4<- function(seq){
  sample <- unlist(strsplit(seq, ""))
  sample=gsub("A",1,sample)
  sample=gsub("C",2,sample)
  sample=gsub("G",3,sample)
  sample=gsub("T",4,sample)
  sample=as.numeric(sample)
  all_combinations <- generate_combinations(4, 4)
  window_size <- 4
  combinations <- embed(sample, window_size)[, window_size:1]
  words <- apply(combinations, 1, paste, collapse = "")
  four_letter_count <- table(factor(words, levels = all_combinations))
  nabcd <- as.vector(four_letter_count)
  if (length(nabcd) < 256) {
    nabcd <- c(nabcd, rep(0, 256 - length(nabcd)))
  }
  four_words_matrix <- matrix(nabcd, nrow = 4^3, ncol = 4,byrow=TRUE)
  nabc.=rep(rowSums(four_words_matrix),rep(4,4^3))
  T3=matrix(nabcd/nabc.,byrow=T,nrow=64)
  T3[is.na(T3)] <- 0
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
  T4=matrix(nabcde/nabcd.,byrow=T,nrow=256)
  T4[is.na(T4)] <- 0
  rows_list4 <- lapply(1:64, function(i) seq(i, 256, by=64))
  abs_dist_list4 <- lapply(1:64, function(i) abs(sweep(T4[rows_list4[[i]], ], 2, T3[[i]])))
  abs_dist_combined4 <- do.call(rbind, abs_dist_list4)
  var_dist <- (1/2)*(1/256)*sum(abs_dist_combined4,na.rm=TRUE)
  return(var_dist)
}
genome_genes_sample$var_dist_3vs4 <- sapply(genome_genes_sample$Sequence, var_dist_3vs4)

#TVD between Order 4 and Order 5
var_dist_4vs5<- function(seq){
  sample <- unlist(strsplit(seq, ""))
  sample=gsub("A",1,sample)
  sample=gsub("C",2,sample)
  sample=gsub("G",3,sample)
  sample=gsub("T",4,sample)
  sample=as.numeric(sample)
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
  T4=matrix(nabcde/nabcd.,byrow=T,nrow=256)
  T4[is.na(T4)] <- 0
  all_combinations <- generate_combinations(4, 6)
  window_size <- 6
  combinations <- embed(sample, window_size)[, window_size:1]
  words <- apply(combinations, 1, paste, collapse = "")
  six_letter_count <- table(factor(words, levels = all_combinations))
  nabcdef <- as.vector(six_letter_count)
  if (length(nabcdef) < 4096) { 
    nabcdef <- c(nabcdef, rep(0, 1024 - length(nabcdef)))
  }
  six_words_matrix <- matrix(nabcdef, nrow = 4^5, ncol = 4,byrow=TRUE)
  nabcde.=rep(rowSums(six_words_matrix),rep(4,4^5))
  T5=matrix(nabcdef/nabcde.,byrow=T,nrow=1024)
  T5[is.na(T5)] <- 0
  rows_list5 <- lapply(1:256, function(i) seq(i, 1024, by=256))
  abs_dist_list5 <- lapply(1:256, function(i) abs(sweep(T5[rows_list5[[i]], ], 2, T4[[i]])))
  abs_dist_combined5 <- do.call(rbind, abs_dist_list5)
  var_dist <- (1/2)*(1/1024)*sum(abs_dist_combined5,na.rm=TRUE)
  return(var_dist)
}
genome_genes_sample$var_dist_4vs5 <- sapply(genome_genes_sample$Sequence, var_dist_4vs5)

#Line graph for Mean TVD 
par(mfrow=c(1,1))
var_dist_data=data.frame(genome_genes_sample$var_dist_indepvs1, genome_genes_sample$var_dist_1vs2, genome_genes_sample$var_dist_2vs3, genome_genes_sample$var_dist_3vs4, genome_genes_sample$var_dist_4vs5)
colnames(var_dist_data)<-c("Order 0 vs Order 1","Order 1 vs Order 2", "Order 2 vs Order 3", "Order 3 vs Order 4", "Order 4 vs Order 5" )
means <- colMeans(var_dist_data, na.rm = TRUE)
std_errors <- apply(var_dist_data, 2, function(x) sd(x, na.rm = TRUE) / sqrt(length(na.omit(x))))
lower_bound <- means - qnorm(0.975) * std_errors
upper_bound <- means + qnorm(0.975) * std_errors
plot(1:5, means, type = "o", pch = 19, col = "#C2185B", ylim = range(c(lower_bound, upper_bound)),
     xlab = "Orders of Markov Chain", ylab = "Mean TVD", main = "Mean TVD Between Consecutive Markov Chain Orders", xaxt = "n")
x_labels <- c("Order 0 vs 1", "Order 1 vs 2", "Order 2 vs 3", "Order 3 vs 4", "Order 4 vs 5")
axis(1, at = 1:5, labels = x_labels)
arrows(1:5, lower_bound, 1:5, upper_bound, angle = 90, code = 3, length = 0.1, col = "purple")