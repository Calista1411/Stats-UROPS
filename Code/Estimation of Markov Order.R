#Creating a dataframe called genome_cds containing all the coding DNA sequences in the genome and their corresponding lengths
file_list <- list.files(path = "~/Downloads/urops/csv files of all cds", pattern = "*.csv", full.names = TRUE)
genome_cds <- do.call(rbind, lapply(file_list, read.csv))
genome_cds <- genome_cds[order(as.numeric(genome_cds$Chr)), ]
genome_cds <- genome_cds[genome_cds$Sequence != "" & !is.na(genome_cds$Sequence), ] #remove rows with no gene name
genome_cds$length <- nchar(genome_cds$Sequence)
genome_cds$length <- sapply(genome_cds$Sequence, function(cds) {
  cds <- gsub("[^ACGT]", "", cds)
  sample <- unlist(strsplit(cds, ""))
  sample <- gsub("A", 1, sample)
  sample <- gsub("C", 2, sample)
  sample <- gsub("G", 3, sample)
  sample <- gsub("T", 4, sample)
  sample <- as.numeric(sample)
  length(sample)
})
genome_cds<- genome_cds[order(as.numeric(genome_cds$length)), ]

#Creating a dataframe called genome_genes containing all the genes sequences in the genome and their corresponding lengths
file_list <- list.files(path = "~/Downloads/urops/csv files of all genes", pattern = "*.csv", full.names = TRUE)
genome_genes <- do.call(rbind, lapply(file_list, read.csv))
genome_genes <- genome_genes[order(as.numeric(genome_genes$Chr)), ]
genome_genes <- genome_genes[genome_genes$Sequence != "" & !is.na(genome_genes$Sequence), ] #remove rows with no gene name
genome_genes$length <- nchar(genome_genes$Sequence)
genome_genes$length <- sapply(genome_genes$Sequence, function(cds) {
  cds <- gsub("[^ACGT]", "", cds)
  sample <- unlist(strsplit(cds, ""))
  sample <- gsub("A", 1, sample)
  sample <- gsub("C", 2, sample)
  sample <- gsub("G", 3, sample)
  sample <- gsub("T", 4, sample)
  sample <- as.numeric(sample)
  length(sample)
})
genome_genes<- genome_genes[order(as.numeric(genome_genes$length)), ]

#Matching the dataset of gene sequences to the dataset of coding DNA sequences to ensure that the dataset consists of only matched pairs of gene sequences and CDS corresponding to the same gene
library(dplyr)
genome_cds <- genome_cds %>% rename(Gene.Name = "CDS.Name")
genome_genes_matched <- genome_genes %>%
  semi_join(genome_cds, by = c("Gene.Name", "Chr")) #match by chromosome number and gene name

#Taking a random sample of of 1,000 matched pairs of gene and coding DNA sequences
set.seed(1) #for reproducibility
genome_cds_sample<-genome_cds[sample(nrow(genome_cds), 1000), ] 
genome_genes_sample<-genome_genes %>%
  semi_join(genome_cds_sample, by = c("Gene.Name", "Chr")) 

#Function to generate combinations (n is the state space, length is the number of letters in each word combination)
generate_combinations <- function(n, length) {
  seq <- 1:n
  result <- list()
  generate <- function(current_combination, remaining_length) {
    if (remaining_length == 0) {
      result <<- c(result, list(paste(current_combination, collapse = "")))
    } else {
      for (i in seq) {
        generate(c(current_combination, i), remaining_length - 1)
      }
    }
  }
  generate(integer(0), length)
  return(unlist(result))
}

#Goodness-of-fit Tests 

#LR Test or Pearson's Chi-Squared Test for Order 0 vs Order 1
first_order <- function(seq) {
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
  test=sum(((nab[valid]-expected[valid])^2)/expected[valid])
  #test=2*sum(nab*(log(nab/na.)-log(nb/n.)))
  pvalue=pchisq(test,df=9,lower.tail=FALSE)
  return(pvalue)
}
genome_cds_sample$pvalue_indep_0.1 <- sapply(genome_cds_sample$Sequence, first_order)
genome_genes_sample$pvalue_indep_0.1 <- sapply(genome_genes_sample$Sequence, first_order)

#LR Test or Pearson's Chi-Squared Test for Order 1 vs Order 2
second_order <- function(seq,length) {
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
  if (length(nabc) < 64) {
    nabcd <- c(nabcd, rep(0, 64 - length(nabcd)))
  }
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
  if (length > 350) {
    pvalue <- pchisq(test, df = 36, lower.tail = FALSE) #Use empirical distribution obtained via simulation
  } else {
    pvalue <- mean(abs(results1vs2) >= abs(test))
  }
  return(pvalue)
}
#for alpha=0.1
genome_cds_sample$pvalue_1vs2_0.1 <- NA
indices <- which(genome_cds_sample$pvalue_indep_0.1 < 0.1)
genome_cds_sample$pvalue_1vs2_0.1[indices] <- mapply(second_order, 
                                                             genome_cds_sample$Sequence[indices],
                                                             genome_cds_sample$length[indices])
#changing the value of alpha
genome_cds_sample$pvalue_1vs2_0.05<-genome_cds_sample$pvalue_1vs2_0.1 #for alpha=0.05
indices <- which(genome_cds_sample$pvalue_indep_0.1 > 0.05)
genome_cds_sample$pvalue_1vs2_0.05[indices] <- NA
#for alpha=0.1
genome_genes_sample$pvalue_1vs2_0.1 <- NA
indices <- which(genome_genes_sample$pvalue_indep_0.1 < 0.1)
genome_genes_sample$pvalue_1vs2_0.1[indices] <- mapply(second_order, 
                                                     genome_genes_sample$Sequence[indices],
                                                     genome_genes_sample$length[indices])
#changing the value of alpha
genome_genes_sample$pvalue_1vs2_0.05<-genome_genes_sample$pvalue_1vs2_0.1 #for alpha=0.05
indices <- which(genome_genes_sample$pvalue_indep_0.1 > 0.05)
genome_genes_sample$pvalue_1vs2_0.05[indices] <- NA

#LR Test or Pearson's Chi-Squared Test for Order 2 vs Order 3
third_order <- function(seq,length) {
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
  valid = !is.na(expected) & expected > 0
  test=sum(((nabcd[valid]-expected[valid])^2)/expected[valid])
  if (length > 5000) {
    pvalue <- pchisq(test, df = 144, lower.tail = FALSE)
  } else {
    pvalue <- mean(abs(results2vs3) >= abs(test))
  }
  pvalue <- pchisq(test, df = 144, lower.tail = FALSE)
  return(pvalue)
}
#for alpha=0.1
genome_cds_sample$pvalue_2vs3_0.1 <- NA
indices <- which(genome_cds_sample$pvalue_1vs2_0.1 < 0.1)
genome_cds_sample$pvalue_2vs3_0.1[indices] <- mapply(third_order, 
                                                      genome_cds_sample$Sequence[indices],
                                                      genome_cds_sample$length[indices]) 
#changing alpha
genome_cds_sample$pvalue_2vs3_0.05<-genome_cds_sample$pvalue_2vs3_0.1 
indices <- which(genome_cds_sample$pvalue_1vs2_0.05 > 0.05 | is.na(genome_cds_sample$pvalue_1vs2_0.05))
genome_cds_sample$pvalue_2vs3_0.05[indices] <- NA
#for alpha=0.1
genome_genes_sample$pvalue_2vs3_0.1 <- NA
indices <- which(genome_genes_sample$pvalue_1vs2_0.1 < 0.1)
genome_genes_sample$pvalue_2vs3_0.1[indices] <- mapply(third_order, 
                                                     genome_genes_sample$Sequence[indices],
                                                     genome_genes_sample$length[indices])
#changing alpha
genome_genes_sample$pvalue_2vs3_0.05<-genome_genes_sample$pvalue_2vs3_0.1 
indices <- which(genome_genes_sample$pvalue_1vs2_0.05 > 0.05 | is.na(genome_genes_sample$pvalue_1vs2_0.05))
genome_genes_sample$pvalue_2vs3_0.05[indices] <- NA

#LR Test or Pearson's Chi-Squared Test for Order 3 vs Order 4
fourth_order <- function(seq,length) {
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
  if (length > 6000) {
   pvalue <- pchisq(test, df = 576, lower.tail = FALSE)
  } else {
   pvalue <- mean(abs(results3vs4) >= abs(test))
  }
  return(pvalue)
}
#for alpha=0.1
genome_cds_sample$pvalue_3vs4_0.1 <- NA
indices <- which(genome_cds_sample$pvalue_2vs3_0.1 < 0.1)
genome_cds_sample$pvalue_3vs4_0.1[indices] <- mapply(fourth_order, 
                                                      genome_cds_sample$Sequence[indices],
                                                      genome_cds_sample$length[indices])
#changing alpha
genome_cds_sample$pvalue_3vs4_0.05<-genome_cds_sample$pvalue_3vs4_0.1 
indices <- which(genome_cds_sample$pvalue_2vs3_0.05 > 0.05 | is.na(genome_cds_sample$pvalue_2vs3_0.05))
genome_cds_sample$pvalue_3vs4_0.05[indices] <- NA
#for alpha=0.1
genome_genes_sample$pvalue_3vs4_0.1 <- NA
indices <- which(genome_genes_sample$pvalue_2vs3_0.1 < 0.1)
genome_genes_sample$pvalue_3vs4_0.1[indices] <- mapply(fourth_order, 
                                                     genome_genes_sample$Sequence[indices],
                                                     genome_genes_sample$length[indices])

#changing alpha
genome_genes_sample$pvalue_3vs4_0.05<-genome_genes_sample$pvalue_3vs4_0.1 
indices <- which(genome_genes_sample$pvalue_2vs3_0.05 > 0.05 | is.na(genome_genes_sample$pvalue_2vs3_0.05))
genome_genes_sample$pvalue_3vs4_0.05[indices] <- NA

#LR Test or Pearson's Chi-Squared Test for Order 3 vs Order 4
fifth_order <- function(seq) {
  sample <- unlist(strsplit(seq, ""))
  sample=gsub("A",1,sample)
  sample=gsub("C",2,sample)
  sample=gsub("G",3,sample)
  sample=gsub("T",4,sample)
  sample=as.numeric(sample)
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
  valid_indices <- expected > 0 & !is.na(expected)
  nabcdef_valid <- nabcdef[valid_indices]
  expected_valid <- expected[valid_indices]
  test=sum(((nabcdef_valid-expected_valid)^2)/expected_valid)
  pvalue <- mean(abs(results24562) >= abs(test)) #use empirical distribution
  return(pvalue)
}
#for alpha=0.1
genome_cds_sample$pvalue_4vs5_0.1 <- NA
indices <- which(genome_cds_sample$pvalue_3vs4_0.1 < 0.1)
genome_cds_sample$pvalue_4vs5_0.1[indices] <- sapply(genome_cds_sample$Sequence[indices], fifth_order)
#changing alpha
genome_cds_sample$pvalue_4vs5_0.05<-genome_cds_sample$pvalue_4vs5_0.1 #for alpha=0.01
indices <- which(genome_cds_sample$pvalue_3vs4_0.05 > 0.05 | is.na(genome_cds_sample$pvalue_3vs4_0.05))
genome_cds_sample$pvalue_4vs5_0.05[indices] <- NA
#for alpha=0.1
genome_cds_sample$pvalue_4vs5_0.1 <- NA
indices <- which(genome_cds_sample$pvalue_3vs4_0.1 < 0.1)
genome_cds_sample$pvalue_4vs5_0.1[indices] <- sapply(genome_cds_sample$Sequence[indices], fifth_order)
#changing alpha
genome_genes_sample$pvalue_4vs5_0.05<-genome_genes_sample$pvalue_4vs5_0.1 #for alpha=0.01
indices <- which(genome_genes_sample$pvalue_3vs4_0.05 > 0.05 | is.na(genome_genes_sample$pvalue_3vs4_0.05))
genome_genes_sample$pvalue_4vs5_0.05[indices] <- NA

#Summary of Order of Coding and Gene Sequences Using Sequential LR Tests

#Number of sequences of order 5 and higher
length(which(genome_cds_sample$pvalue_4vs5_0.1<0.1))
length(which(genome_genes_sample$pvalue_4vs5_0.1<0.1))

#Number of sequences of order 0,1,2,3 and 4 respectively
order <- c("Order 0", "Order 1", "Order 2", "Order 3","Order 4")
number_of_genes <- c(sum(is.na(genome_genes_sample$pvalue_1vs2_0.1)),
                     sum(is.na(genome_genes_sample$pvalue_2vs3_0.1)) - sum(is.na(genome_genes_sample$pvalue_1vs2_0.1)),
                     sum(is.na(genome_genes_sample$pvalue_3vs4_0.1)) - sum(is.na(genome_genes_sample$pvalue_2vs3_0.1)),
                     sum(is.na(genome_genes_sample$pvalue_4vs5_0.1)) - sum(is.na(genome_genes_sample$pvalue_3vs4_0.1)),
                     length(which(genome_genes_sample$pvalue_4vs5_0.1>=0.1))
)
order <- c("Order 0", "Order 1", "Order 2", "Order 3","Order 4")
number_of_cds <- c(sum(is.na(genome_cds_sample$pvalue_1vs2_0.1)),
                     sum(is.na(genome_cds_sample$pvalue_2vs3_0.1)) - sum(is.na(genome_cds_sample$pvalue_1vs2_0.1)),
                     sum(is.na(genome_cds_sample$pvalue_3vs4_0.1)) - sum(is.na(genome_cds_sample$pvalue_2vs3_0.1)),
                     sum(is.na(genome_cds_sample$pvalue_4vs5_0.1)) - sum(is.na(genome_cds_sample$pvalue_3vs4_0.1)),
                     length(which(genome_cds_sample$pvalue_4vs5_0.1>=0.1))
)

#Comparative Bar Plot for Markov Chain Orders of Coding and Gene Sequences
par(mfrow=c(1,2))
order <- c(0:4,"5 or higher")
alpha_levels <- c(0.1, 0.05, 0.01, 0.001,0.00001)
frequencies <- matrix(c(26, 285, 205, 420, 57, 7,
                        33, 379, 212, 322, 49, 5,
                        65, 503, 183, 212, 32, 5, 
                        104, 593, 140, 134, 24, 5,
                        189, 643, 84, 61, 18, 5), 
                      nrow = 5, byrow = TRUE)
colors <- c("#A7C7E7", "#B5D3A7", "#F7D08A", "#FAB4B4", "#D3B5E5")
alphas <- c(0.1, 0.05, 0.01, 0.001, 0.0001)

barplot(
  frequencies, beside = TRUE, col = colors, 
  names.arg = order, 
  xlab = "Order of Markov Chain", ylab = "Frequency", 
  main = "Markov Chain Order of Coding DNA Sequences",
  ylim=c(0,1000)
)
legend("topright", legend = c("0.1", "0.05", "0.01", "0.001", "0.00001"),
       fill = colors, cex = 0.85,title="Significance Level") 

order <- c(0:4,"5 or higher")
alpha_levels <- c(0.1, 0.05, 0.01, 0.001,0.00001)
frequencies <- matrix(c(12, 46, 54, 77, 130, 681, 
                        12, 65, 63, 65, 130, 665,
                        17, 86, 65, 73, 117, 642,
                        20, 100, 88, 78, 91, 623,
                        25, 125, 115, 80, 46, 609), 
                      nrow = 5, byrow = TRUE)
barplot(
  frequencies, beside = TRUE, col = colors, 
  names.arg = order, 
  xlab = "Order of Markov Chain", ylab = "Frequency", 
  main = "Markov Chain Order of Gene Sequences",
  ylim=c(0,1000)
)
legend_x <- par("usr")[1] + 1  
legend_y <- par("usr")[4] - 0.5  
legend(legend_x, legend_y, legend = c("0.1", "0.05", "0.01", "0.001", "0.00001"),
       fill = colors, cex = 0.85, title = "Significance Level")

#Distribution of p-values for the LR test of Order 0 vs Order 1
par(mfrow=c(1,2))
hist(genome_cds_sample$pvalue_indep_0.1,col="#FAB4B4",main="Distribution of p-values of Coding DNA Sequences",prob=TRUE,xlab="p-value")
curve(dunif(x, min = 0, max = 1), 
      col = "#9A6FB0",
      lwd = 2, 
      add = TRUE)
ks.test(genome_cds_sample$pvalue_indep_0.1,punif,0,1)
hist(genome_genes_sample$pvalue_indep_0.1,col="#FAB4B4",main="Distribution of p-values of Gene Sequences",prob=TRUE,xlab="p-value",xlim=c(0,0.8))
curve(dunif(x, min = 0, max = 1), 
      col = "#9A6FB0",
      lwd = 2, 
      add = TRUE)
ks.test(genome_genes_sample$pvalue_indep_0.1,punif,0,1)

#Information Criteria 

#BIC  
BIC0<-function(seq){
  sample <- unlist(strsplit(seq, ""))
  sample=gsub("A",1,sample)
  sample=gsub("C",2,sample)
  sample=gsub("G",3,sample)
  sample=gsub("T",4,sample)
  sample=as.numeric(sample)
  counts <- table(factor(sample, levels = 1:4))
  proportions <- counts / length(sample)
  log_likelihood <- sum(counts * log(proportions))
  bic <- 3 * log(length(sample)) - 2 * log_likelihood
  return(bic)
}
genome_cds_sample$BIC0 <- sapply(genome_cds_sample$Sequence, BIC0)
genome_genes_sample$BIC0 <- sapply(genome_genes_sample$Sequence, BIC0)

BIC1 <- function(seq) {
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
  valid_indices=nab>0
  nab_valid=nab[valid_indices]
  na._valid=na.[valid_indices]
  bic=12*log(length(sample))-2*sum(nab_valid*log(nab_valid/na._valid))
  return(bic)
}
genome_cds_sample$BIC1 <- sapply(genome_cds_sample$Sequence, BIC1)
genome_genes_sample$BIC1 <- sapply(genome_genes_sample$Sequence, BIC1)

BIC2 <- function(seq) {
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
  valid=nabc>0
  bic=48*log(length(sample))-2*sum(nabc[valid]*log(nabc[valid]/nab.[valid]))
  return(bic)
}
genome_cds_sample$BIC2 <- sapply(genome_cds_sample$Sequence, BIC2)
genome_genes_sample$BIC2 <- sapply(genome_genes_sample$Sequence, BIC2)

BIC3 <- function(seq) {
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
  if (length(nabcd) < 64) {
    nabcd <- c(nabcd, rep(0, 64 - length(nabcd)))
  }
  four_words_matrix <- matrix(nabcd, nrow = 4^3, ncol = 4,byrow=TRUE)
  nabc.=rep(rowSums(four_words_matrix),rep(4,4^3))
  valid=nabcd>0 & nabc.>0
  bic=192*log(length(sample))-2*sum(nabcd[valid]*log(nabcd[valid]/nabc.[valid]))
  return(bic)
}
genome_cds_sample$BIC3 <- sapply(genome_cds_sample$Sequence, BIC3)
genome_genes_sample$BIC3 <- sapply(genome_genes_sample$Sequence, BIC3)

BIC4 <-function(seq) {
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
  valid=nabcde>0 & nabcd.>0
  bic=768*log(length(sample))-2*sum(nabcde[valid]*log(nabcde[valid]/nabcd.[valid]))
  return(bic)
}
genome_cds_sample$BIC4 <- sapply(genome_cds_sample$Sequence, BIC4)
genome_genes_sample$BIC4 <- sapply(genome_genes_sample$Sequence, BIC4)

BIC5<- function(seq) {
  sample <- unlist(strsplit(seq, ""))
  sample=gsub("A",1,sample)
  sample=gsub("C",2,sample)
  sample=gsub("G",3,sample)
  sample=gsub("T",4,sample)
  sample=as.numeric(sample)
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
  valid=nabcdef>0 & nabcde. >0
  bic=3072*log(length(sample))-2*sum(nabcdef[valid]*log(nabcdef[valid]/nabcde.[valid]))
  return(bic)
}
genome_cds_sample$BIC5 <- sapply(genome_cds_sample$Sequence, BIC5)
genome_genes_sample$BIC5 <- sapply(genome_genes_sample$Sequence, BIC5)

#Order with the minimum BIC for each sequence
genome_cds_sample$min_BIC_order <- apply(genome_cds_sample[, c("BIC0","BIC1","BIC2","BIC3","BIC4","BIC5")], 1, function(x) {
  order <- which.min(x) - 1  
  return(order)  
})
genome_genes_sample$min_BIC_order <- apply(genome_genes_sample[, c("BIC0","BIC1","BIC2","BIC3","BIC4","BIC5")], 1, function(x) {
  order <- which.min(x) - 1  
  return(order)  
})

#AIC
AIC0<-function(seq){
  sample <- unlist(strsplit(seq, ""))
  sample=gsub("A",1,sample)
  sample=gsub("C",2,sample)
  sample=gsub("G",3,sample)
  sample=gsub("T",4,sample)
  sample=as.numeric(sample)
  counts <- table(factor(sample, levels = 1:4))
  proportions <- counts / length(sample)
  log_likelihood <- sum(counts * log(proportions))
  aic <- 2*3 - 2 * log_likelihood
  return(aic)
}
genome_cds_sample$AIC0 <- sapply(genome_cds_sample$Sequence, AIC0)
genome_genes_sample$AIC0 <- sapply(genome_genes_sample$Sequence, AIC0)

AIC1 <- function(seq) {
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
  valid_indices=nab>0
  nab_valid=nab[valid_indices]
  na._valid=na.[valid_indices]
  aic=2*12-2*sum(nab_valid*log(nab_valid/na._valid))
  return(aic)
}
genome_cds_sample$AIC1 <- sapply(genome_cds_sample$Sequence, AIC1)
genome_genes_sample$AIC1 <- sapply(genome_genes_sample$Sequence, AIC1)

AIC2 <- function(seq) {
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
  valid=nabc>0
  aic=2*48-2*sum(nabc[valid]*log(nabc[valid]/nab.[valid]))
  return(aic)
}
genome_cds_sample$AIC2 <- sapply(genome_cds_sample$Sequence, AIC2)
genome_genes_sample$AIC2 <- sapply(genome_genes_sample$Sequence, AIC2)

AIC3 <- function(seq) {
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
  if (length(nabcd) < 64) {
    nabcd <- c(nabcd, rep(0, 64 - length(nabcd)))
  }
  four_words_matrix <- matrix(nabcd, nrow = 4^3, ncol = 4,byrow=TRUE)
  nabc.=rep(rowSums(four_words_matrix),rep(4,4^3))
  valid=nabcd>0 & nabc.>0
  aic=2*192-2*sum(nabcd[valid]*log(nabcd[valid]/nabc.[valid]))
  return(aic)
}
genome_cds_sample$AIC3 <- sapply(genome_cds_sample$Sequence, AIC3)
genome_genes_sample$AIC3 <- sapply(genome_genes_sample$Sequence, AIC3)

AIC4 <-function(seq) {
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
  valid=nabcde>0 & nabcd.>0
  aic=2*768-2*sum(nabcde[valid]*log(nabcde[valid]/nabcd.[valid]))
  return(aic)
}
genome_cds_sample$AIC4 <- sapply(genome_cds_sample$Sequence, AIC4)
genome_genes_sample$AIC4 <- sapply(genome_genes_sample$Sequence, AIC4)

AIC5<- function(seq) {
  sample <- unlist(strsplit(seq, ""))
  sample=gsub("A",1,sample)
  sample=gsub("C",2,sample)
  sample=gsub("G",3,sample)
  sample=gsub("T",4,sample)
  sample=as.numeric(sample)
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
  valid=nabcdef>0 & nabcde. >0
  aic=2*3072-2*sum(nabcdef[valid]*log(nabcdef[valid]/nabcde.[valid]))
  return(aic)
}
genome_cds_sample$AIC5 <- sapply(genome_cds_sample$Sequence, AIC5)
genome_genes_sample$AIC5 <- sapply(genome_genes_sample$Sequence, AIC5)

#Order with the minimum AIC for each sequence
genome_cds_sample$min_AIC_order <- apply(genome_cds_sample[, c("AIC0","AIC1", "AIC2", "AIC3","AIC4","AIC5")], 1, function(x) {
  order <- which.min(x)-1 
  return(order)  
})
genome_genes_sample$min_AIC_order <- apply(genome_genes_sample[, c("AIC0","AIC1", "AIC2", "AIC3","AIC4","AIC5")], 1, function(x) {
  order <- which.min(x)-1 
  return(order)  
})

#Contingency Table for BIC and AIC
table(genome_cds_sample$min_BIC_order,genome_cds_sample$min_AIC_order)
table(genome_genes_sample$min_BIC_order,genome_genes_sample$min_AIC_order)

#Save files
write.csv(genome_genes, file = "all human genome genes.csv", row.names = FALSE)
write.csv(genome_cds, file = "All Human Genome Coding DNA Sequences.csv", row.names = FALSE)
write.csv(genome_genes_sample, file = "Sample of Gene Sequences (Processed).csv", row.names = FALSE)
write.csv(genome_cds_sample, file = "Sample of Coding DNA Sequences (Processed).csv", row.names = FALSE)