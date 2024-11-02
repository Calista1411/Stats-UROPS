#fonts
install.packages("extrafont")
library(extrafont)
font_import() 
loadfonts(device = "pdf")
par(family = "CM Roman") 

#for chr1
all_cds=read.csv("all_cds.csv",sep=",",header=TRUE)
all_genes<-read.csv("chr1_all_genes.csv",sep=",",header=TRUE)

#for whole genome
file_list <- list.files(path = "~/Downloads/urops/csv files of all genes", pattern = "*.csv", full.names = TRUE)
genome_genes <- do.call(rbind, lapply(file_list, read.csv))
head(genome_genes)
genome_genes <- genome_genes[order(as.numeric(genome_genes$Chr)), ]
names(genome_genes)
genome_genes <- genome_genes[genome_genes$Sequence != "" & !is.na(genome_genes$Sequence), ] #remove rows with no gene name

#sequence lengths 
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
max(genome_genes$length)
min(genome_genes$length)
median(genome_genes$length)
mean(genome_genes$length)
genome_genes <- genome_genes[order(as.numeric(genome_genes$length)), ]

#remove short segments
genome_genes_trimmed<-genome_genes[which(genome_genes$length>32758),]
set.seed(1)
genome_genes_trimmed2<-genome_genes_trimmed[sample(nrow(genome_genes_trimmed), 1000), ]#random sample of 1k rows
genome_genes_trimmed3<-genome_genes[sample(nrow(genome_genes_trimmed),1000),]
genome_cds_trimmed2 <- genome_cds[genome_cds$CDS.Name %in% genome_genes_trimmed2$Gene.Name, ]
genome_cds_trimmed2 <- merge(genome_cds, genome_genes_trimmed2, by = c("genome_genes_trimmed2$Gene.Name", "genome_genes_trimmed2$Chr"))
genome_cds <- genome_cds %>% rename(Gene.Name = "CDS.Name")
library(dplyr)
genome_genes_matched <- genome_genes %>%
  semi_join(genome_cds, by = c("Gene.Name", "Chr")) #match chr number and name
set.seed(1)
genome_cds_sample<-genome_cds[sample(nrow(genome_cds), 1000), ] #take random sample of 1k
genome_genes_sample<-genome_genes %>%
  semi_join(genome_cds_sample, by = c("Gene.Name", "Chr")) #take the same random sample of 1k

#function to generate combinations 
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

#summary of order using likelihood ratio test
sum(is.na(genome_genes_trimmed2$pvalue_1vs2)) #no of zero order=41 / 0
sum(is.na(genome_genes_trimmed2$pvalue_2vs3))-sum(is.na(genome_genes_trimmed2$pvalue_1vs2)) #no of first order=718 / whole genome: 704
sum(is.na(genome_genes_trimmed2$pvalue_3vs4))-sum(is.na(genome_genes_trimmed2$pvalue_2vs3)) #no of second order=466 / 2190
sum(is.na(genome_genes_trimmed2$pvalue_4vs5))-sum(is.na(genome_genes_trimmed2$pvalue_3vs4)) #no of third order=727 / 3017
#adds up to 1952. 
length(which(genome_genes_sample$pvalue_4vs5_0.00001<0.00001)) #16 remaining with p-value still <0.05
#no of fourth order=81 / 941
order <- c("Order 0", "Order 1", "Order 2", "Order 3","Order 4")
number_of_genes <- c(sum(is.na(genome_genes_sample$pvalue_1vs2_0.00001)),
                   sum(is.na(genome_genes_sample$pvalue_2vs3_0.00001)) - sum(is.na(genome_genes_sample$pvalue_1vs2_0.00001)),
                   sum(is.na(genome_genes_sample$pvalue_3vs4_0.00001)) - sum(is.na(genome_genes_sample$pvalue_2vs3_0.00001)),
                   sum(is.na(genome_genes_sample$pvalue_4vs5_0.00001)) - sum(is.na(genome_genes_sample$pvalue_3vs4_0.00001)), 46
)
# Bar plot
# Load fonts for use with pdf

# Set the font family to Computer Modern
library(showtext)
font_add("CM Roman", regular = "/devinacalista/Downloads/computer-modern/cmunrm.ttf")
showtext_auto()

# Generate the barplot
par(mar = c(5, 4, 4, 2) + 0.1)
barplot(number_of_CDS, 
        names.arg = order,       # Category labels
        col = "pink",         # Color of the bars
        xlab = "Order",          # X-axis label
        ylab = "Number of CDS",  # Y-axis label
        border = "black",
        family="CM Roman" 
)
mtext("Markov Chain Order of Coding DNA Sequences", side = 3, line = 2, cex = 1.2, family = "CM Roman", font = 2)

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
legend_x <- par("usr")[1] + 1  # Adjust x position slightly right from the left edge
legend_y <- par("usr")[4] - 0.5  # Keep y position near the top
legend(legend_x, legend_y, legend = c("0.1", "0.05", "0.01", "0.001", "0.00001"),
       fill = colors, cex = 0.85, title = "Significance Level")

#p values distribution
ks.test(genome_cds_sample$pvalue_indep_0.1,punif,0,1)  
ks.test(genome_cds_sample$cds_1vs2,punif,0,1,exact=TRUE) 
ks.test(genome_cds_sample$cds_2vs3,punif,0,1,exact=TRUE) 
ks.test(genome_cds_sample$cds_3vs4,punif,0,1,exact=TRUE) 
ks.test(genome_cds_sample$cds_4vs5,punif(0,1)) 
install.packages("nortest")
library(nortest)
ad.test(genome_cds_sample$pvalue_indep_0.1)

ks.test(genome_genes_sample$pvalue_indep_0.1,punif,0,1)
ks.test(genome_genes_sample$gene_1vs2,punif,0,1) #0.004
ks.test(genome_genes_sample$gene_2vs3,punif(0,1)) 
ks.test(genome_genes_sample$gene_3vs4,punif(0,1)) 
ks.test(genome_genes_sample$gene_4vs5,punif(0,1))

cds_1vs2<-genome_cds_sample$pvalue_1vs2_0.05
indices=is.na(cds_1vs2)
cds_1vs2[indices]<-0
ks.test(cds_1vs2,punif(0,1))

gene_1vs2<-genome_genes_sample$pvalue_1vs2_0.05
indices=is.na(gene_1vs2)
gene_1vs2[indices]<-punif(0,1)
ks.test(gene_1vs2,punif(0,1))

ks.test(genome_cds_sample$pvalue_indep_0.1,punif,0,1)  
ks.test(genome_cds_sample$pvalue_1vs2_0.05,punif,0,1,exact=TRUE) 
ks.test(genome_cds_sample$pvalue_2vs3_0.05,punif,0,1,exact=TRUE) 
ks.test(genome_cds_sample$pvalue_3vs4_0.05,punif,0,1,exact=TRUE) 
ks.test(genome_cds_sample$pvalue_4vs5_0.05,punif,0,1,exact=TRUE) 

#hypo tests for order of markov chain
#first order #all cds and genes asymptotically ok 
library(dplyr)
first_order <- function(cds) {
  sample <- unlist(strsplit(cds, ""))
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
  #pvalue <- mean(abs(results) >= abs(test))
  pvalue=pchisq(test,df=9,lower.tail=FALSE)
  return(pvalue)
}
genome_cds_sample$pvalue_indep_0.1 <- sapply(genome_cds_sample$Sequence, first_order)

#second order #those below 350, use empirical 
second_order <- function(cds,length) {
  sample <- unlist(strsplit(cds, ""))
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
    pvalue <- pchisq(test, df = 36, lower.tail = FALSE)
  } else {
    pvalue <- mean(abs(results1vs2) >= abs(test))
  }
  #pvalue <- pchisq(test, df = 36, lower.tail = FALSE)
  return(pvalue)
}
genome_genes_sample$gene_1vs2<-mapply(second_order, 
                                   genome_genes_sample$Sequence,
                                   genome_genes_sample$length)
genome_cds_sample$pvalue_1vs2_0.1 <- NA
indices <- which(genome_cds_sample$pvalue_indep_0.1 < 0.1)
#genome_cds_sample$pvalue_1vs2_0.1 <- sapply(genome_cds_sample$Sequence, second_order)
genome_cds_sample$pvalue_1vs2_0.1[indices] <- mapply(second_order, 
                                                             genome_cds_sample$Sequence[indices],
                                                             genome_cds_sample$length[indices])

genome_genes_sample$pvalue_1vs2_0.05<-genome_genes_sample$pvalue_1vs2_0.1 #for alpha=0.05
indices <- which(genome_genes_sample$pvalue_indep_0.1 > 0.05)
genome_genes_sample$pvalue_1vs2_0.05[indices] <- NA

genome_cds_sample$pvalue_1vs2_0.00001<-genome_cds_sample$pvalue_1vs2_0.1 
indices <- which(genome_cds_sample$pvalue_indep_0.1 > 0.00001)
genome_cds_sample$pvalue_1vs2_0.00001[indices] <- NA

#third order #below 5000 use empirical
third_order <- function(cds,length) {
  sample <- unlist(strsplit(cds, ""))
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
genome_cds_sample$pvalue_2vs3_0.1 <- NA
indices <- which(genome_cds_sample$pvalue_1vs2_0.1 < 0.1)
#genome_cds_sample$pvalue_2vs3_0.1 <- sapply(genome_cds_sample$Sequence, third_order)
genome_cds_sample$pvalue_2vs3_0.1[indices] <- mapply(third_order, 
                                                      genome_cds_sample$Sequence[indices],
                                                      genome_cds_sample$length[indices]) #for simulation
genome_genes_sample$cds_2vs3<-mapply(third_order, 
                                   genome_cds_sample$Sequence,
                                   genome_cds_sample$length)
genome_genes_sample$pvalue_2vs3_0.00001<-genome_genes_sample$pvalue_2vs3_0.1 #changing alpha
indices <- which(genome_genes_sample$pvalue_1vs2_0.00001 > 0.00001 | is.na(genome_genes_sample$pvalue_1vs2_0.00001))
genome_genes_sample$pvalue_2vs3_0.00001[indices] <- NA

#4th order #below 6000 use empirical
fourth_order <- function(cds,length) {
  sample <- unlist(strsplit(cds, ""))
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
  if (length(nabcde) < 1024) { # Ensure it has 1024 counts
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
  if (length(nbcde) < 256) { # Ensure it has 256 counts
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
  #pvalue <- pchisq(test, df = 576, lower.tail = FALSE)
  return(pvalue)
}
genome_cds_sample$pvalue_3vs4_0.1 <- NA
indices <- which(genome_cds_sample$pvalue_2vs3_0.1 < 0.1)
#genome_genes_sample$pvalue_3vs4_0.1 <- sapply(genome_genes_sample$Sequence, fourth_order)
genome_cds_sample$pvalue_3vs4_0.1[indices] <- mapply(fourth_order, 
                                                      genome_cds_sample$Sequence[indices],
                                                      genome_cds_sample$length[indices])
genome_cds_sample$cds_3vs4<-mapply(fourth_order, 
                                   genome_cds_sample$Sequence,
                                   genome_cds_sample$length)
genome_genes_sample$pvalue_3vs4_0.00001<-genome_genes_sample$pvalue_3vs4_0.1 #changing alpha
indices <- which(genome_genes_sample$pvalue_2vs3_0.00001 > 0.00001 | is.na(genome_genes_sample$pvalue_2vs3_0.00001))
genome_genes_sample$pvalue_3vs4_0.00001[indices] <- NA

#5th order
fifth_order <- function(cds) {
  sample <- unlist(strsplit(cds, ""))
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
  if (length(nabcdef) < 4096) { # Ensure it has 4096 counts
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
  if (length(nbcdef) < 1024) { # Ensure it has 1024 counts
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
  pvalue <- mean(abs(results24562) >= abs(test))
  #pvalue=pchisq(test,df=2304,lower.tail=FALSE)
  return(pvalue)
}
genome_cds_sample$pvalue_4vs5_0.1 <- NA
indices <- which(genome_cds_sample$pvalue_3vs4_0.1 < 0.1)
genome_cds_sample$pvalue_4vs5_0.1[indices] <- sapply(genome_cds_sample$Sequence[indices], fifth_order)

genome_genes_sample$pvalue_4vs5_0.00001<-genome_genes_sample$pvalue_4vs5_0.1 #for alpha=0.01
indices <- which(genome_genes_sample$pvalue_3vs4_0.00001 > 0.00001 | is.na(genome_genes_sample$pvalue_3vs4_0.00001))
genome_genes_sample$pvalue_4vs5_0.00001[indices] <- NA

genome_cds_sample$cds_4vs5<-sapply(genome_cds_sample$Sequence, fifth_order)

#BIC 
BIC0<-function(cds){
  sample <- unlist(strsplit(cds, ""))
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

BIC1 <- function(cds) {
  sample <- unlist(strsplit(cds, ""))
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

BIC2 <- function(cds) {
  sample <- unlist(strsplit(cds, ""))
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

BIC3 <- function(cds) {
  sample <- unlist(strsplit(cds, ""))
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

BIC4 <-function(cds) {
  sample <- unlist(strsplit(cds, ""))
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
  if (length(nabcde) < 1024) { # Ensure it has 1024 counts
    nabcde <- c(nabcde, rep(0, 1024 - length(nabcde)))
  }
  five_words_matrix <- matrix(nabcde, nrow = 4^4, ncol = 4,byrow=TRUE)
  nabcd.=rep(rowSums(five_words_matrix),rep(4,4^4))
  valid=nabcde>0 & nabcd.>0
  bic=768*log(length(sample))-2*sum(nabcde[valid]*log(nabcde[valid]/nabcd.[valid]))
  return(bic)
}
genome_cds_sample$BIC4 <- sapply(genome_cds_sample$Sequence, BIC4)

BIC5<- function(cds) {
  sample <- unlist(strsplit(cds, ""))
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
  if (length(nabcdef) < 4096) { # Ensure it has 4096 counts
    nabcdef <- c(nabcdef, rep(0, 4096 - length(nabcdef)))
  }
  six_words_matrix <- matrix(nabcdef, nrow = 4^5, ncol = 4,byrow=TRUE)
  nabcde.=rep(rowSums(six_words_matrix),rep(4,4^5))
  valid=nabcdef>0 & nabcde. >0
  bic=3072*log(length(sample))-2*sum(nabcdef[valid]*log(nabcdef[valid]/nabcde.[valid]))
  return(bic)
}
genome_cds_sample$BIC5 <- sapply(genome_cds_sample$Sequence, BIC5)

genome_cds_sample$min_BIC_order <- apply(genome_cds_sample[, c("BIC0","BIC1","BIC2","BIC3","BIC4","BIC5")], 1, function(x) {
  order <- which.min(x) - 1  # Subtract 1 to match the order (0 for BIC0, 1 for BIC1, etc.)
  return(order)  # Return 0, 1, 2, 3, 4, or 5 as the minimum BIC order
})
count_first_order <- sum(genome_genes_sample$min_BIC_order == 1) #2032 have first-order as min BIC // whole genome: 6760 // first 1000 sequences of genes length 32758+: 137
count_second_order <- sum(genome_genes_sample$min_BIC_order == 2) #10 // 129 // 849
count_third_order <- sum(genome_genes_sample$min_BIC_order == 3) #7 // 30 // 14
count_fourth_order <- sum(genome_genes_sample$min_BIC_order == 4) # // // 0

#valid=nabcde > 0 & nabcd. > 0
#fourth=768*log(length(sample))-2*sum(nabcde[valid]*log(nabcde[valid]/nabcd.[valid]))
#fifth=3072*log(length(sample))-2*sum(nabcdef[valid]*log(nabcdef[valid]/nabcde.[valid]))

#AIC
AIC0<-function(cds){
  sample <- unlist(strsplit(cds, ""))
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

AIC1 <- function(cds) {
  sample <- unlist(strsplit(cds, ""))
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

AIC2 <- function(cds) {
  sample <- unlist(strsplit(cds, ""))
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

AIC3 <- function(cds) {
  sample <- unlist(strsplit(cds, ""))
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

AIC4 <-function(cds) {
  sample <- unlist(strsplit(cds, ""))
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
  if (length(nabcde) < 1024) { # Ensure it has 1024 counts
    nabcde <- c(nabcde, rep(0, 1024 - length(nabcde)))
  }
  five_words_matrix <- matrix(nabcde, nrow = 4^4, ncol = 4,byrow=TRUE)
  nabcd.=rep(rowSums(five_words_matrix),rep(4,4^4))
  valid=nabcde>0 & nabcd.>0
  aic=2*768-2*sum(nabcde[valid]*log(nabcde[valid]/nabcd.[valid]))
  return(aic)
}
genome_cds_sample$AIC4 <- sapply(genome_cds_sample$Sequence, AIC4)

AIC5<- function(cds) {
  sample <- unlist(strsplit(cds, ""))
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
  if (length(nabcdef) < 4096) { # Ensure it has 4096 counts
    nabcdef <- c(nabcdef, rep(0, 4096 - length(nabcdef)))
  }
  six_words_matrix <- matrix(nabcdef, nrow = 4^5, ncol = 4,byrow=TRUE)
  nabcde.=rep(rowSums(six_words_matrix),rep(4,4^5))
  valid=nabcdef>0 & nabcde. >0
  aic=2*3072-2*sum(nabcdef[valid]*log(nabcdef[valid]/nabcde.[valid]))
  return(aic)
}
genome_cds_sample$AIC5 <- sapply(genome_cds_sample$Sequence, AIC5)

genome_cds_sample$min_AIC_order <- apply(genome_cds_sample[, c("AIC0","AIC1", "AIC2", "AIC3","AIC4","AIC5")], 1, function(x) {
  order <- which.min(x)-1 # Get the index of the minimum value (1 for BIC1, 2 for BIC2, 3 for BIC3)
  return(order)  # Return 1, 2, or 3 corresponding to the order
})
count_first_order <- sum(genome_cds_sample$min_AIC_order == 1) #1395 have first-order as min AIC // 2445 // 0
count_second_order <- sum(genome_cds_sample$min_AIC_order == 2) #594 // 3855 // 142
count_third_order <- sum(genome_genes_trimmed1$min_AIC_order == 3) #60 // 619 // 488
count_fourth_order <- sum(genome_genes_trimmed1$min_AIC_order == 4)# // 255
count_fifth_order <- sum(genome_genes_trimmed1$min_AIC_order == 5)#115

library(openintro)
table(genome_genes_trimmed2$min_BIC_order,genome_genes_trimmed2$min_AIC_order)
table(genome_cds_sample$min_BIC_order,genome_cds_sample$min_AIC_order)

all_cds=all_cds[,-c(17,18)]

#save file

all_cds <- all_cds[, c("CDS.Name", "Sequence", "GC_Content", "pvalue_indep", "pvalue_1vs2",
                       "pvalue_2vs3","pvalue_3vs4","pvalue_4vs5","BIC1","BIC2","BIC3",
                       "min_BIC_order","AIC1","AIC2","AIC3","min_AIC_order","var_dist_indepvs1",
                       "var_dist_1vs2","var_dist_2vs3","var_dist_3vs4","max_var_dist")]
write.csv(genome_cds_trimmed, file = "genome_cds_alldata.csv", row.names = FALSE)
write.csv(genome_genes, file = "all human genome genes.csv", row.names = FALSE)

#############################################################################

AIC4 <- function(cds) {
  sample <- unlist(strsplit(cds, ""))
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
  if (length(nabcde) < 1024) { # Ensure it has 1024 counts
    nabcde <- c(nabcde, rep(0, 1024 - length(nabcde)))
  }
  five_words_matrix <- matrix(nabcde, nrow = 4^4, ncol = 4,byrow=TRUE)
  nabcd.=rep(rowSums(five_words_matrix),rep(4,4^4))
  valid=nabcde>0 & nabcd. >0
  aic=-2*sum(nabcde[valid]*log(nabcde[valid]/nabcd.[valid]))
  return(aic)
}
all_cds$AIC4 <- sapply(all_cds$Sequence, AIC4)

AIC5<- function(cds) {
  sample <- unlist(strsplit(cds, ""))
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
  if (length(nabcdef) < 4096) { # Ensure it has 4096 counts
    nabcdef <- c(nabcdef, rep(0, 4096 - length(nabcdef)))
  }
  six_words_matrix <- matrix(nabcdef, nrow = 4^5, ncol = 4,byrow=TRUE)
  nabcde.=rep(rowSums(six_words_matrix),rep(4,4^5))
  valid=nabcdef>0 & nabcde. >0
  aic=-2*sum(nabcdef[valid]*log(nabcdef[valid]/nabcde.[valid]))
  return(aic)
}
all_cds$AIC5 <- sapply(all_cds$Sequence, AIC5)

#6th order
all_combinations <- generate_combinations(4, 7)
window_size <- 7
combinations <- embed(sample, window_size)[, window_size:1]
words <- apply(combinations, 1, paste, collapse = "")
seven_letter_count <- table(factor(words, levels = all_combinations))
nabcdefg <- as.vector(seven_letter_count)
seven_words_matrix <- matrix(nabcdefg, nrow = 4^6, ncol = 4,byrow=TRUE)
nabcdef.=rep(rowSums(seven_words_matrix),rep(4,4^6))

all_combinations <- generate_combinations(4, 6)
window_size2 <- 6
combinations2 <- embed(sample[-1], window_size2)[, window_size2:1]
words2 <- apply(combinations2, 1, paste, collapse = "")
six_letter_count <- table(factor(words2, levels = all_combinations))
nbcdefg <- as.vector(six_letter_count)
six_words_matrix=matrix(six_letter_count,nrow=4^5,ncol=4,byrow=TRUE)
nbcdef.=rep(rep(rowSums(six_words_matrix),rep(4,4^5)),4)
nbcdefg=rep(as.vector(six_letter_count),4)
expected=(nabcdef.*nbcdefg)/(nbcdef.)
valid_indices <- expected > 0
nabcdefg_valid <- nabcdefg[valid_indices]
expected_valid <- expected[valid_indices]
test=((nabcdefg_valid-expected_valid)^2)/expected_valid
valid_indices=!is.na(test)
test_valid=test[valid_indices]
sum(test_valid)
pchisq(sum(test_valid),df=9216,lower.tail=FALSE)
4^6*3-4^5*3

#7th order 
all_combinations <- generate_combinations(4, 8)
window_size <- 8
combinations <- embed(sample, window_size)[, window_size:1]
words <- apply(combinations, 1, paste, collapse = "")
eight_letter_count <- table(factor(words, levels = all_combinations))
nabcdefgh <- as.vector(eight_letter_count)
eight_words_matrix <- matrix(nabcdefgh, nrow = 4^7, ncol = 4,byrow=TRUE)
nabcdefg.=rep(rowSums(eight_words_matrix),rep(4,4^7))

all_combinations <- generate_combinations(4, 7)
window_size2 <- 7
combinations2 <- embed(sample[-1], window_size2)[, window_size2:1]
words2 <- apply(combinations2, 1, paste, collapse = "")
seven_letter_count <- table(factor(words2, levels = all_combinations))
nbcdefgh <- as.vector(seven_letter_count)
seven_words_matrix=matrix(seven_letter_count,nrow=4^6,ncol=4,byrow=TRUE)
nbcdefg.=rep(rep(rowSums(seven_words_matrix),rep(4,4^6)),4)
nbcdefgh=rep(as.vector(seven_letter_count),4)
expected=(nabcdefg.*nbcdefgh)/(nbcdefg.)
valid_indices <- expected > 0
nabcdefgh_valid <- nabcdefgh[valid_indices]
expected_valid <- expected[valid_indices]
test=((nabcdefgh_valid-expected_valid)^2)/expected_valid
valid_indices=!is.na(test)
test_valid=test[valid_indices]
sum(test_valid)
pchisq(sum(test_valid),df=36864,lower.tail=FALSE)
4^7*3-4^6*3

#8th order
all_combinations <- generate_combinations(4, 9)
window_size <- 9
combinations <- embed(sample, window_size)[, window_size:1]
words <- apply(combinations, 1, paste, collapse = "")
nine_letter_count <- table(factor(words, levels = all_combinations))
nabcdefghi <- as.vector(nine_letter_count)
nine_words_matrix <- matrix(nabcdefghi, nrow = 4^8, ncol = 4,byrow=TRUE)
nabcdefgh.=rep(rowSums(nine_words_matrix),rep(4,4^8))

all_combinations <- generate_combinations(4, 8)
window_size2 <- 8
combinations2 <- embed(sample[-1], window_size2)[, window_size2:1]
words2 <- apply(combinations2, 1, paste, collapse = "")
eight_letter_count <- table(factor(words2, levels = all_combinations))
nbcdefghi <- as.vector(eight_letter_count)
eight_words_matrix=matrix(eight_letter_count,nrow=4^7,ncol=4,byrow=TRUE)
nbcdefgh.=rep(rep(rowSums(eight_words_matrix),rep(4,4^7)),4)
nbcdefghi=rep(as.vector(eight_letter_count),4)
expected=(nabcdefgh.*nbcdefghi)/(nbcdefgh.)
valid_indices <- expected > 0
nabcdefghi_valid <- nabcdefghi[valid_indices]
expected_valid <- expected[valid_indices]
test=((nabcdefghi_valid-expected_valid)^2)/expected_valid
valid_indices=!is.na(test)
test_valid=test[valid_indices]
sum(test_valid)
pchisq(sum(test_valid),df=147456,lower.tail=FALSE)
4^8*3-4^7*3

#9th order
all_combinations <- generate_combinations(4, 10)
window_size <- 10
combinations <- embed(sample, window_size)[, window_size:1]
words <- apply(combinations, 1, paste, collapse = "")
ten_letter_count <- table(factor(words, levels = all_combinations))
nabcdefghij <- as.vector(ten_letter_count)
ten_words_matrix <- matrix(nabcdefghij, nrow = 4^9, ncol = 4,byrow=TRUE)
nabcdefghi.=rep(rowSums(ten_words_matrix),rep(4,4^9))

all_combinations <- generate_combinations(4, 9)
window_size2 <- 8
combinations2 <- embed(sample[-1], window_size2)[, window_size2:1]
words2 <- apply(combinations2, 1, paste, collapse = "")
nine_letter_count <- table(factor(words2, levels = all_combinations))
nbcdefghij <- as.vector(nine_letter_count)
nine_words_matrix=matrix(nine_letter_count,nrow=4^8,ncol=4,byrow=TRUE)
nbcdefghi.=rep(rep(rowSums(nine_words_matrix),rep(4,4^8)),4)
nbcdefghij=rep(as.vector(nine_letter_count),4)
expected=(nabcdefghi.*nbcdefghij)/(nbcdefghi.)
valid_indices <- expected > 0
nabcdefghij_valid <- nabcdefghij[valid_indices]
expected_valid <- expected[valid_indices]
test=((nabcdefghij_valid-expected_valid)^2)/expected_valid
valid_indices=!is.na(test)
test_valid=test[valid_indices]
sum(test_valid)
pchisq(sum(test_valid),df=147456,lower.tail=FALSE)


#AIC
first_order=2*12-2*sum(nij*log(nij/ni.))
second_order=2*48-2*sum(nijk*log(nijk/nij.))
third=2*192-2*sum(nabcd*log(nabcd/nabc.))
valid=nabcde > 0 & nabcd. > 0
fourth=2*768-2*sum(nabcde[valid]*log(nabcde[valid]/nabcd.[valid]))
valid=nabcdef > 0 & nabcde. >0
fifth=2*3072-2*sum(nabcdef[valid]*log(nabcdef[valid]/nabcde.[valid]))
valid=nabcdefg > 0 & nabcdef. >0
sixth=2*12288-2*sum(nabcdefg[valid]*log(nabcdefg[valid]/nabcdef.[valid]))
valid=nabcdefgh > 0 & nabcdefg. >0
seven=2*49152-2*sum(nabcdefgh[valid]*log(nabcdefgh[valid]/nabcdefg.[valid]))
valid=nabcdefghi > 0 & nabcdefgh. >0
eight=2*4^8*3-2*sum(nabcdefghi[valid]*log(nabcdefghi[valid]/nabcdefgh.[valid]))

#results
#ORF of a gene (around 92k bases)
  #7th order for chi-squared test, 3rd order for BIC, 5th order for AIC
  #BIC:1025951,256089,257088,261994,...
  #AIC:1025838,255635,255269,254717,253964,272819,...

#first 100k bases
  #6th order for chi-squared test, 2nd order for BIC, 4th order for AIC
  #BIC: 1074885,267376,268056,272909,...
  #AIC: 1074771,266920,266229,265604,266397,...

#first 1 million bases
  #4th order for BIC
  #BIC: 10837920,2685718,2674992,2666973,2674673,...
  #AIC: 10837778,2685150,2672723,2657898,2638376,2603777,2562099,...

#ORF of a gene (around 92k bases)
  #7th order for chi-squared test, 3rd order for BIC, 5th order for AIC
#first 100k bases
  #6th order for chi-squared test, 2nd order for BIC, 4th order for AIC
#first 1 million bases
  #4th order for BIC
#non coding region from position 5135311-5422642
  #2nd order for BIC, 5th order for AIC
  

#sample=ori_sample[123000000:123500000] #stable cg content part


