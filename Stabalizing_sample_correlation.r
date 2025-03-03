# ---------------------------------------------------------------------
# Simulation Study: Sampling Distribution of Sample Correlation (r_n)
# ---------------------------------------------------------------------
# Objective:
# This script performs a simulation to study the sampling distribution of 
# the sample correlation coefficient (r_n) and its variance-stabilized 
# transformation g(r_n) = arctanh(r_n).
#
# Steps:
# 1. Generate 1000 samples of size n ∈ {15, 25} from a bivariate normal 
#    distribution with mean (0,0), variances (1,1), and correlation ρ = 0.6.
# 2. Compute the sample correlation coefficient r_n for each sample.
# 3. Apply Fisher’s transformation g(r_n) = arctanh(r_n) to stabilize variance.
# 4. Store results and plot histograms of r_n and g(r_n) to compare their 
#    distributions.
# 5. Calculate the coverage probability of the confidence interval constructed 
# using the variance stabalizing function for ρ ∈ {0, 0.2, 0.4, 0.6, 0.8} 
#
# Expected Outcome:
# - The histogram of g(r_n) should be more symmetric and approximately normal, 
#   confirming the variance-stabilizing effect of Fisher’s transformation.
# ---------------------------------------------------------------------



# ________________________________________________set the coding environment
rm(list = objects()); gc()   # empty the global environment

options(
  java.parameters = c(
    "-XX:+UseG1GC",
    "-Xms3072m"
  )
) # Allocates more heap space for Java and uses a different garbage collector.


if(!require( "pacman")){
  install.packages("pacman")
}    # prepare package which helps package loading

pacman::p_load(
  MASS
  ) # load necessary packages

# __________________________________________________set the seed for reproducibility
set.seed(1133)


#____________________________________________________________sampling distributions
N <- 1000 # Number of random samples
n <-c(15,25) # Sample sizes

# Target parameters for univariate normal distributions
rho <- 0.6
mu1 <- 0; s1 <- 1
mu2 <- 0; s2 <- 1

# Parameters for bivariate normal distribution
mu <- c(mu1,mu2) # Mean
sigma <- matrix(c(s1^2, s1*s2*rho, s1*s2*rho, s2^2),
                2) # Covariance matrix

# Simulation of data
rn<-matrix(nrow=1000,ncol=2)
for (i in 1:N){
  for(k in 1:length(n)){
    bvn1 <- mvrnorm(n[k], mu = mu, Sigma = sigma ) # from MASS package
    rn[i,k] = cor(bvn1)[1,2]
  }
}

colnames(rn) <- c("rn15","rn25")

grn<-matrix(nrow=1000,ncol=2)
for (i in 1:N){
  for(k in 1:length(n)){
    bvn1 <- mvrnorm(n[k], mu = mu, Sigma = sigma ) # from MASS package
    grn[i,k] = atanh(cor(bvn1)[1,2])
  }
}

colnames(grn) <- c("grn15","grn25")

# plot histograms
par(mfrow=c(2, 2), family="Times", ps=9, oma = c(0,1,1,1) + 0.1,
    mar = c(4,3.8,2,2) + 0.1)
hist(rn[,1], main = "sample size = 15",xlab = "Rn")
hist(rn[,2], main = "sample size = 25",xlab = "Rn")
hist(grn[,1], main = "sample size = 15",xlab = "g(Rn)")
hist(grn[,2], main = "sample size = 25",xlab = "g(Rn)",xlim=c(0,1.5))

#_________________________________________________coverage probability
N <- 1000 # Number of random samples
n <-c(15,25) # Sample sizes
rho_v <- c(0,0.2,0.4,0.6,0.8)
mu1 <- 0; s1 <- 1
mu2 <- 0; s2 <- 1
mu <- c(mu1,mu2) # Mean

coverage<-matrix(nrow=N,ncol=length(n))
coverage.prob<-matrix(nrow=length(rho_v),ncol=length(n))

alpha.level = 0.05
z = qnorm(1-alpha.level/2)

for (j in 1:length(rho_v)){
  sigma <- matrix(c(s1^2, s1*s2*rho_v[j], s1*s2*rho_v[j], s2^2),2)
  rho <- rho_v[j]
  for (i in 1:N){
    for(k in 1:length(n)){
      bvn1 <- mvrnorm(n[k], mu = mu, Sigma = sigma ) # from MASS package
      grn = atanh(cor(bvn1)[1,2])
      u.ci = tanh(grn+z/sqrt(n[k]))
      l.ci = tanh(grn-z/sqrt(n[k]))
      coverage[i,k]= rho>=l.ci & rho<=u.ci

    }
  }
coverage.prob[j,1] = mean(coverage[,1])
coverage.prob[j,2]= mean(coverage[,2])
}

coverage.prob
