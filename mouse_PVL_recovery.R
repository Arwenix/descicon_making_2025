
# PVL recovery
install.packages("pacman")
pacman::p_load(extraDistr, R2jags, parallel, ggpubr)

set.seed(1983)

### NB! Don't forget to set your working directory
setwd("/work/desccionmaking/mousedata/copy_module_3(1)/PVL")

# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}

print("done")
#------ create task environment -------------------
# NB! mod(ntrials, nstruct) (aka. ntrials %% nstruct) must be 0
ntrials <- 100 # total number of trials in our payoff structure
nstruct <- 10 # makes shure the probabilities of wins and losses are correct for ever 10 trials

#probabilities of loss
p1_prob_loss <- 0.1 
p2_prob_loss <- 0.2
p3_prob_loss <- 0.5 
p4_prob_loss <- 0.6 

#probabilities of win
p1_prob_win <- 0.9 
p2_prob_win <- 0.8 
p3_prob_win <- 0.5 
p4_prob_win <- 0.4 

good_freq_w <- 1 # 
good_infreq_w <- 2 
bad_freq_w <- 3 # p3
bad_infreq_w <- 4 # 

#good_freq_l <- -1 #
#good_infreq_l <-2  # 
#bad_freq_l <- -4 # 
#bad_infreq_l <- -3 # 

#testing out if the model is better with the actual punishments 
good_freq_l <- -0.363401247337694 # p1
good_infreq_l <-0.726802494675388  # p2
bad_freq_l <- 	-2.90720997870155 # p4
bad_infreq_l <- -2.18040748402616 # p3

# p3 - 50-50 loss win 
A_R <- c(rep(bad_freq_w, nstruct*p3_prob_win), rep(0,nstruct*(1-p3_prob_win))) # we win half of the time
A_L <- c(rep(bad_infreq_l, nstruct*p3_prob_loss),rep(0,nstruct*(1-p3_prob_loss))) # we have losses half of the time

# P4 - Frequent bad loss, high infrequent win 
B_R <- c(rep(bad_infreq_w, nstruct*p4_prob_win), rep(0,nstruct*(1-p4_prob_win))) #we win 40% of the time
B_L <- c(rep(bad_freq_l, nstruct*p4_prob_loss),rep(0,nstruct*(1-p4_prob_loss))) # loss 60% of the time

# p2 - frequent win, infrequent smol loss
C_R <- c(rep(good_infreq_w, nstruct*p2_prob_win), rep(0,nstruct*(1-p2_prob_win))) #wins 80% of the time 
C_L <- c(rep(good_infreq_l, nstruct*p2_prob_loss),rep(0,nstruct*(1-p2_prob_loss))) 

# Good frequent loss = P1
D_R <- c(rep(good_freq_w, nstruct*p1_prob_win), rep(0,nstruct*(1-p1_prob_win))) # wins 90% of the time
D_L <- c(rep(good_freq_l, nstruct*p1_prob_loss),rep(0,nstruct*(1-p1_prob_loss))) #losses

###

# create the pseudorandomized full payoff structure
A <- array(NA,ntrials) # setting up and empty array to be filled
B <- array(NA,ntrials)
C <- array(NA,ntrials)
D <- array(NA,ntrials)
for (i in 1:(ntrials/nstruct)) {
  A[(1+(i-1)*nstruct):(i*nstruct)] <- (A_R + sample(A_L)) # randomly shuffling the loss-array for every ten trials (and adding those losses to the winnings)
  B[(1+(i-1)*nstruct):(i*nstruct)] <- (B_R + sample(B_L))
  C[(1+(i-1)*nstruct):(i*nstruct)] <- (C_R + sample(C_L))
  D[(1+(i-1)*nstruct):(i*nstruct)] <- (D_R + sample(D_L))
}


payoff <- cbind(A,B,C,D)/100 # combining all four decks as columns with each 100 trials - dividing our payoffs by 100 to make the numbers a bit easier to work with

# let's look at the payoff
colSums(payoff) # the two bad decks should sum to -25 (i.e. -2500), and the two good ones to 25 (i.e. 2500)

#-------test PVL delta function and jags script ---------

#---set params

w <- 1.25 # weighting parameter (aka. loss aversion)
A <- .18 # shape parameter (aka. risk aversion)
theta <- 0.55 # inverse heat parameter (aka. choice consitency)
a <- .12 # learning rate parameter (aka. prediction error weighting)

#ntrials <- 100 # we have already specified this earlier (see line 15)

source("PVL.R")
PVL_sims <- PVL(payoff,ntrials,w,A,a,theta)

par(mfrow=c(2,2))
plot(PVL_sims$Ev[,1], ylim=c(-1,1))
plot(PVL_sims$Ev[,2], ylim=c(-1,1))
plot(PVL_sims$Ev[,3], ylim=c(-1,1))
plot(PVL_sims$Ev[,4], ylim=c(-1,1))
title(paste("Traces of expeceted value (Ev) for all four decks over", ntrials, "trials"), line = -1, outer = TRUE)

x <- PVL_sims$x
X <- PVL_sims$X

# set up jags and run jags model
data <- list("x","X","ntrials") 
params<-c("w","A","theta","a")
temp_samples <- jags.parallel(data, inits=NULL, params,
                              model.file ="PVL.txt",
                              n.chains=3, n.iter=5000, n.burnin=1000, n.thin=1,
                              n.cluster=3)

recov_w <- temp_samples$BUGSoutput$sims.list$w
recov_A <- temp_samples$BUGSoutput$sims.list$A
recov_a <- temp_samples$BUGSoutput$sims.list$a
recov_theta <- temp_samples$BUGSoutput$sims.list$theta

par(mfrow=c(4,1))

################ run this code instead since its complaining about the margin size
# Set up 2x2 grid and adjust margins
par(mfrow = c(2, 2), mar = c(2, 2, 2, 1))  

# Plot density distributions
plot(density(recov_w), main = "Density of w", xlab = "", ylab = "")
plot(density(recov_A), main = "Density of A", xlab = "", ylab = "")
plot(density(recov_a), main = "Density of a", xlab = "", ylab = "")
plot(density(recov_theta), main = "Density of theta", xlab = "", ylab = "")

# Add an overall title
mtext(paste("Density plots (for recovered w, A, a & theta) with ntrials =", ntrials), 
      outer = TRUE, cex = 1.2, line = -1)




###--------------Run full parameter recovery -------------
niterations <- 100 # fewer because it takes too long
true_w <- array(0,c(niterations))
true_A <- array(0,c(niterations))
true_a <- array(0,c(niterations))
true_theta <- array(0,c(niterations))

infer_w <- array(0,c(niterations))
infer_A <- array(0,c(niterations))
infer_a <- array(0,c(niterations))
infer_theta <- array(0,c(niterations))

# checking the runtime on our parameter recovery
start_time = Sys.time()

for (i in 1:niterations) {
  
  # let's see how robust the model is. Does it recover all sorts of values?
  w <- runif(1, 0.04372868, 18.52924) #inputting the 95% values we recovered from the data
  A <- runif(1, 0.00000001, 1.1969)
  a <- runif(1, 0.00240881, 0.8441233)
  theta <- runif(1, 0.05673509, 1.237673) 

  PVL_sims <- PVL(payoff,ntrials,w,A,a,theta)
  
  x <- PVL_sims$x
  X <- PVL_sims$X
  
  # set up jags and run jags model
  data <- list("x","X","ntrials") 
  params<-c("w","A","a","theta")
  samples <- jags.parallel(data, inits=NULL, params,
                           model.file ="PVL.txt",
                           n.chains=3, n.iter=5000, n.burnin=1000, n.thin=1,
                           n.cluster=3)
  
  
  true_w[i] <- w
  true_A[i] <- A
  true_theta[i] <- theta
  true_a[i] <- a
  
  # find maximum a posteriori
  Q <- samples$BUGSoutput$sims.list$w
  infer_w[i] <- MPD(Q)
  # infer_w[i] <-density(Q)$x[which(density(Q)$y==max(density(Q)$y))]
  
  Q <- samples$BUGSoutput$sims.list$A
  infer_A[i] <- MPD(Q)
  # infer_A_A[i] <-density(Q)$x[which(density(Q)$y==max(density(Q)$y))]
  
  Q <- samples$BUGSoutput$sims.list$a
  infer_a[i] <- MPD(Q)
  # infer_a[i] <-density(Q)$x[which(density(Q)$y==max(density(Q)$y))]
  
  Q <- samples$BUGSoutput$sims.list$theta
  infer_theta[i] <- MPD(Q)
  # infer_theta[i] <-density(Q)$x[which(density(Q)$y==max(density(Q)$y))]
  
  print(i)
  
}

end_time = Sys.time()
print("Runtime for RW model: ")
end_time - start_time

# let's look at some scatter plots
par(mfrow=c(2,2))
plot(true_w,infer_w)
plot(true_A,infer_A, ylim=c(0,1))
plot(true_a,infer_a)
plot(true_theta,infer_theta)


# plotting code courtesy of Lasse
source('/work/desccionmaking/mousedata/copy_module_3(1)/recov_plot.R')
pl1 <- recov_plot(true_w, infer_w, c("true w", "infer w"), 'smoothed linear fit')
pl2 <- recov_plot(true_A, infer_A, c("true A", "infer A"), 'smoothed linear fit')
pl3 <- recov_plot(true_a, infer_a, c("true a", "infer a"), 'smoothed linear fit')
pl4 <- recov_plot(true_theta, infer_theta, c("true theta", "infer theta"), 'smoothed linear fit')
ggarrange(pl1, pl2, pl3, pl4)