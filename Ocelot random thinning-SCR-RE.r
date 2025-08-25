#==============================================================================#
#                                                                              #
#               Ocelot (Leopardus pardalis) density estimation                 #
#    El Cielo-Sierra de Tamalave Biological Corridor (Tamaulipas, México)      #
#                 RANDOM THINNING-SPATIAL CAPTURE-RECAPTURE                    #
#               Gabriela Mendoza-Gutiérrez, Leroy Soria-Díaz,                  #
#      Zavdiel A. Manuel-De la Rosa, José Jiménez, Nayeli Martínez-González,   #
#            Claudia C. Astudillo-Sánchez, Carlos Barriga-Vallejo              #
#                             19:28 24/08/2025                                 #
#                                                                              #
#==============================================================================#

setwd('C:/...')

# Load required libraries
library(scrbook)
library(secr)
library(lattice)
library(coda)
library(mcmcOutput)
library(nimble)
library(terra)

# Load custom SCR functions and samplers
source('Functions/SCR_functions.R')

# Load capture-recapture history
ocelote.ch <- read.capthist("BData/capt.txt", "BData/traps.txt", detector='count', cov="sex", noccasions=130)
summary(ocelote.ch)
secr::closure.test(ocelote.ch)

a<-attributes(ocelote.ch)
sex<-a$covariates[]
sex.indicator <- as.numeric(as.factor(sex$sex))
sex<- sex.indicator

# Rearrange capture history to 3D array: individuals x traps x occasions
y3d <- aperm(ocelote.ch, c(1, 3, 2))

# Extract trap locations and center coordinates
traplocs <- as.matrix(secr::traps(ocelote.ch))
X <- data.matrix(traplocs) / 1000
X[,1] <- X[,1] - mean(X[,1])
X[,2] <- X[,2] - mean(X[,2])
rownames(X) <- 1:52
colnames(X) <- c("X", "Y")

# Sampling parameters
(J <- nrow(X))             # Number of traps
(K <- dim(y3d)[3])         # Number of sampling occasions
(nind <- dim(y3d)[1])      # Number of detected individuals
n0 <- nind
(detections <- sum(y3d))   # Total number of detections

# Data augmentation
M <- 200                   # Augmented population size


################################################################################

# Calculate buffered state space size for irregular trap array
trapShape <- vect("C:/GIS/traps.shp")
buff_trap <- buffer(trapShape, width = 6267) # 3 * sigma buffer size
buffTrap <- aggregate(buff_trap)  # Combine geometries into single polygon
plot(buffTrap)                    # Plot buffered trap array
points(trapShape)                 # Add original trap locations to the plot
areaS <- expanse(buffTrap) / 1e6  # Calculate area in square kilometers


# Define scaled state space using buffer = 3*sigma
buff <- 6267/1000
xl <- min(X[,1]) - buff
xu <- max(X[,1]) + buff
yl <- min(X[,2]) - buff
yu <- max(X[,2]) + buff
(xlim = c(xl, xu))
(ylim = c(yl, yu))
area <- diff(xlim) * diff(ylim)

################################################################################

# Collapse 3D capture history to 2D: individuals x traps
y <- apply(y3d, c(1, 2), sum); sum(y)
yaug <- array(0, c(M, J))
yaug[1:nind, ] <- y

# Plot trap locations and capture intensity
plot(X, xlim=c(-11,11), ylim=c(-35,35), pch="+", asp=TRUE)
datn <- apply(y3d, c(2, 3), sum)
tot <- apply(datn, 1, sum)
symbols(X, circles=tot/1, inches=FALSE, bg="#00000022", fg=NULL, add=TRUE)
points(X, pch="+", cex=1)

# Load trap operation matrix
KT <- read.csv("BData/traps.csv", sep=",")
KT <- KT[,4:133] 
KT <- data.matrix(KT) # KT is a J x K matrix: 1 = active, 0 = inactive
colnames(KT) <- 1:130
KT                            # Trap activity across 130 days
(K <- ncol(KT))               # Number of occasions
(J <- nrow(KT))               # Number of traps

# Visualize trap activity: light blue = active, dark = inactive
image(1:K, 1:J, t(KT), yaxt = "n", xlab="Occasion", ylab="", cex.lab=1.25, col=topo.colors(2))
mtext(side = 2, "Camera trap", line = 2.5, cex=1.25)
axis(2, rev(seq(1, J, by=2)))
KT <- apply(KT, 1, sum)       # Total active days per trap

# Load non-identified (non-ID) capture frequencies
ocelote.un <- read.capthist("BData/NID.txt", "BData/traps.txt", detector='count', noccasions=130)
summary(ocelote.un)
y.un <- aperm(ocelote.un, c(1, 3, 2))

# Collapse non-ID data to trap x occasion
nnid <- apply(y.un, c(2, 3), sum)
sum(nnid)                     # Total non-ID events
nnidd <- apply(nnid, 1, sum)  # Non-ID events per trap


#### Model definition
NimModel <- nimbleCode({

  # Priors
  sig ~ dunif(0,10)       # Movement parameter
  psi ~ dbeta(1,1)        # Inclusion probability
  id.prob ~ dunif(0,1)    # Probability of individual identification
  # Random effects hyperparameter
  sigma.p ~ dunif(0, 10)
  mu0 ~ dnorm(0, 0.1)
  # Baseline detection rate
  for(j in 1:J){
    lp[j] ~ dnorm(mu0, sd=sigma.p)
    log(lam0[j]) <- lp[j]
  }
  
  for(i in 1:M) {
    z[i] ~ dbern(psi)       # Latent inclusion indicator
    s[i,1] ~ dunif(xlim[1], xlim[2])  # Activity center X
    s[i,2] ~ dunif(ylim[1], ylim[2])  # Activity center Y
	d2[i,1:J] <- (s[i,1] - X[1:J,1])^2 + (s[i,2] - X[1:J,2])^2
    # Detection rate per trap
    lam[i,1:J] <- GetDetectionRate(s = s[i,1:2], 
                                   X = X[1:J,1:2],
                                   J = J,
                                   sigma = sig,
                                   lam0 = lam0[1:J],
                                   z = z[i])

    # Full capture history (latent)
    y.full[i,1:J] ~ dPoissonVector(lambda = lam[i,1:J] * KT[1:J])
    
	for(j in 1:J) {
      # Observed capture history (ID only)
      y.obs[i,j] ~ dbin(id.prob, y.full[i,j])
      # Implement zero-trick for irregular trap array
      # https://groups.google.com/g/spatialcapturerecapture/c/NzqUovn8jF0/m/Plg2g6O6AgAJ   
      outj[i,j] <- sqrt(d2[i,j]) > buffer  # avoid MCMC sampling beyond buffer distance
    }
    out[i] <- equals(sum(outj[i,1:J]), J)  # Zero-trick condition
    zeros[i] ~ dbern(out[i])               # Zero-trick implementation
  }

  # Derived quantities
  N <- sum(z[1:M])          # Estimated population size
  D <- 100 * N / area       # Density per 100 km²
})


#### CONSTANTS
str(constants <- list(nnid = nnid,      # Non-ID events
                      J = J,            # Number of traps
                      M = M,            # Augmented individuals
                      KT = KT,          # Trap effort
                      buffer=buff,      # State space buffer
                      area = areaS      # State space area
))

#### DATA
SEX <- c(sex - 1, rep(NA, (M-nind)))
yred <- apply(yaug, c(1,2), sum)
str(data   <-    list(y.obs = yred,         # Observed ID histories
                      xlim = c(xl,xu),      # X limits of state space
                      ylim = c(yl,yu),      # Y limits of state space
                      zeros = rep(0, M),    # Zero-trick vector
                      X = X                 # Trap coordinates
))

#### INITIAL VALUES
s.start <- X[sample(1:J, size = M, replace = TRUE), ]
for(i in 1:nind){
  s.start[i,1] <- mean( X[yaug[i,]>0,1] )
  s.start[i,2] <- mean( X[yaug[i,]>0,2] )
}
d <- e2dist(s.start[1:M,], X)
lam0s <- runif(1, 0, 0.8)
sigs <- runif(1, 1, 3)
lam <- lam0s * exp(-(d^2) / (2 * sigs^2))

yi <- array(0, c(M, J, K)) # Latent resighting array
for (j in 1:J) {
  for (k in 1:K) {
    if (nnid[j, k] > 0) {
      probs <- lam[,j]
      probs <- probs / sum(probs)
      latent.id <- sample(1:M, nnid[j,k], prob = probs, replace = FALSE)
      yi[latent.id , j, k] <- 1
    }
  }
}

yis <- apply(yi, c(1,2), sum) + apply(yaug, c(1,2), sum)
zst <- apply(yis, 1, sum); zst[zst > 0] <- 1
id.prob.s <- sum(yaug) / (sum(yaug) + sum(nnid))

set.seed(1960)
str(inits   <-   list(mu0=runif(1,-2,0),
                      sigma.p=runif(1,0,1),
                      z = zst,
                      s = s.start,
                      sig = sigs,
                      id.prob = id.prob.s, 
                      psi = runif(1, 0.3, 0.9),
                      y.full = yis
))
#### Parameters to monitor
params <- c('psi', 'mu0','sigma.p', 'sig', 'N', 'D', 'id.prob')
params2 <- c('z','s')

#### Compile and run model
start.time <- Sys.time()
Rmodel <- nimbleModel(code = NimModel,
                      constants = constants,
                      data = data,
                      inits = inits,
                      check = FALSE)
# initial values for "complex" quantities (by O. Gimenez in
# https://gist.github.com/oliviergimenez/e41e9cb99174f2124f948308e19ca7ec)
simNodes <- 'lp'
simNodeScalar <- Rmodel$expandNodeNames(simNodes)
allNodes <- Rmodel$getNodeNames()
nodesSorted <- allNodes[allNodes %in% simNodeScalar]
set.seed(1) # This makes the simulations here reproducible
for(n in nodesSorted) {
  Rmodel$simulate(n)
  depNodes <- Rmodel$getDependencies(n)
  Rmodel$calculate(depNodes)
}

Rmodel$lp
Rmodel$calculate()
Rmodel$initializeInfo()
Cmodel <- compileNimble(Rmodel)

conf <- configureMCMC(Rmodel, monitors = params, monitors2=params2, 
                      thin = 5, thin2=10, enableWAIC=TRUE, useConjugacy = TRUE)

#### Add custom samplers
conf$removeSampler("y.full")
for(j in 1:J){
  conf$addSampler(target = paste("y.full[1:", M, ",", j, "]", sep=""),
                  type = 'IDSampler',
                  control = list(nnidd = nnidd[j], j = j, M = M),
                  silent = TRUE)
}

# Add block sampler for activity centers
conf$removeSampler("s")
ACnodes <- paste0("s[", 1:constants$M, ", 1:2]")
for(node in ACnodes) {
  conf$addSampler(target = node,
                  type = "RW_block",
                  control = list(adaptScaleOnly = TRUE),
                  silent = TRUE)
}

# Build and compile MCMC
Rmcmc <- buildMCMC(conf)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run MCMC
start.time2 <- Sys.time()
outNim <- runMCMC(Cmcmc, niter = 60000, nburnin = 10000, nchains = 3, 
                  inits = inits, setSeed = TRUE, progressBar = TRUE, 
                  samplesAsCodaMCMC = TRUE, WAIC=TRUE)
end.time <- Sys.time()
end.time - start.time2  # Runtime

summary(mcmcOutput(outNim$samples))



