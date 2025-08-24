

GetDetectionRate <- nimbleFunction(
  run = function(s = double(1), lam0=double(1), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
     d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
     ans <- lam0[1:J]*exp(-d2/(2*sigma^2))
     return(ans)
    }
  }
)

dPoissonVector <- nimbleFunction(
  run = function(x = double(1), lambda = double(1),
  log = integer(0, default = 0)) {
    J <- length(x)
    ans <- 0.0
    for(j in 1:J)
      ans <- ans + dpois(x[j], lambda[j], 1)
    returnType(double())
    if(log) return(ans)
    else return(exp(ans))
  })

rPoissonVector  <- nimbleFunction(
  run = function(n = integer(), lambda = double(1)) {
    J <- length(lambda)
    ans<- numeric(J)
    for(j in 1:J)
      ans[j] <- rpois(1, lambda[j])
    returnType(double(1))
    return(ans)
  })

registerDistributions(list(
  dPoissonVector = list(
    BUGSdist = "dPoissonVector(lambda)",
    Rdist = "dPoissonVector(lambda)",
    discrete = TRUE,
    range = c(0, Inf),
    types = c('value = double(1)', 'lambda = double(1)'))
))


#### Custom M-H sampler to jointly update y.un[1:M,j] so that they sum to n[j]
IDSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    nnidd <- control$nnidd
    j <- control$j
    M <- control$M
    calcNodes <- model$getDependencies(target)
  },

  run = function() {
    lam.curr <- model$lam[1:M,j]
    switch.probs <- lam.curr / sum(lam.curr)

    y.latent.curr <- model$y.full[1:M,j] - model$y.obs[1:M,j]
    y.latent.prop <- rmulti(1, nnidd, switch.probs)
    model$y.full[1:M,j] <<- model$y.obs[1:M,j] + y.latent.prop

    model_lp_initial <- model$getLogProb(calcNodes)
    model_lp_proposed <- model$calculate(calcNodes)

    log_MH_ratio <- (model_lp_proposed + dmulti(y.latent.curr, nnidd, switch.probs, log=TRUE)) -
                    (model_lp_initial  + dmulti(y.latent.prop, nnidd, switch.probs, log=TRUE))

    accept <- decide(log_MH_ratio)
    if(accept) {
      copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    } else {
      copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
    }
  },
  methods = list(reset = function () {})
)

