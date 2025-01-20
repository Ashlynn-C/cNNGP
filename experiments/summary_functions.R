library(LaplacesDemon)
library(coda)
library(spBayes)

crps <- function(y, y.hat, y.var){
    sd <- sqrt(y.var)
    y.std <- (y-y.hat)/sd
    mean(-sd*(1/sqrt(pi) - 2*dnorm(y.std) - y.std*(2*pnorm(y.std) - 1)))
}

rmspe <- function(y, y.hat){
    sqrt(mean((y-y.hat)^2))
}

results.df = function(true.y.pred, y.pred.samples, lp.pred.samples, corrected.intercept.samples, tau.sq.samples, sigma.sq.samples, phi.samples, beta.samples){
    
    res.names = c("crps", "waic", "rmspe", "mae", "cov","med.widths",
                  "tau.sq", "tau.sq.lb", "tau.sq.ub", "ESS.tau.sq",
                  "sigma.sq", "sigma.sq.lb", "sigma.sq.ub", "ESS.sigma.sq",
                  "phi", "phi.lb", "phi.ub", "ESS.phi",
                  "beta0", "beta0.lb", "beta0.ub", "ESS.beta0",
                  "beta0.corrected", "beta0.corrected.lb", "beta0.corrected.ub", "ESS.beta0.corrected",
                  "beta1", "beta1.lb", "beta1.ub", "ESS.beta1")
    
    res.df = data.frame(matrix(NA, nrow = 1, ncol = length(res.names),
                               dimnames = list(1, res.names)))
    
    model.y.pred = rowMeans(y.pred.samples)
    model.ints = apply(y.pred.samples,1,quantile,prob = c(0.025,0.975))
    
    res.df[,"crps"] = crps(true.y.pred, model.y.pred, apply(y.pred.samples,1,var))
    res.df[,"waic"] = WAIC(t(sapply(1:length(true.y.pred), function(i) 
        dnorm(true.y.pred[i], mean = lp.pred.samples[i,], sd=sqrt(tau.sq.samples),log = TRUE))))$WAIC
    res.df[,"rmspe"] = rmspe(true.y.pred, model.y.pred)
    res.df[,"mae"] = mean(abs(true.y.pred-model.y.pred))
    res.df[,"cov"] = mean(true.y.pred >= model.ints[1,] & true.y.pred <= model.ints[2,])
    res.df[,"med.widths"] = median(model.ints[2,] - model.ints[1,])
    res.df[,"tau.sq"] = mean(tau.sq.samples)
    res.df[,c("tau.sq.lb", "tau.sq.ub")] = quantile(tau.sq.samples, prob = c(0.025, 0.975))
    res.df[,"ESS.tau.sq"] = effectiveSize(tau.sq.samples)
    res.df[,"sigma.sq"] = mean(sigma.sq.samples)
    res.df[,c("sigma.sq.lb", "sigma.sq.ub")] = quantile(sigma.sq.samples, prob = c(0.025, 0.975))
    res.df[,"ESS.sigma.sq"] = effectiveSize(sigma.sq.samples)
    res.df[,"phi"] = mean(phi.samples)
    res.df[,c("phi.lb", "phi.ub")] = quantile(phi.samples, prob = c(0.025, 0.975))
    res.df[,"ESS.phi"] = effectiveSize(phi.samples)
    res.df[,"beta0"] = mean(beta.samples[1,])
    res.df[,c("beta0.lb", "beta0.ub")] = quantile(beta.samples[1,], prob = c(0.025, 0.975))
    res.df[,"ESS.beta0"] = effectiveSize(beta.samples[1,])
    res.df[,"beta0.corrected"] = mean(corrected.intercept.samples)
    res.df[,c("beta0.corrected.lb", "beta0.corrected.ub")] = quantile(corrected.intercept.samples, prob = c(0.025, 0.975))
    res.df[,"ESS.beta0.corrected"] = effectiveSize(corrected.intercept.samples)
    res.df[,"beta1"] = mean(beta.samples[2,])
    res.df[,c("beta1.lb", "beta1.ub")] = quantile(beta.samples[2,], prob = c(0.025, 0.975))
    res.df[,"ESS.beta1"] = effectiveSize(beta.samples[2,])
    
    return(res.df)
}

blockNNGP.summary = function(inla.model){
    
    blockNNGP.pred.time0 = proc.time()
    
    # find blocks new points belong to
    blocks.pred = matrix(NA, ncol=1, length(y.pred))
    k = 1
    for (j in 2:(n.partition+1)){
        for (i in 2:(n.partition+1)){
            liminfh <- hb[i-1]
            limsuph <- hb[i]
            liminfv <- vb[j-1]
            limsupv <- vb[j]
            indblock <- which(loc.pred[,1] >= liminfh & 
                                  loc.pred[,1] <= limsuph & 
                                  loc.pred[,2] >= liminfv & 
                                  loc.pred[,2] <= limsupv)
            blocks.pred[indblock] <- k
            k <- k + 1
        }
    }
    
    # get posterior samples
    blockNNGP.sample.time0 = proc.time()
    samples = inla.posterior.sample(length(keep.idx), inla.model, intern = F)
    blockNNGP.sample.time = proc.time() - blockNNGP.sample.time0
    hyperpar.samples = bind_rows(lapply(samples,'[[',"hyperpar"))
    colnames(hyperpar.samples) = c("tausq.inv", "sigma.sq", "phi")
    hyperpar.samples = hyperpar.samples %>% mutate(tau.sq = 1/tausq.inv, 
                                                   sigma.sq = exp(-sigma.sq),
                                                   phi = 30 - (29)/ (1 + exp(phi)))
    latent.samples = do.call(cbind,lapply(samples,'[[',"latent"))
    beta.samples = data.frame(latent.samples[rownames(latent.samples) %in% c("(Intercept):1","x:1"),])

    ### locations were ordering during model fitting, reordering here 
    ind1 = sort.int(blocks.fit, index.return=TRUE)
    w.samples = data.frame(t(latent.samples[rownames(latent.samples) %in% paste0("idx:",1:nrow(loc)),]))
    w.samples = w.samples[,order(ind1$ix)]
    
    blockNNGP.w.pred = matrix(NA, ncol = length(samples), nrow = length(y.pred))
    blockNNGP.y.pred = matrix(NA, ncol = length(samples), nrow = length(y.pred))
    blockNNGP.log.lik = matrix(NA, ncol = length(samples), nrow = length(y.pred))
    
    # get posterior samples for each point in holdout set
    for(block in 1:(n.partition^2)){
        
        ref.block.idx = which(blocks.fit == block)
        ref.block = loc[ref.block.idx,]
        block.D = iDist(ref.block)
        
        pred.block.idx = which(blocks.pred == block)
        
        if(length(pred.block.idx) == 0) {next}
        
        for(pred.point in 1:length(pred.block.idx)){
            
            pred.idx = pred.block.idx[pred.point]
            d.iN = iDist(loc.pred[pred.idx,,drop=F], loc[ref.block.idx,])
            
            w.pred.samples = sapply(1:ncol(blockNNGP.y.pred), function(i) {
                block.C = hyperpar.samples$sigma.sq[i]*exp(-hyperpar.samples$phi[i]*block.D)
                C.iN = hyperpar.samples$sigma.sq[i]*exp(-hyperpar.samples$phi[i]*t(d.iN))
                Bmat = t(C.iN)%*%solve(block.C)
                Fmat = hyperpar.samples$sigma.sq[i] - Bmat%*%C.iN
                rnorm(1, mean = Bmat%*%t(as.matrix(w.samples[i,ref.block.idx])), sd = sqrt(Fmat))
            })
            
            blockNNGP.w.pred[pred.idx,] = w.pred.samples
            blockNNGP.y.pred[pred.idx,] = rnorm(ncol(blockNNGP.y.pred), 
                                                mean = x.pred[pred.idx,]%*%as.matrix(beta.samples) + w.pred.samples, 
                                                sd = sqrt(hyperpar.samples$tau.sq))

        }
    }
    
    blockNNGP.pred.time = proc.time() - blockNNGP.pred.time0
    
    blockNNGP.res.df = results.df(y.pred,
                                  blockNNGP.y.pred,
                                  x.pred%*%as.matrix(beta.samples) + blockNNGP.w.pred,
                                  drop(as.matrix(beta.samples[1,] + colMeans(blockNNGP.w.pred))),
                                  hyperpar.samples$tau.sq,
                                  hyperpar.samples$sigma.sq,
                                  hyperpar.samples$phi, 
                                  as.matrix(beta.samples))
    
    blockNNGP.res.df["pred.time"] = blockNNGP.pred.time["elapsed"]
    blockNNGP.res.df["post.sample.time"] = blockNNGP.sample.time["elapsed"]
    
    return(blockNNGP.res.df)
}
