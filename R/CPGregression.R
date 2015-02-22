CPGregression <-
function (formula,data = environment(formula),twoEquations=FALSE,calculateErrors=TRUE,verbose=calculateErrors) {
	
library(tweedie)
library(statmod)
library(stats)
library(maxLik)
library(numDeriv)
	
llik.tweedie.terms <- function(par,whichterms){
	
	covariatecount <- dim(X)[2]
		
	zeta <- par[1]
	if (twoEquations == TRUE){
		gammavec <- par[2:(1+covariatecount)]
		betavec <- par[(2+covariatecount):(1+2*covariatecount)]
	} else {
		gammavec <- c(par[2],rep(0,covariatecount-1))
		betavec <- par[3:(2+covariatecount)]	
	}
	
		LL <- rep(NA,length(whichterms))
		for (i in 1:length(whichterms)){
			mui <- exp(sum(betavec * X[whichterms[i],]))
			phii <- (exp(((1-zeta)/(2))*sum(X[whichterms[i],]*(betavec - gammavec)))*exp(((2-zeta)/(2))*sum(X[whichterms[i],]*(betavec + gammavec))))/(2-zeta)
			LL[i] <- log(dtweedie(Y[whichterms[i]],zeta,mui,phii))
			}
		return(LL)
		
	}
	
llik.tweedie <- function(par){
	return(-sum(llik.tweedie.terms(par,1:N)))
	}
	
llik.tweedie.ci <- function(x,whichpar,par.est,robustratio=1){
	par <- replace(par.est,1:length(par.est) == whichpar,x)
	return(2*(llik.tweedie(par.est) - llik.tweedie(par))/robustratio + qchisq(0.95,1))
	}

	
	## Get starting values and model matrix using glm ##
	if (verbose) print("Finding MLE starting values using GLM...")
	power <-  suppressWarnings(tweedie.profile(formula,data=data,p.vec=seq(1.05,1.95,by=0.1),do.smooth=TRUE,method="interpolation"))
	zeta.start  <- power$p.max
	phi.start <- power$phi.max
	glmfit <- suppressWarnings(glm(formula,data=data,family=tweedie(var.power=zeta.start,link.power=0),x=TRUE))
	X <- glmfit$x
	Y <- glmfit$y	
	N <- length(Y)	
	K <- dim(X)[2]
	
	## need to calculate appropriate start value for gamma here
	gamma.start <- 2*log(phi.start)
	
	if (twoEquations){			
	start.values <- c(zeta.start,gamma.start,rep(0,K-1),glmfit$coef)
	lower.parameter.bounds <- c(1+1e-3,rep(-3e1,2*dim(X)[2]))
	upper.parameter.bounds <- c(2-1e-3,rep(3e1,2*dim(X)[2]))
	parameter.names <- c("Zeta",paste("Gamma",colnames(X)),paste("Beta",colnames(X)))
	} else {
	start.values <- c(zeta.start,gamma.start,glmfit$coef)
	lower.parameter.bounds <- c(1+1e-3,rep(-3e1,1+dim(X)[2]))
	upper.parameter.bounds <- c(2-1e-3,rep(3e1,1+dim(X)[2]))	
	parameter.names <- c("Zeta","Gamma (Intercept)",paste("Beta",colnames(X)))
	}
	
	if (verbose) print("Finding MLE...")	   
	nlm.out <- nlminb(start.values, llik.tweedie,lower=lower.parameter.bounds,upper=upper.parameter.bounds,control = list(trace=verbose))
	par.est <- nlm.out$par
	par.lci <- par.lci.rob <- rep(NA,length(par.est))
	par.hci <- par.hci.rob <- rep(NA,length(par.est))
	model.se <- robust.se <- rep(NA,length(par.est))

if (calculateErrors){
	
	if (verbose) print("Calculating Wald Standard Errors...")	   	
	jacobian.matrix <- jacobian(llik.tweedie.terms,par.est,whichterms=1:N)
	hessian.matrix <- hessian(llik.tweedie,par.est)
	D.matrix <- matrix(0,length(par.est),length(par.est))
	for (i in 1:N){
		D.matrix <- D.matrix + (1/N)*outer(jacobian.matrix[i,],jacobian.matrix[i,])
		}
	I.matrix <- (1/N)*hessian.matrix
	model.se <- sqrt(diag(solve(hessian.matrix)))
	robust.se <- sqrt(diag(solve(I.matrix) %*% D.matrix %*% solve(I.matrix) * (1/N)))
	C.matrix <- solve(I.matrix) %*% D.matrix %*% solve(I.matrix) / solve(I.matrix)
		
	if (verbose) print("Calculating Profile Confidence Intervals...")	
	
	lower.parameter.bounds <- c(1+1e-3,par.est[2:length(par.est)]-4*model.se[2:length(par.est)])
	upper.parameter.bounds <- c(2-1e-3,par.est[2:length(par.est)]+4*model.se[2:length(par.est)])
	
	for (k in 1:length(par.est)){
		ll.lower <- llik.tweedie.ci(lower.parameter.bounds[k],k,par.est)
		ll.est <- llik.tweedie.ci(par.est[k],k,par.est)
		ll.upper <- llik.tweedie.ci(upper.parameter.bounds[k],k,par.est)
		#print(c(ll.lower,ll.est,ll.upper))
		
		if (ll.lower < 0) {
		par.lci[k] <- suppressWarnings(uniroot(llik.tweedie.ci,interval=c(lower.parameter.bounds[k],par.est[k]),whichpar=k,par.est=par.est,robustratio=1))$root
		par.lci.rob[k] <- suppressWarnings(uniroot(llik.tweedie.ci,interval=c(lower.parameter.bounds[k],par.est[k]),whichpar=k,par.est=par.est,robustratio=C.matrix[k,k]))$root
		} else {
		par.lci[k] <- lower.parameter.bounds[k]
		}
		
		if (ll.upper < 0) {
		par.hci[k] <- suppressWarnings(uniroot(llik.tweedie.ci,interval=c(par.est[k],upper.parameter.bounds[k]),whichpar=k,par.est=par.est,robustratio=1))$root
		par.hci.rob[k] <- suppressWarnings(uniroot(llik.tweedie.ci,interval=c(par.est[k],upper.parameter.bounds[k]),whichpar=k,par.est=par.est,robustratio=C.matrix[k,k]))$root
		} else {
		par.hci[k] <- upper.parameter.bounds[k]	
		}
		
		}
	}
	
	Estimates <- cbind(par.est,model.se,robust.se,par.lci,par.hci,par.lci.rob,par.hci.rob)
	colnames(Estimates) <- c("Estimate","Model SE","Robust SE","95% CI Lower","95% CI Upper","Robust 95% CI Lower","Robust 95% CI Upper")
	rownames(Estimates) <- parameter.names
	out <- list(Estimates = Estimates,LogL=-llik.tweedie(par.est))
	if (verbose) print("Done.")  
    return(out)
}
