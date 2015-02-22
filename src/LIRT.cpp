#include <RcppArmadillo.h>
#include <eliot.h>

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List LIRTcpp(SEXP LambdaMatrix, SEXP VoteMatrix, SEXP IRTprec, SEXP StartIdealPoints, SEXP SaveChains, SEXP BurntIterations, SEXP SavedIterations){

Rcpp::RNGScope scope;
	
Rcpp::NumericMatrix Lambda(LambdaMatrix);  
Rcpp::IntegerMatrix Y(VoteMatrix);  
int ntopics = Lambda.ncol(); 
double prec = as<double>(IRTprec);
Rcpp::NumericVector initIdeal(StartIdealPoints);
Rcpp::LogicalVector store(SaveChains);
int burn = as<int>(BurntIterations);
int mcmc = as<int>(SavedIterations);
int nvotes = Y.nrow(); 
int nvoters = Y.ncol(); 

Rprintf("%d voters\n",nvoters);
Rprintf("%d votes\n",nvotes);
    R_FlushConsole();
    R_ProcessEvents(); 

// SET UP CHAIN STORAGE OBJECTS

int idealchainstorageloc = 0;	
int idealchainlength = nvoters*ntopics;
if (store[2]) idealchainlength = idealchainlength*mcmc;
Rcpp::NumericVector idealchain(idealchainlength);
		
int discrimchainstorageloc = 0;			
int discrimchainlength = nvotes*2;
if (store[3]) discrimchainlength = discrimchainlength*mcmc;
Rcpp::NumericVector discrimchain(discrimchainlength);

Rcpp::NumericVector rhochain(mcmc); // chain for rho
	
Rcpp::NumericVector llirtchain(mcmc); // log-likelihood chain for IRT model
	
Rcpp::NumericMatrix IRTystar(nvotes,nvoters);					// working latent utilities for each voter vote
Rcpp::NumericMatrix IRTtheta(nvoters,ntopics); 					// working voter ideal points for each topic
arma::mat IRTthetaarma(nvoters,ntopics);
Rcpp::NumericMatrix IRTbeta(nvotes,2);							// working discrimination parameters for each vote

// INITIALIZE WORKING VARIABLES	FOR IRT

double muim = 0;
Rcpp::NumericVector drawtemp(1);

arma::mat Xs(nvoters,2);
arma::mat Ys(nvoters,1);
arma::mat Tinvs(2,2);
	Tinvs(0,0) = prec;
	Tinvs(0,1) = 0;
	Tinvs(1,0) = 0;
	Tinvs(1,1) = prec;
arma::mat discrimvcovmat(2,2);
arma::mat discrimmumat(2,1); 
arma::mat discrimcholmat(2,2);
arma::mat abmat(2,1);
Rcpp::NumericVector IRTtempAB(2);

arma::mat Bs(nvotes,ntopics);
arma::mat Ws(nvotes,1);
arma::mat V(ntopics,ntopics);
double rho = 0.9;
V.fill(rho);
for (int k = 0; k < ntopics; k++) V(k,k) = 1;
arma::mat Vinvs = inv(V);
arma::mat idealmumat(ntopics,1);
arma::mat idealvcovmat(ntopics,ntopics);
arma::mat idealcholmat(ntopics,ntopics);
arma::mat thetamat(ntopics,1);
Rcpp::NumericVector IRTtempTheta(ntopics);

double proposalprec = 101;
double rhoproposal = 0.9;
double metropolisratio = 1;
double forwardjumplogdensity = 1;
double backwardjumplogdensity = 1;
double llcurrent = 1;
double llproposal = 1;
double burnacceptrate = 0;
double mcmcacceptrate = 0;
arma::mat llkerneltemp(1,1);
Rcpp::IntegerVector BurnAccepts(burn);
Rcpp::IntegerVector MCMCAccepts(mcmc);

double IRTll = 0;
Rcpp::NumericVector IRTllterm(1);
Rcpp::NumericVector zerovec(1);
zerovec(0) = 0;

Rprintf("Initializing Ideal Point Model...\n");
    R_FlushConsole();
    R_ProcessEvents();    

// INITIALIZE IDEAL POINT MODEL

for (int i = 0; i < nvoters; i++){ // for all voters...
    IRTtheta(i,0) = initIdeal(i);
	for (int k = 1; k < ntopics; k++){ // for all topics...
		IRTtheta(i,k) = IRTtheta(i,0);
	}
}

for (int m = 0; m < nvotes; m++){ // for all votes...
	IRTbeta(m,0) = Rf_rnorm(0,1);
	IRTbeta(m,1) = Rf_rnorm(0,1); 
}


Rprintf("Initiating Markov Chain Monte Carlo...\n");    
    R_FlushConsole();
    R_ProcessEvents();   

// BEGIN MARKOV CHAIN MONTE CARLO

for (int iter = 0; iter < burn + mcmc; iter ++){
 
	if (iter % 100 == 0) Rprintf("Beginning iteration %d of %d... \n",iter,burn+mcmc); // print simulation status
		
	// Draw truncated random normal for each vote for each voter

	for (int m = 0; m < nvotes; m++){ // for all votes...		
		for (int i = 0; i < nvoters; i++){ // for all voters...
			muim = -1.0*IRTbeta(m,0);
			for (int k=0; k < ntopics; k++) muim = muim + IRTbeta(m,1)* Lambda(m,k)*IRTtheta(i,k);
			if (Y(m,i) == 1){ // yes vote observed
				IRTystar(m,i) = eliot::rltnorm(muim,1.0,0.0);
			}
			if (Y(m,i) == 0){ // no vote observed
				IRTystar(m,i) = eliot::rutnorm(muim,1.0,0.0);
			}	
			if (Y(m,i) == -9){ // vote missing
				IRTystar(m,i) = Rf_rnorm(muim,1);
			}	
		}
	}
	
	
	// Draw normal for two discrimination parameters for each vote
	
	for (int m = 0; m < nvotes; m++){ // for all votes...	
		for (int i = 0; i < nvoters; i++){ // for all voters...			
			Xs(i,0) = -1;
			Xs(i,1) = 0;
			for (int k=0; k < ntopics; k++) Xs(i,1) = Xs(i,1) + Lambda(m,k)*IRTtheta(i,k);
			Ys(i,0) = IRTystar(m,i);
		}
		discrimvcovmat = inv_sympd(trans(Xs)*Xs + Tinvs);
		discrimmumat = discrimvcovmat*trans(Xs)*Ys;
		discrimcholmat = chol(discrimvcovmat);
		IRTtempAB = rnorm(2);
		for (int el = 0; el < 2; el++) abmat(el,0) = IRTtempAB(el);
		abmat = discrimmumat + discrimcholmat*abmat;
		for (int el = 0; el < 2; el++) IRTbeta(m,el) = abmat(el,0);
	}
	
	// Draw normal for each topic dimensions ideal points for each voter
	
	for (int m = 0; m < nvotes; m++){
		for (int k = 0; k < ntopics; k++){
			Bs(m,k) = Lambda(m,k)*IRTbeta(m,1);	
		}	
	}		
	for (int i = 0; i < nvoters; i++){ // for all voters...	
		for (int m = 0; m < nvotes; m++){
			Ws(m,0) = IRTystar(m,i) + IRTbeta(m,0);
		}
		idealvcovmat = inv_sympd(trans(Bs)*Bs + Vinvs);
		idealmumat = idealvcovmat*trans(Bs)*Ws;
		idealcholmat = chol(idealvcovmat);
		IRTtempTheta = rnorm(ntopics);
		for (int k = 0; k < ntopics; k++) thetamat(k,0) = IRTtempTheta(k);
		thetamat = idealmumat + idealcholmat*thetamat;
		for (int k = 0; k < ntopics; k++) {
			IRTtheta(i,k) = thetamat(k,0);
			IRTthetaarma(i,k) = IRTtheta(i,k);	
		}		
	}
	
	// Draw new value of correlation between ideal point dimensions via adaptive Metropolis step with beta proposal distribution
	// Note: adaptation of the proposal distribution begins at iteration 100 and continues until the end of burn in.
	
	rhoproposal = Rcpp::as<double>(rbeta(1,proposalprec*rho,proposalprec*(1-rho)));
	
	llcurrent = 0; 
	for (int i = 0; i < nvoters; i++){ // for all voters...	
			llkerneltemp = IRTthetaarma.row(i)*Vinvs*trans(IRTthetaarma.row(i));
			llcurrent = llcurrent - (1.0/2.0)*log(det(V)) - (1.0/2.0)*llkerneltemp(0,0);
	}
	
	V.fill(rhoproposal);
	for (int k = 0; k < ntopics; k++) V(k,k) = 1.0;
	Vinvs = inv(V);
	
	llproposal = 0;  
	for (int i = 0; i < nvoters; i++){ // for all voters...	
			llkerneltemp = IRTthetaarma.row(i)*Vinvs*trans(IRTthetaarma.row(i));
			llproposal = llproposal - (1.0/2.0)*log(det(V)) - (1.0/2.0)*llkerneltemp(0,0);
	}
	
	forwardjumplogdensity = Rf_lgammafn(proposalprec) - Rf_lgammafn(proposalprec*rho) - Rf_lgammafn(proposalprec*(1-rho)) + (proposalprec*rho - 1.0)*log(rhoproposal) + (proposalprec*(1-rho) - 1.0)*log(1-rhoproposal);
	backwardjumplogdensity = Rf_lgammafn(proposalprec) - Rf_lgammafn(proposalprec*rhoproposal) - Rf_lgammafn(proposalprec*(1-rhoproposal)) + (proposalprec*rhoproposal - 1.0)*log(rho) + (proposalprec*(1-rhoproposal) - 1.0)*log(1-rho);
	
	metropolisratio = exp(llproposal + backwardjumplogdensity - llcurrent - forwardjumplogdensity); 
	
	if (iter < burn){
		if (Rf_runif(0,1) < metropolisratio){
			BurnAccepts(iter) = 1;
			rho = rhoproposal;
		} else {
			BurnAccepts(iter) = 0;
		}
		
		if (iter > 100){
			burnacceptrate = 0;
			for (int iterate = iter-100; iterate < iter; iterate++) burnacceptrate = burnacceptrate + double(BurnAccepts(iterate))/100.0;
			if (burnacceptrate > 0.44){
				proposalprec = proposalprec*0.98;
			} else {
				proposalprec = proposalprec*1.02;
			}
		}	
		
	} else {
		if (Rf_runif(0,1) < metropolisratio){
			MCMCAccepts(iter - burn) = 1;
			rho = rhoproposal;
		} else {
			MCMCAccepts(iter - burn) = 0;
		}			
	}
			
	V.fill(rho);
	for (int k = 0; k < ntopics; k++) V(k,k) = 1;
	Vinvs = inv(V);
		
	// Compute log-likelihood for the voting data
	
	IRTll = 0;
	for (int m = 0; m < nvotes; m++){ // for all votes...		
		for (int i = 0; i < nvoters; i++){ // for all voters...
			muim = -1.0*IRTbeta(m,0);
			for (int k=0; k < ntopics; k++) muim = muim + IRTbeta(m,1)*Lambda(m,k)*IRTtheta(i,k);
			if (Y(m,i) == 1){ // yes vote observed
				IRTllterm = pnorm(zerovec,-muim,1.0);
				IRTll = IRTll + log(IRTllterm(0));
			}
			if (Y(m,i) == 0){ // no vote observed
				IRTllterm = pnorm(zerovec,muim,1.0);
				IRTll = IRTll + log(IRTllterm(0));
			}		
		}
	}
	
	
	
	if (iter >= burn){
				
		if (store[2]){
			for (int i = 0; i < nvoters; i++){
				for (int k = 0; k < ntopics; k++){
					idealchain(idealchainstorageloc) = IRTtheta(i,k);
					idealchainstorageloc++;
				}
			}	
		} else {
			for (int i = 0; i < nvoters; i++){
				for (int k = 0; k < ntopics; k++){
					idealchain(idealchainstorageloc) = idealchain(idealchainstorageloc) + IRTtheta(i,k)/mcmc;
					idealchainstorageloc++;
				}
			}
			idealchainstorageloc = 0;
		}
			
		if (store[3]){
			for (int j = 0; j < 2; j++){
				for (int m = 0; m < nvotes; m++){			
					discrimchain(discrimchainstorageloc) = IRTbeta(m,j);
					discrimchainstorageloc++;
				}
			}	
		} else {
				for (int j = 0; j < 2; j++){
					for (int m = 0; m < nvotes; m++){
					discrimchain(discrimchainstorageloc) = discrimchain(discrimchainstorageloc) + IRTbeta(m,j)/mcmc;
					discrimchainstorageloc++;
				}
			}
			discrimchainstorageloc = 0;
		}
			
	rhochain(iter - burn) = rho;		
			
	llirtchain(iter - burn) = IRTll;
	
	} 
    
    R_FlushConsole();
    R_ProcessEvents();
	
} // END MARKOV CHAIN MONTE CARLO

Rprintf("Markov Chain Monte Carlo Completed...\n");
	
// Load Mean Posterior Estimates into state variables for IRT model

idealchainstorageloc = 0;
if (store[2]){
	for (int i = 0; i < nvoters; i++){
		for (int k = 0; k < ntopics; k++){
			IRTtheta(i,k) = 0;
		}
	}	
	for (int iter = 0; iter < mcmc; iter++){	
		for (int i = 0; i < nvoters; i++){
			for (int k = 0; k < ntopics; k++){
				IRTtheta(i,k) = IRTtheta(i,k) + idealchain(idealchainstorageloc)/mcmc;
				idealchainstorageloc++;
			}
		}	
	}	
} else {
	for (int i = 0; i < nvoters; i++){
		for (int k = 0; k < ntopics; k++){
			IRTtheta(i,k) = idealchain(idealchainstorageloc);
			idealchainstorageloc++;
		}
	}	
}

discrimchainstorageloc = 0;
if (store[3]){
	for (int j = 0; j < 2; j++){
		for (int m = 0; m < nvotes; m++){
			IRTbeta(m,j) = 0;
		}
	}
	for (int iter = 0; iter < mcmc; iter++){		
		for (int j = 0; j < 2; j++){
			for (int m = 0; m < nvotes; m++){
				IRTbeta(m,j) = IRTbeta(m,j) + discrimchain(discrimchainstorageloc)/mcmc;
				discrimchainstorageloc++;
			}
		}
	}
} else {
	for (int j = 0; j < 2; j++){
		for (int m = 0; m < nvotes; m++){
			IRTbeta(m,j) = discrimchain(discrimchainstorageloc);
			discrimchainstorageloc++;
		}
	}	
}

// Calculate log-likelihood of IRT mean posterior in order to calculate the model DIC...
	
IRTll = 0;
for (int m = 0; m < nvotes; m++){ // for all votes...		
	for (int i = 0; i < nvoters; i++){ // for all voters...
		muim = -1.0*IRTbeta(m,0);
		for (int k=0; k < ntopics; k++) muim = muim + IRTbeta(m,1)*Lambda(m,k)*IRTtheta(i,k);
		if (Y(m,i) == 1){ // yes vote observed
			IRTllterm = pnorm(zerovec,-muim,1.0);
			IRTll = IRTll + log(IRTllterm(0));
		}
		if (Y(m,i) == 0){ // no vote observed
			IRTllterm = pnorm(zerovec,muim,1.0);
			IRTll = IRTll + log(IRTllterm(0));
		}		
	}
}

Rcpp::NumericVector dbar(1);	
Rcpp::NumericVector pd(1);	
Rcpp::NumericVector dic(1);	
dbar(0) = -2.0*mean(llirtchain);
pd(0) = dbar(0) + 2.0*IRTll;
dic(0) = dbar(0) + pd(0);

mcmcacceptrate = 0;
for (int iterate = 0; iterate < mcmc; iterate++) mcmcacceptrate = mcmcacceptrate + double(MCMCAccepts(iterate))/(1.0 + mcmc);	

return Rcpp::List::create( Rcpp::Named("ideal.chain") = idealchain, Rcpp::Named("discrim.chain") = discrimchain, Rcpp::Named("rho.chain") = rhochain, Rcpp::Named("ll.irt.chain") = llirtchain, Rcpp::Named("dbar") = dbar,Rcpp::Named("pd") = pd,Rcpp::Named("dic") = dic, Rcpp::Named("accept.rho") = mcmcacceptrate);

}



