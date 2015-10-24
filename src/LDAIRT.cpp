#include <RcppArmadillo.h>
#include <RcppTN.h>
#include <eliot.h>

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List LDAIRTcpp( SEXP WordMatrixDocs, SEXP WordMatrixTerms, SEXP WordMatrixCounts, SEXP VoteMatrix, SEXP LDAtopics, SEXP LDAterms, SEXP LDAwords, SEXP LDAdocterms, SEXP LDAalpha, SEXP LDAbeta, SEXP IRTprec, SEXP StartIdealPoints, SEXP SaveChains, SEXP BurntIterations, SEXP SavedIterations ){

Rcpp::RNGScope scope;
	
Rcpp::IntegerVector WMD(WordMatrixDocs); 
Rcpp::IntegerVector WMT(WordMatrixTerms); 
Rcpp::IntegerVector WMC(WordMatrixCounts); 
Rcpp::IntegerMatrix Y(VoteMatrix);  
int ntopics = as<int>(LDAtopics); 
int nterms = as<int>(LDAterms);
int nwords = as<int>(LDAwords);
int ndocterms = as<int>(LDAdocterms); 
double alpha = as<double>(LDAalpha);
double beta = as<double>(LDAbeta);
double prec = as<double>(IRTprec);
Rcpp::NumericVector initIdeal(StartIdealPoints);
Rcpp::LogicalVector store(SaveChains);
int burn = as<int>(BurntIterations);
int mcmc = as<int>(SavedIterations);
int nvotes = Y.nrow(); 
int nvoters = Y.ncol(); 

Rprintf("%d voters\n",nvoters);
Rprintf("%d votes\n",nvotes);
Rprintf("%d terms\n",nterms);
Rprintf("%d total words\n",nwords);
    R_FlushConsole();
    R_ProcessEvents(); 

// SET UP CHAIN STORAGE OBJECTS

int wordchainstorageloc = 0;
int wordchainlength = nterms*ntopics;
if (store[0]) wordchainlength = wordchainlength*mcmc;
Rcpp::NumericVector wordchain(wordchainlength);
	
int topicchainstorageloc = 0;		
int topicchainlength = ntopics*nvotes;			
if (store[1]) topicchainlength = topicchainlength*mcmc;
Rcpp::NumericVector topicchain(topicchainlength);

int idealchainstorageloc = 0;	
int idealchainlength = nvoters*ntopics;
if (store[2]) idealchainlength = idealchainlength*mcmc;
Rcpp::NumericVector idealchain(idealchainlength);
		
int discrimchainstorageloc = 0;			
int discrimchainlength = nvotes*2;
if (store[3]) discrimchainlength = discrimchainlength*mcmc;
Rcpp::NumericVector discrimchain(discrimchainlength);

Rcpp::NumericVector rhochain(mcmc); // chain for rho
	
Rcpp::NumericVector llldachain(mcmc); // log-likelihood chain for LDA model
Rcpp::NumericVector llirtchain(mcmc); // log-likelihood chain for IRT model
	
// SET UP STATE VARIABLES	

Rcpp::NumericMatrix LDAnmz(nvotes,ntopics);						// working topic counts for each vote
Rcpp::IntegerVector LDAnm(nvotes);								// total number of terms in each document
Rcpp::NumericMatrix LDAnzt(ntopics,nterms);						// working term use by topic
Rcpp::IntegerVector LDAnz(ntopics);								// total number of term appearances for each topic
Rcpp::IntegerVector LDAzmn(nwords);								// working topic assignment for each word
Rcpp::NumericMatrix LDAphi(ntopics,nterms);						// multinomial parameters for term probabilities within topics

Rcpp::NumericMatrix LDAIRTlambda(nvotes,ntopics);				// working topic weights for each vote

Rcpp::NumericMatrix IRTystar(nvotes,nvoters);					// working latent utilities for each voter vote
Rcpp::NumericMatrix IRTtheta(nvoters,ntopics); 					// working voter ideal points for each topic
arma::mat IRTthetaarma(nvoters,ntopics);
Rcpp::NumericMatrix IRTbeta(nvotes,2);							// working discrimination parameters for each vote

// INITIALIZE WORKING VARIABLES	FOR LDA

double pvecdenom = 0;
double phidenom = 0;
double lambdadenom = 0;

double LDAll = 0;

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

// INITIALIZE TOPIC MODEL

Rprintf("Initializing Topic Model...\n");  
    R_FlushConsole();
    R_ProcessEvents();     

Rcpp::NumericVector pvec = rep(1.0/ntopics,ntopics);
int wordcounter = 0;
for (int dtc = 0; dtc < ndocterms; dtc++){  // loop through each sparse matrix entry...
	for (int l = 0; l < WMC(dtc); l++){  // for all appearances of the term...
		LDAzmn(wordcounter) = eliot::rcat(pvec); 
		LDAnmz(WMD(dtc),LDAzmn(wordcounter))++; 
		LDAnm(WMD(dtc))++; 
		LDAnzt(LDAzmn(wordcounter),WMT(dtc))++; 
		LDAnz(LDAzmn(wordcounter))++;
		wordcounter++;
	}
}			

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
	
	// Reassign	each term appearance to a topic
	
	LDAll = 0;
	wordcounter = 0;
	for (int dtc = 0; dtc < ndocterms; dtc++){
	
		for (int l = 0; l < WMC(dtc); l++){  // draw topic assignment for each appearances of the term...
		
			// unassign the current appearance of the term...
			LDAnmz(WMD(dtc),LDAzmn(wordcounter))--; 
			LDAnm(WMD(dtc))--; 
			LDAnzt(LDAzmn(wordcounter),WMT(dtc))--; 
			LDAnz(LDAzmn(wordcounter))--;
		
			for (int k = 0; k < ntopics; k++){
				pvec(k) = (LDAnmz(WMD(dtc),k) + alpha)*(LDAnzt(k,WMT(dtc)) + beta)/(LDAnz(k) + nterms*beta);
			}   
		
			// reassign the current appearance of the term...								
			LDAzmn(wordcounter) = eliot::rcat(pvec); 
			
			// save the current conditional log-likelihood of the term...	
			// This yields an approximation to the final iteration log-likelihood, 
			// but avoids doing all these loops twice for each iteration						
			LDAll = LDAll + log(pvec(LDAzmn(wordcounter))/sum(pvec)); 
	
			LDAnmz(WMD(dtc),LDAzmn(wordcounter))++; 
			LDAnm(WMD(dtc))++; 
			LDAnzt(LDAzmn(wordcounter),WMT(dtc))++; 
			LDAnz(LDAzmn(wordcounter))++;
			wordcounter++;
		}
	}


	// Calculate word mixtures for each topic
	
	for (int k = 0; k < ntopics; k++){
		for (int t = 0; t < nterms; t++){
			LDAphi(k,t) = (LDAnzt(k,t) + beta)/(LDAnz(k) + nterms*beta);
		}
	}
	
	
	// Caculate topic mixtures for each document
	
	for (int m = 0; m < nvotes; m++){
		for (int k = 0; k < ntopics; k++){
			LDAIRTlambda(m,k) = (LDAnmz(m,k) + alpha)/(LDAnm(m) + ntopics*alpha);
		}
	}
	
	
	// Draw truncated random normal for each vote for each voter

	for (int m = 0; m < nvotes; m++){ // for all votes...		
		for (int i = 0; i < nvoters; i++){ // for all voters...
			muim = -1.0*IRTbeta(m,0);
			for (int k=0; k < ntopics; k++) muim = muim + IRTbeta(m,1)*LDAIRTlambda(m,k)*IRTtheta(i,k);
			if (Y(m,i) == 1){ // yes vote observed
				// IRTystar(m,i) = eliot::rltnorm(muim,1.0,0.0);
				IRTystar(m,i) = RcppTN::rtn1(muim,1.0,0.0,INFINITY);				
			}
			if (Y(m,i) == 0){ // no vote observed
			  // IRTystar(m,i) = eliot::rutnorm(muim,1.0,0.0);
			  IRTystar(m,i) = RcppTN::rtn1(muim,1.0,-INFINITY,0.0);	
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
			for (int k=0; k < ntopics; k++) Xs(i,1) = Xs(i,1) + LDAIRTlambda(m,k)*IRTtheta(i,k);
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
			Bs(m,k) = LDAIRTlambda(m,k)*IRTbeta(m,1);	
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
			for (int k=0; k < ntopics; k++) muim = muim + IRTbeta(m,1)*LDAIRTlambda(m,k)*IRTtheta(i,k);
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
		
		if (store[0]){
			for (int t = 0; t < nterms; t++){
				for (int k = 0; k < ntopics; k++){
					wordchain(wordchainstorageloc) = LDAphi(k,t);
					wordchainstorageloc++;
				}
			}	
		} else {
			for (int t = 0; t < nterms; t++){
				for (int k = 0; k < ntopics; k++){
					wordchain(wordchainstorageloc) = wordchain(wordchainstorageloc) + LDAphi(k,t)/mcmc;
					wordchainstorageloc++;
				}
			}
			wordchainstorageloc = 0;
		}
			
		if (store[1]){
			for (int m = 0; m < nvotes; m++){
				for (int k = 0; k < ntopics; k++){
					topicchain(topicchainstorageloc) = LDAIRTlambda(m,k);
					topicchainstorageloc++;
				}
			}	
		} else {
			for (int m = 0; m < nvotes; m++){
				for (int k = 0; k < ntopics; k++){
					topicchain(topicchainstorageloc) = topicchain(topicchainstorageloc) + LDAIRTlambda(m,k)/mcmc;
					topicchainstorageloc++;
				}
			}
			topicchainstorageloc = 0;
		}
			
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
			
	llldachain(iter - burn) = LDAll;
	llirtchain(iter - burn) = IRTll;
	
	} 
    
    R_FlushConsole();
    R_ProcessEvents();
	
} // END MARKOV CHAIN MONTE CARLO

Rprintf("Markov Chain Monte Carlo Completed...\n");
	
// Load Mean Posterior Estimates into state variables for IRT model

topicchainstorageloc = 0;
if (store[1]){
	for (int m = 0; m < nvotes; m++){
		for (int k = 0; k < ntopics; k++){
			LDAIRTlambda(m,k) = 0;
		}
	}
	for (int iter = 0; iter < mcmc; iter++){
		for (int m = 0; m < nvotes; m++){
			for (int k = 0; k < ntopics; k++){
				LDAIRTlambda(m,k) = LDAIRTlambda(m,k) + topicchain(topicchainstorageloc)/mcmc;
				topicchainstorageloc++;
			}
		}	
	}
} else {
	for (int m = 0; m < nvotes; m++){
		for (int k = 0; k < ntopics; k++){
			LDAIRTlambda(m,k) = topicchain(topicchainstorageloc);
			topicchainstorageloc++;
		}
	}		
}

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
		for (int k=0; k < ntopics; k++) muim = muim + IRTbeta(m,1)*LDAIRTlambda(m,k)*IRTtheta(i,k);
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

return Rcpp::List::create(Rcpp::Named("word.chain") = wordchain, Rcpp::Named("topic.chain") = topicchain, Rcpp::Named("ideal.chain") = idealchain, Rcpp::Named("discrim.chain") = discrimchain, Rcpp::Named("rho.chain") = rhochain, Rcpp::Named("ll.lda.chain") = llldachain, Rcpp::Named("ll.irt.chain") = llirtchain, Rcpp::Named("dbar") = dbar,Rcpp::Named("pd") = pd,Rcpp::Named("dic") = dic, Rcpp::Named("accept.rho") = mcmcacceptrate);


}



