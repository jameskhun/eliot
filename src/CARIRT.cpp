#include <RcppArmadillo.h>
#include <RcppTN.h>
#include <eliot.h>

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List CARIRTcpp( SEXP ConnectMatrix, SEXP VoteMatrix, SEXP nonUnanimous, SEXP InitialAveragePrefs, SEXP InitialAlphas, SEXP InitialBetas,  SEXP SaveChains, SEXP BurntIterations, SEXP SavedIterations, SEXP ThinningFactor){

Rcpp::RNGScope scope;

// LOAD DATA
	
arma::mat A = Rcpp::as<arma::mat>(ConnectMatrix);
Rcpp::IntegerMatrix Y(VoteMatrix);  
int nvotes = Y.nrow(); 
int nvoters = Y.ncol();  
Rcpp::IntegerVector NonUnanimous(nonUnanimous);
Rcpp::LogicalVector store(SaveChains);
int burn = as<int>(BurntIterations);
int mcmc = as<int>(SavedIterations);
int thin = as<int>(ThinningFactor);

// SET UP WORKING VARIABLES
	
Rcpp::NumericMatrix mu(nvotes,nvoters);
arma::mat psi(nvotes,nvoters);
Rcpp::NumericVector theta(InitialAveragePrefs);
Rcpp::NumericVector alpha(InitialAlphas);
Rcpp::IntegerVector beta(InitialBetas);
Rcpp::NumericVector lambda(1);
Rcpp::NumericVector ll(1);
arma::vec llj(nvotes);

double alphadagger;
double alphaproposaltempsum;
double alphaproposaltempcount;
double betadagger;
Rcpp::NumericVector psidagger(nvoters);

double psilb;
double psiub;
Rcpp::NumericMatrix psiprec(nvotes,nvoters);
double psinum;
double alphamean;
double alphasd;
double MHratiolognum;
double MHratiologdenom;
double MHdraw;
Rcpp::NumericVector lambdatemp(1);
arma::mat psinumtemp(1,1);
arma::mat bpsitemp(1,1);
Rcpp::NumericVector Ajplus(nvotes);
arma::mat Wmat(nvotes,nvotes);
for (int j = 0; j < nvotes; j++){ // for all votes...
	Ajplus(j) = 0;
	for (int jj = 0; jj < nvotes; jj++){ // for all votes...
		Ajplus(j) = Ajplus(j) + A(j,jj);
		Wmat(j,jj) = - A(j,jj);
	}
	Wmat(j,j) = Ajplus(j);
}

double apsi;
double bpsi;
double apsiprior = 0;
double bpsiprior = 0; 

Rprintf("%d voters\n",nvoters);
Rprintf("%d votes\n",nvotes);
    R_FlushConsole();
    R_ProcessEvents(); 
      
// SET UP CHAIN STORAGE OBJECTS

int thincounter = 1;

int alphachainstorageloc = 0;	
int alphachainlength = nvotes;
if (store[0]) alphachainlength = alphachainlength*mcmc;
Rcpp::NumericVector alphachain(alphachainlength);

int betachainstorageloc = 0;	
int betachainlength = nvotes;
if (store[1]) betachainlength = betachainlength*mcmc;
Rcpp:: NumericVector betachain(betachainlength);

int lambdachainstorageloc = 0;	
int lambdachainlength = 1;
if (store[2]) lambdachainlength = lambdachainlength*mcmc;
Rcpp::NumericVector lambdachain(lambdachainlength);

int averageprefchainstorageloc = 0;	
int averageprefchainlength = nvoters;
if (store[3]) averageprefchainlength = averageprefchainlength*mcmc;
Rcpp::NumericVector averageprefchain(averageprefchainlength);

int expectedprefchainstorageloc = 0;	
int expectedprefchainlength = nvotes*nvoters;
if (store[4]) expectedprefchainlength = expectedprefchainlength*mcmc;
Rcpp::NumericVector expectedprefchain(expectedprefchainlength);

int latentprefchainstorageloc = 0;	
int latentprefchainlength = nvotes*nvoters;
if (store[5]) latentprefchainlength = latentprefchainlength*mcmc;
Rcpp::NumericVector latentprefchain(latentprefchainlength);
	
int llchainstorageloc = 0;		
Rcpp::NumericVector llchain(mcmc);      

// GENERATE INITIAL VALUES

Rprintf("Initializing Model...\n");
    R_FlushConsole();
    R_ProcessEvents(); 

// Draw initial psis from truncated standard normal, conditional on initial alpha and beta

for (int j = 0; j < nvotes; j++){ // for all votes...
	for (int i = 0; i < nvoters; i++){ // for all voters...
		if (beta[j] == 1){
			if (Y(j,i) == 1){ // yes vote observed
				// psi(j,i) = eliot::rltnorm(theta(i),0.5,alpha(j));
				psi(j,i)  = RcppTN::rtn1(theta(i),0.5,alpha(j),INFINITY);	
			}
			if (Y(j,i) == -1){ // no vote observed
				// psi(j,i) = eliot::rutnorm(theta(i),0.5,alpha(j));
				psi(j,i)  = RcppTN::rtn1(theta(i),0.5,-INFINITY,alpha(j));	
			}	
			if (Y(j,i) == 0){ // vote missing
				psi(j,i) = Rf_rnorm(theta(i),0.5);
			}		
		} else {
			if (Y(j,i) == 1){ // yes vote observed
				// psi(j,i) = eliot::rutnorm(theta(i),0.5,alpha(j));
				psi(j,i)  = RcppTN::rtn1(theta(i),0.5,-INFINITY,alpha(j));		
			}
			if (Y(j,i) == -1){ // no vote observed
				// psi(j,i) = eliot::rltnorm(theta(i),0.5,alpha(j));
			  psi(j,i)  = RcppTN::rtn1(theta(i),0.5,alpha(j),INFINITY);	
			}	
			if (Y(j,i) == 0){ // vote missing
				psi(j,i) = Rf_rnorm(theta(i),0.5);
			}
			
		}		
	}
}


// Set initial lambdas
lambda[0] = 1; 

Rprintf("Initiating Markov Chain Monte Carlo...\n");    
    R_FlushConsole();
    R_ProcessEvents();   
    

// BEGIN MARKOV CHAIN MONTE CARLO

for (int iter = 0; iter < (burn + mcmc*thin); iter++){
 
	// BEGIN LOOP OVER VOTES
	for (int j = 0; j < nvotes; j++){  
		
		// IF ANY VOTES IN CASE FOR ANY VOTER
		if (NonUnanimous(j) != -1){
					
			// PROPOSE BETA FLIP, ALPHA REFLECTION OVER MEAN LATENT VOTE FOR M-H STEP
			
			betadagger = -beta(j);
			
			alphaproposaltempsum = 0;
			alphaproposaltempcount = 0;
			for (int i = 0; i < nvoters; i++) {
				if (Y(j,i) != 0) {
					alphaproposaltempsum = alphaproposaltempsum + psi(j,i);
					alphaproposaltempcount++;
				}
			}
			alphadagger = 2.0*(alphaproposaltempsum/alphaproposaltempcount) - alpha(j);
					
			MHratiolognum = -0.5*alphadagger*alphadagger; // alpha prior
			MHratiologdenom = -0.5*alpha(j)*alpha(j); // alpha prior
			for (int i = 0; i < nvoters; i++) {
				psiprec(j,i) = lambda(0) * Ajplus(j);  // can this go outside the loop?
				psinumtemp = A.row(j)*psi.col(i);
				psinum = lambda(0)*psinumtemp(0,0);
				mu(j,i) = psinum / psiprec(j,i);
				if (Y(j,i) != 0){
					MHratiolognum = MHratiolognum + R::pnorm(alphadagger, mu(j,i), sqrt(1.0 / psiprec(j,i)), (Y(j,i)*(-betadagger)+1)/2, 1);
					MHratiologdenom = MHratiologdenom + R::pnorm(alpha(j), mu(j,i), sqrt(1.0 / psiprec(j,i)), (Y(j,i)*(-beta(j))+1)/2, 1);
				}
			}	
										
			MHdraw = Rf_runif(0.0,1.0);
			if (MHdraw < exp(MHratiolognum - MHratiologdenom)){ 
			 	alpha(j) = alphadagger;
				beta(j) = betadagger;
			} 
				
				
			// FOR VOTE, DRAW PSI FOR EACH VOTER	
			
			for (int i = 0; i < nvoters; i++){
				if (Y(j,i)*beta(j) == 1){
					// psi(j,i) = eliot::rltnorm(mu(j,i), sqrt(1.0 / psiprec(j,i)), alpha(j));
					psi(j,i) = RcppTN::rtn1(mu(j,i), sqrt(1.0 / psiprec(j,i)), alpha(j),INFINITY);
				}
				if (Y(j,i)*beta(j) == -1){
					// psi(j,i) = eliot::rutnorm(mu(j,i), sqrt(1.0 / psiprec(j,i)), alpha(j));
				  psi(j,i) = RcppTN::rtn1(mu(j,i), sqrt(1.0 / psiprec(j,i)), -INFINITY,alpha(j));
				}
				if (Y(j,i) == 0){
					psi(j,i) = Rf_rnorm(mu(j,i), sqrt(1.0 / psiprec(j,i)));
				}
			}									
			
			// FOR VOTE, DRAW ALPHA
			
			if (NonUnanimous(j) == 1){
			  psilb = -INFINITY;
			  psiub = INFINITY;
				if (beta(j) == 1){
					for (int i = 0; i < nvoters; i++){
						if (psi(j,i) > psilb && Y(j,i) == -1) psilb = psi(j,i);
						if (psi(j,i) < psiub && Y(j,i) == 1) psiub = psi(j,i);
					}	
				} else {
					for (int i = 0; i < nvoters; i++){
						if (psi(j,i) > psilb && Y(j,i) == 1) psilb = psi(j,i);
						if (psi(j,i) < psiub && Y(j,i) == -1) psiub = psi(j,i);
					}
				}	
				// alpha(j) = eliot::rdtnorm(0.0, 1.0, psilb, psiub);
				alpha(j) = RcppTN::rtn1(0.0, 1.0, psilb, psiub);
			} 
			
			if (NonUnanimous(j) == 0){
				if (beta(j) == 1){
					psiub = 1000.0;
					for (int i = 0; i < nvoters; i++){
						if (psi(j,i) < psiub && Y(j,i) == 1) psiub = psi(j,i);
					}	
					alpha(j) = eliot::rutnorm(0.0, 1.0, psiub);
				} else {
					psilb = -1000.0;
					for (int i = 0; i < nvoters; i++){
						if (psi(j,i) > psilb && Y(j,i) == 1) psilb = psi(j,i);
					}	
					alpha(j) = eliot::rltnorm(0.0, 1.0, psilb);
				}	
			}	
		}
		
		// IF NO VOTES IN CASE FOR ANY VOTER
		if (NonUnanimous(j) == -1){
			
			MHdraw = Rf_runif(0.0,1.0);		
			if (MHdraw > 0.5){
				beta(j) = -beta(j);
			}
			
			for (int i = 0; i < nvoters; i++){
				psiprec(j,i) = lambda(0) * Ajplus(j);
				psinumtemp = A.row(j)*psi.col(i);
				psinum = lambda(0)*psinumtemp(0,0);
				mu(j,i) = psinum / psiprec(j,i);
				psi(j,i) = Rf_rnorm(mu(j,i), sqrt(1.0 / psiprec(j,i)));
			}
					
			alpha(j) = Rf_rnorm(0.0, 1.0);	
		}		
					
	}	// END LOOP OVER VOTES
	
		
	// DRAW LAMBDA_PSI
	
	apsi = apsiprior + nvotes*nvoters/2.0;
	bpsi = bpsiprior;
	for (int i = 0; i < nvoters; i++){   
		bpsitemp = trans(psi.col(i))*Wmat*psi.col(i);
		bpsi = bpsi + 0.5*bpsitemp(0,0);	
	}
	lambda(0) = Rf_rgamma(apsi,1.0/bpsi);
	
	
	// COMPUTE LOG-LIKELIHOOD	
	
	ll = 0;
	for (int j = 0; j < nvotes; j++){ // for all votes...
		if (NonUnanimous(j) != -1){
			for (int i = 0; i < nvoters; i++){ // for all voters...
				if (Y(j,i) != 0) ll = ll + R::pnorm(alpha(j), mu(j,i), sqrt(1.0 / psiprec(j,i)), (Y(j,i)*(-beta(j))+1)/2, 1);
			}
		}
	}	
	
	// RENORMALIZE TO HAVE ALPHA MEAN ZERO, SD ONE
	
	alphamean = 0;
	for (int j = 0; j < nvotes; j++) alphamean = alphamean + alpha(j);
	alphamean = alphamean/nvotes;	
	
	alphasd = 0;	
	for (int j = 0; j < nvotes; j++) alphasd = alphasd + (alpha(j) - alphamean)*(alpha(j) - alphamean);
	alphasd = sqrt(alphasd/nvotes);
	
	for (int j = 0; j < nvotes; j++){
		alpha(j) = (alpha(j) - alphamean)/alphasd;
		for (int i = 0; i < nvoters; i++){  
			psi(j,i) = (psi(j,i) - alphamean)/alphasd;
			mu(j,i) = (mu(j,i) - alphamean)/alphasd;
		}
	}	
	lambda(0) = lambda(0)* alphasd* alphasd;

	// SAVE CURRENT VALUES OF THE CHAIN
	
	if (iter >= burn){
		
		if (thincounter == thin){ // if storing this iteration...
			thincounter = 0;
		
			if (store[0]){
				for (int j = 0; j < nvotes; j++){
					alphachain(alphachainstorageloc) = alpha(j);
					alphachainstorageloc++;
				}	
			} else {
				for (int j = 0; j < nvotes; j++){
					alphachain(alphachainstorageloc) = alphachain(alphachainstorageloc) + alpha(j)/mcmc;
					alphachainstorageloc++;
				}
				alphachainstorageloc = 0;
			}
				
			if (store[1]){
				for (int j = 0; j < nvotes; j++){
					betachain(betachainstorageloc) = double(beta(j));
					betachainstorageloc++;
				}	
			} else {
				for (int j = 0; j < nvotes; j++){
					betachain(betachainstorageloc) = betachain(betachainstorageloc) + double(beta(j))/mcmc;
					betachainstorageloc++;
				}
				betachainstorageloc = 0;
			}
				
			if (store[2]){
				for (int k = 0; k < 1; k++){
					lambdachain(lambdachainstorageloc) = lambda(k);
					lambdachainstorageloc++;
				}	
			} else {
				for (int k = 0; k < 1; k++){
					lambdachain(lambdachainstorageloc) = lambdachain(lambdachainstorageloc) + lambda(k)/mcmc;
					lambdachainstorageloc++;
				}
				lambdachainstorageloc = 0;
			}
			
			if (store[3]){
				for (int i = 0; i < nvoters; i++){
					averageprefchain(averageprefchainstorageloc) = sum(psi.col(i)) / nvotes;
					averageprefchainstorageloc++;
				}	
			} else {
				for (int i = 0; i < nvoters; i++){
					averageprefchain(averageprefchainstorageloc) = averageprefchain(averageprefchainstorageloc) + (sum(psi.col(i)) / nvotes)/mcmc;
					averageprefchainstorageloc++;
				}
				averageprefchainstorageloc = 0;
			}
			
			if (store[4]){
				for (int i = 0; i < nvoters; i++){
					for (int j = 0; j < nvotes; j++){
						expectedprefchain(expectedprefchainstorageloc) = mu(j,i);
						expectedprefchainstorageloc++;
					}
				}	
			} else {
				for (int i = 0; i < nvoters; i++){
					for (int j = 0; j < nvotes; j++){
						expectedprefchain(expectedprefchainstorageloc) = expectedprefchain(expectedprefchainstorageloc) + mu(j,i)/mcmc;
						expectedprefchainstorageloc++;
					}
				}
				expectedprefchainstorageloc = 0;
			}
			
			
			if (store[5]){
				for (int i = 0; i < nvoters; i++){
					for (int j = 0; j < nvotes; j++){
						latentprefchain(latentprefchainstorageloc) = psi(j,i);
						latentprefchainstorageloc++;
					}
				}	
			} else {
				for (int i = 0; i < nvoters; i++){
					for (int j = 0; j < nvotes; j++){
						latentprefchain(latentprefchainstorageloc) = latentprefchain(latentprefchainstorageloc) + psi(j,i)/mcmc;
						latentprefchainstorageloc++;
					}
				}
				latentprefchainstorageloc = 0;
			}
							
			llchain(llchainstorageloc) = ll(0);
			llchainstorageloc++;
			
		} // end mcmc storage loop
		
		thincounter++;	
	
	}
 
	Rprintf("Finished iteration %d of %d (log likelihood %f)... \n",iter+1,burn+mcmc*thin,ll(0)); // print simulation status 
	
	R_FlushConsole();
    R_ProcessEvents();
	
} // END MARKOV CHAIN MONTE CARLO

Rprintf("Markov Chain Monte Carlo Completed...\n");
	
return Rcpp::List::create(Rcpp::Named("alpha.chain") = alphachain, Rcpp::Named("beta.chain") = betachain, Rcpp::Named("lambda.chain") = lambdachain,Rcpp::Named("average.pref.chain") = averageprefchain,Rcpp::Named("expected.pref.chain") = expectedprefchain,Rcpp::Named("latent.pref.chain") = latentprefchain,  Rcpp::Named("ll.chain") = llchain);


}



