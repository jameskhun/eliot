// includes from the plugin
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::interfaces(cpp)]]

// [[Rcpp::export]]
	int rcat(Rcpp::NumericVector pvec){
		double total = sum(pvec);
		double udraw = Rf_runif(0,total);
		double cumsum = pvec[0];
		int draw = 0;
		while (cumsum < udraw) {
			cumsum = cumsum + pvec[draw+1];
			draw ++;
			}
		return draw;  
	}

// [[Rcpp::export]]	
	double rutnorm(double mean, double sd, double ub){
		double pi = M_PI;
	    int method;  
		double alphastar;
		NumericVector z(1);
		double draw;
		NumericVector acceptprob(1);
		NumericVector acceptdraw(1);
		
		// rescale to standard normal
		double muplus = (ub - mean)/sd;	
		
		if (muplus <= 0) {
			method = 1; 
		} else {
			method = 2;	
		}
			
		switch(method){
	 	
	 	case 1: // exponential proposal using upper bound
	 	
			alphastar = (-muplus + sqrt(muplus* muplus + 4))/2;
			do { 
				z = muplus - rexp(alphastar);
				acceptprob = exp(-1*(z[0]+alphastar)*(z[0]+alphastar)/2);
				acceptdraw = runif(1);
			} while (acceptdraw[0] > acceptprob[0]);
		    draw = z[0]*sd + mean;  
			return draw;
		 		
	 	case 2: // normal proposal	
		
			do { 
				z = rnorm(1);
	      	} while (z[0] > muplus);
	      	draw = z[0]*sd + mean; 
			return draw;
						
		} 
	}
	
// [[Rcpp::export]]	
	double rltnorm(double mean, double sd, double lb){
		double pi = M_PI;
	    int method;  
		double alphastar;
		NumericVector z(1);
		double draw;
		NumericVector acceptprob(1);
		NumericVector acceptdraw(1);
		
		// rescale to standard normal
		double muminus = (lb - mean)/sd;
	
		if (muminus >= 0) {
			method = 1; 
		} else {
			method = 2;	
		}
			
		switch(method){		
	
		case 1: // exponential proposal using lower bound
				
			alphastar = (muminus + sqrt(muminus*muminus + 4))/2;
			do { 
				z = muminus + rexp(alphastar);
				acceptprob = exp(-1*(z[0]-alphastar)*(z[0]-alphastar)/2);
				acceptdraw = runif(1);
			} while (acceptdraw[0] > acceptprob[0]);
		    draw = z[0]*sd + mean; 
			return draw;
	 	
	 	case 2: // normal proposal	
		
			do { 
				z = rnorm(1);
	      	} while (z[0] < muminus);
	      	draw = z[0]*sd + mean; 
			return draw;
						
		} 
	}
	
	
	
	
// [[Rcpp::export]]
	double rdtnorm(double mean, double sd, double lb, double ub){
		double pi = M_PI;
	    int method;  
		double alphastar;
		NumericVector z(1);
		double draw;
		NumericVector acceptprob(1);
		NumericVector acceptdraw(1);
		
		// rescale to standard normal
		double muminus = (lb - mean)/sd;
		double muplus = (ub - mean)/sd;	
		
		// uniform proposal method	
		
		double critp = ((2*sqrt(exp(1)))/(muminus + sqrt(muminus*muminus+4)))*exp((muminus*muminus - muminus*sqrt(muminus*muminus+4))/(4));
		double critm = ((2*sqrt(exp(1)))/(-muplus + sqrt(muplus* muplus +4)))*exp((muplus* muplus + muplus*sqrt(muplus* muplus +4))/(4));
		
		if (muplus*muminus < 0){
			if (muplus - muminus < sqrt(2*pi)) {
				method = 1; 
				} else {
				method = 4;
				}
			}
	
		if (muminus >= 0) {
			if (muplus > muminus + critp) {
				method = 2; 
			} else {
				method = 1;	
			}
		}
		
		if (muplus <= 0) {
			if (muplus > muminus + critm) {
				method = 3; 
			} else {
				method = 1;	
			}
		}
			
		switch(method){
		
		case 1: // uniform proposal
			
			do { 
				z = runif(1,muminus,muplus);
				if (muplus <= 0){
					acceptprob = exp((muplus*muplus - z[0]*z[0])/2);
					}	
				if (muminus >= 0){
					acceptprob = exp((muminus*muminus - z[0]*z[0])/2);
					}
				if (muplus > 0 && muminus < 0){
					acceptprob = exp(-z[0]*z[0]/2);
					}	
				acceptdraw = runif(1);	
			} while (acceptdraw[0] > acceptprob[0]);	
		    draw = z[0]*sd + mean; 
			return draw;			
	
		case 2: // exponential proposal using lower bound
				
			alphastar = (muminus + sqrt(muminus*muminus + 4))/2;
			do { 
				z = muminus + rexp(alphastar);
				acceptprob = exp(-1*(z[0]-alphastar)*(z[0]-alphastar)/2);
				acceptdraw = runif(1);
			} while (acceptdraw[0] > acceptprob[0] || z[0] > muplus);
		    draw = z[0]*sd + mean; 
			return draw;
	 	
	 	case 3: // exponential proposal using upper bound
	 	
			alphastar = (-muplus + sqrt(muplus* muplus + 4))/2;
			do { 
				z = muplus - rexp(alphastar);
				acceptprob = exp(-1*(z[0]+alphastar)*(z[0]+alphastar)/2);
				acceptdraw = runif(1);
			} while (acceptdraw[0] > acceptprob[0] || z[0] < muminus);
		    draw = z[0]*sd + mean;  
			return draw;
		 		
	 	case 4: // normal proposal	
		
			do { 
				z = rnorm(1);
	      	} while (z[0] < muminus || z[0] > muplus);
	      	draw = z[0]*sd + mean; 
			return draw;
						
		} 	
	}



