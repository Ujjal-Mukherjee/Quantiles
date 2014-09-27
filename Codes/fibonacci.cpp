#include <algorithm>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <R.h>
#include <Rmath.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double InnerProduct(NumericVector x){

	int n=x.size();
	double p=0;

	for(int i=0; i<n; i++)
		p+=x[i]*x[i];

	return p;
}


using namespace Rcpp;

// [[Rcpp::export]]

double VectMult(NumericVector x, NumericVector y){

	int nx=x.size();
	int ny=y.size();
	int m;

	if(nx!=ny){
		std::cout << "Non-conforming Arrays\n";
		exit;
	} 
	else {
		for(int i=0;i<nx;i++){
			m+=x[i]*y[i];
		}	
		return m;
	}
}

using namespace Rcpp;

// [[Rcpp::export]]

NumericMatrix MatrixMult(NumericMatrix x, NumericMatrix y){

	int nx=x.nrow();
	int kx=x.ncol();
	int ny=y.nrow();
	int ky=y.nrow();
	
	
	if(kx!=ny){
		exit;
	}
	else {
		Rcpp::NumericMatrix M(nx,ky);

		for(int i=0;i<nx;i++){
			for(int j=0;j<ky;j++){
				M(i,j)=VectMult(x(i,_),y(_,j));
			}
		}
		return M;
	}
}
