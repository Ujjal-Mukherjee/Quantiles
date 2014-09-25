/**************************************************************************************************
 *This program calculates the data depth for a point px in a data cloud dx                        *
 *------------------------------------------------------------------------                        *
 *The program contains the following function                                                     *
 *(1) IsInQuadrant(): Returns a logical value whether a point is in a specified hyper-quadrant    *
 *(2) QuadrantFrac(): Returns the fraction of total points in dx that are in a specified quadrant *
 *(3) ecdf(): Returns the multivariate empirical CDF for a point px given a data cloud dx         *
 *(4) ecdfGrid(): Returns a matrix containing ecdf values for a frid of points given input data dx*
 *(5) GenQuadIndicator(): Returns all possible combinations of hyper-quadrant indicator for a     *
 *                      multi-variate data matrix for a given point vector px                     *
 *(6) depthx(): Returns the depth of a point px given a multivariate data cloud dx                *
 *(7) medianp(): Returns the median for a multivariate data cloud dx                              *
 *(8) depthGrid(): Returns depth values over a grid of points given data clous dx                 *
 **************************************************************************************************/
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <R.h>
#include <Rmath.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>


/**************************************************************************************************
 *IsInQuadrant(): Returns a logical value whether a point is in a specified hyper-quadrant        *
 **************************************************************************************************/
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

int IsInQuadrant(NumericVector px, NumericVector xx, NumericVector ex){

	arma::colvec e(ex.begin(), ex.size(), false);
	arma::colvec x(xx.begin(), xx.size(), false);
	arma::colvec p(px.begin(), px.size(), false);
	
	arma::colvec e_sign = arma::floor(2*e-0.5);

	arma::colvec x_signed = e_sign % x;
	arma::colvec p_signed = e_sign % p;
	arma::colvec diff_p_x = p_signed-x_signed;
	
	float i = diff_p_x.min();
	
	if(i>=0)
	{
		return 1;
	}
	else
	{
		return 0;
	}

}

/**************************************************************************************************
 *QuadrantFrac(): Returns the fraction of total points in dx that are in a specified quadrant     *
 **************************************************************************************************/
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

float QuadrantFrac(NumericMatrix dx, NumericVector x, NumericVector e){

	arma::mat d(dx.begin(), dx.nrow(), dx.ncol(), false);
	unsigned int n = d.n_rows;
	int count = 0;

	for(unsigned int i=0; i<n; i++){
		count += IsInQuadrant(Rcpp::wrap(arma::strans(d.row(i))), x, e);
	}

	double p = (static_cast<float>(count))/(static_cast<float>(n));

	return p;

}

/**************************************************************************************************
 *ecdf(): Returns the multivariate empirical CDF for a point px given a data cloud dx             *
 **************************************************************************************************/

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

float ecdf(NumericMatrix dx,  NumericVector x){

	int n = x.size();

	Rcpp::NumericVector e(n,0.0);

	float p = QuadrantFrac(dx,x,e);

	return p;

}

/**************************************************************************************************
 *ecdfGrid(): Returns a matrix containing ecdf values for a frid of points given input data dx    *
 **************************************************************************************************/
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::colvec ecdfGrid(NumericMatrix dx, NumericMatrix xx){

	arma::mat x(xx.begin(), xx.nrow(), xx.ncol(), false);
	
	unsigned int n = x.n_rows;

	arma::colvec p(n);

	for(unsigned int i=0; i<n; i++){
		
		p(i) = ecdf(dx, Rcpp::wrap(arma::strans(x.row(i))));

	}

	return p;
}

/**************************************************************************************************
 *GenQuadIndicator(): Returns all possible combinations of hyper-quadrant indicator for a         *
 *                    multi-variate data matrix for a given point vector px                       *
 **************************************************************************************************/
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::NumericMatrix GenQuadIndicator(NumericVector x){

	int n = x.size();

	int s = std::pow(2,n);

	Rcpp::NumericMatrix e = Rcpp::wrap(arma::zeros(s,n));

	for(int i=0; i<s; i++){

		int temp = i;
		
		for(int j=0; j<n; j++){

			e(i,(n-1-j)) = temp % 2;
			temp = temp / 2;
	
		}	

	}

	return e;
}

/**************************************************************************************************
 *depthx(): Returns the depth of a point px given a multivariate data cloud dx                    *
 **************************************************************************************************/
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double depthx(NumericMatrix dx, NumericVector px){

	Rcpp::NumericMatrix ex = GenQuadIndicator(px);

	int n = ex.nrow();

	/*int nd = dx.nrow();*/

	arma::mat e(ex.begin(), ex.nrow(), ex.ncol(), false);

	Rcpp::NumericVector p = Rcpp::wrap(arma::zeros(n));

	for(int i=0; i<n; i++){
		
		Rcpp::NumericVector et = Rcpp::wrap(arma::strans(e.row(i)));
		p(i) = QuadrantFrac(dx,px,et);

	}

	Rcpp::NumericVector PX = p;
	arma::colvec P(PX.begin(), PX.size(), false);

	double varp = arma::var(P);

	double depth = exp(0.0-varp);

	return depth;
}

/**************************************************************************************************
 *medianp(): Returns the median for a multivariate data cloud dx                                  *
 **************************************************************************************************/
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector medianp(NumericMatrix dx){

	arma::mat d(dx.begin(), dx.nrow(), dx.ncol(), false);

	int n = d.n_rows;

	Rcpp::NumericVector dep(n);

	for(int i=0; i<n; i++){

		Rcpp::NumericVector vtemp = Rcpp::wrap(arma::strans(d.row(i)));

		dep(i) = depthx(dx, vtemp);

	}

	arma::colvec depx(dep.begin(), dep.size(), false);

	arma::uvec vindex = arma::sort_index(depx, "descend");

	int index = vindex[0];

	Rcpp::NumericVector med = Rcpp::wrap(arma::strans(d.row(index)));

	return med;

}

/**************************************************************************************************
 *depthGrid(): Returns depth values over a grid of points given data clous dx                     *
 **************************************************************************************************/
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector DepthGrid(NumericMatrix dx){

	arma::mat d(dx.begin(), dx.nrow(), dx.ncol(), false);

	int n = d.n_rows;

	Rcpp::NumericVector dep(n);

	for(int i=0; i<n; i++){

		Rcpp::NumericVector vtemp = Rcpp::wrap(arma::strans(d.row(i)));

		dep(i) = depthx(dx, vtemp);

	}

	return dep;
}


/**************************************************************************************************
 *CenterScale(): Returns depth values over a grid of points given data clous dx                     *
 **************************************************************************************************/

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericMatrix CenterScale(NumericMatrix dx, NumericVector mx, double S){

	arma::mat d(dx.begin(), dx.nrow(), dx.ncol(), false);

	arma::colvec m(mx.begin(), mx.size(), false);

	arma::mat d_m = d - arma::ones(d.n_rows) * m;

	d_m = d_m / S;

	return Rcpp::wrap(d_m);

}



/**************************************************************************************************
 *                                              END OF FILE                                       *
 **************************************************************************************************/












