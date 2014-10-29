/**********************************************************************************************************
 *This Program Computes Weighted Projection Quantile From Data Cloud.                                     *
 *The Program Has Several Function to Compute The Following.                                              *
 *<pq> Projection Quantile for a direction vector u.                                                      *
 *<genU> Generate cartesian vector of points from polar coordinates of size n.                            *
 *<ProjQuant> Generate a complete projection quantile profile from a cartesian vector. Calls <pq>.        *
 *<KMeanDist> Calculates the mean distance of the k max. dist points from a point.                        *
 *<OrthoProjVecNorm> Computes a vector of orthogonal norms for a data cloud along a direction vector.     *
 *<WProjQuant> Computes the weighted projection quantile along a direction u, given the k-mean dist,      *
 *orthogonal projection norms and a vector of tuning parameters.                                          *
 *<WtProjQuantProfile> Generates a weighted projection quantile profile along a polar coordinate vector.  *
 *NOTE: THIS IS A ALPHA VERSION. NEEDS MORE GENERALIZABILITY FOR THE FINAL VERSION.                       *
 **********************************************************************************************************/
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <R.h>
#include <Rmath.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>

#ifndef PI
#define PI 3.1415926535897932
#endif

/****************************************************
 *<pq> Projection Quantile for a direction vector u.*
 ****************************************************/

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::colvec pq(NumericMatrix dx, NumericVector ux){

	arma::colvec u(ux.begin(), ux.size(), false);
	arma::mat d(dx.begin(), dx.nrow(), dx.ncol(), false);

	double Normu = std::inner_product(u.begin(),u.end(),u.begin(),0.0);
	arma::colvec U = u/sqrt(Normu);

	arma::colvec P = arma::sort(d*U,"ascend");
	double alpha = (1+sqrt(Normu))/2;
	
	int pos = abs(alpha*P.size()-1);
	arma::colvec ProjQuant = U*P[pos];

	return ProjQuant;
}

/******************************************************************************
 *<genU> Generate cartesian vector of points from polar coordinates of size n.*
 ******************************************************************************/

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::Mat<double> genU(int n, NumericVector ux){

	arma::colvec u(ux.begin(), ux.size(), false);
	double Normu = sqrt(std::inner_product(u.begin(),u.end(),u.begin(),0.0));

	double theta = 0; 
	
	arma::Mat<double> umat(n,2);
	
	double incr = 2*PI/n;

	for(int i=0; i<n; i++){

	umat(i,0) = Normu*cos(theta);
	umat(i,1) = Normu*sin(theta);
	theta += incr;

	}	
	
	return umat;

}


/**************************************************************************************************
 *<ProjQuant> Generate a complete projection quantile profile from a cartesian vector. Calls <pq>.*
 **************************************************************************************************/
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]


arma::Mat<double> ProjQuant(NumericMatrix dx, NumericVector ux, int n){

	arma::mat um = genU(n, ux);

	arma::Mat<double> Quantiles(n,ux.size());

	for(int i=0; i<n; i++){


		ux[0] = um(i,0);
		ux[1] = um(i,1);

		arma::colvec temp = pq(dx, ux);
		Quantiles(i,0) = temp[0];
		Quantiles(i,1) = temp[1];

	}

	return Quantiles;
}


/**********************************************************************************
 *<KMeanDist> Calculates the mean distance of the k max. dist points from a point.*
 **********************************************************************************/

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double KMeanDist(NumericMatrix dx, NumericVector px, unsigned int k){

	arma::colvec p(px.begin(), px.size(), false);
	arma::mat d(dx.begin(), dx.nrow(), dx.ncol(), false);

	arma::Mat<double> one = arma::ones(dx.nrow(), 1);
	arma::Mat<double> PMat = one*(arma::strans(p));

	arma::Mat<double> DiffMat = arma::square(d-PMat);

	arma::colvec DiffVec = DiffMat*arma::ones(DiffMat.n_cols,1);

	arma::colvec Sqrt_DiffVec = arma::sqrt(DiffVec);

	arma::colvec Sorted_Sqrt_DiffVec = arma::sort(Sqrt_DiffVec, "descend");

	unsigned int n = Sorted_Sqrt_DiffVec.n_elem;

	double KMean = 0;

	if (k <= n)
	{
		Sorted_Sqrt_DiffVec = Sorted_Sqrt_DiffVec.subvec(0,(k-1));
		KMean = mean(Sorted_Sqrt_DiffVec);
	}	
	else 
	{
		KMean = mean(Sorted_Sqrt_DiffVec);
	}
	
	return KMean;
}

/*****************************************************************************************************
 *<OrthoProjVecNorm> Computes a vector of orthogonal norms for a data cloud along a direction vector.*
 *****************************************************************************************************/
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::colvec OrthoProjVecNorm(NumericMatrix dx, NumericVector ux){

	arma::colvec u(ux.begin(), ux.size(), false);
	arma::mat d(dx.begin(), dx.nrow(), dx.ncol(), false);

	double Normu = std::inner_product(u.begin(),u.end(),u.begin(),0.0);
	arma::colvec U = u/sqrt(Normu);
	arma::colvec ProjNorm = arma::square(d*U);

	arma::colvec NormVec = arma::square(d)*(arma::ones(d.n_cols,1));

	arma::colvec InvNormVec = arma::pow(NormVec,-1);

	arma::colvec SqrtInvNormVec = arma::sqrt(InvNormVec);
	
	arma::colvec OrthoNorm = arma::sqrt(NormVec-ProjNorm);

	arma::colvec OrthoNormRat = OrthoNorm % SqrtInvNormVec;	

	return OrthoNormRat;

}

/****************************************************************************************************
 *<WProjQuant> Computes the weighted projection quantile along a direction u, given the k-mean dist,*
 *orthogonal projection norms and a vector of tuning parameters.                                    *
 ****************************************************************************************************/

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::colvec WtProjQuant(NumericMatrix dx, NumericVector ux, NumericVector we, double di, float a, float b, float l){

	arma::colvec u(ux.begin(), ux.size(), false);
	arma::mat d(dx.begin(), dx.nrow(), dx.ncol(), false);
	arma::colvec w(we.begin(), we.size(), false);

	double Normu = std::inner_product(u.begin(),u.end(),u.begin(),0.0);
	arma::colvec U = u/sqrt(Normu);
	arma::colvec Proj = d*U;

	arma::colvec aa = arma::pow(w,b);

	arma::colvec bb = arma::ones(w.n_elem)*di;

	arma::colvec dd = -aa-bb;

	arma::colvec wx = arma::exp(dd);

	arma::colvec WtProj = Proj % wx;

	arma::Col<unsigned int> w_index = arma::sort_index(wx);

	arma::Col<unsigned int> w_index_sort = w_index.subvec(floor(l*(we.size())),(we.size()-1));

	arma::colvec wx_sort = wx.elem(w_index_sort);

	arma::colvec WtProj_sub = WtProj.elem(w_index_sort);

	arma::colvec Proj_sub = Proj.elem(w_index_sort);

	arma::umat WtProj_sub_sort = sort_index(WtProj_sub);

	int alpha = abs(floor((WtProj_sub.n_rows)*((1+sqrt(Normu))/2))-1);

	/*int alpha_index = WtProj_sub_sort(alpha);

	double WtProjQuantile = WtProj_sub(alpha_index)/wx_sort(alpha_index);

	 double WtProjQuantile = Proj_sub.elem(alpha_index);*/

	arma::colvec Proj_sub_sort = arma::sort(Proj_sub, "ascend");
	
	double WtProjQuantile = Proj_sub_sort(alpha);

	arma::colvec WtProjQuantileVec = WtProjQuantile*U;
	
	return WtProjQuantileVec;

}

/********************************************************************************************************
 *<WtProjQuantProfile> Generates a weighted projection quantile profile along a polar coordinate vector.*
 ********************************************************************************************************/

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::NumericVector WtProjQuantProfile(NumericMatrix dx, NumericVector ux, int n, int k, float a, float b, float l){

	arma::mat um = genU(n, ux);

	arma::Mat<double> Quantiles(n,ux.size());

	for(int i=0; i<n; i++){

		Rcpp::NumericVector u = Rcpp::wrap(arma::strans(um.row(i)));
		Rcpp::NumericVector q = Rcpp::wrap(pq(dx,u));
		double dk = KMeanDist(dx,q,k);
		Rcpp::NumericVector w = Rcpp::wrap(OrthoProjVecNorm(dx,u));
		Rcpp::NumericVector Quant = Rcpp::wrap(arma::strans(WtProjQuant(dx, u, w, dk, a, b, l)));
		Quantiles.row(i) = Rcpp::as<arma::rowvec>(Quant);

	}

	return Rcpp::wrap(Quantiles);


}

/********************************************************************************************************
 *<ProjQuantileDepth> Generates a weighted projection quantile profile along a polar coordinate vector.*
 ********************************************************************************************************/
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double ProjQuantileDepth(NumericMatrix dx, NumericVector ux, int k, float a, float b, float l){

	arma::colvec u(ux.begin(), ux.size(), false);
	arma::mat d(dx.begin(), dx.nrow(), dx.ncol(), false);

	double Normu = std::inner_product(u.begin(),u.end(),u.begin(),0.0);
	arma::colvec U = u/sqrt(Normu);

	if(Normu>=1.0)
	{
		ux = Rcpp::wrap(U);
	}
	
	Rcpp::NumericVector q = Rcpp::wrap(pq(dx,ux));
	double dk = KMeanDist(dx,q,k);
	Rcpp::NumericVector w = Rcpp::wrap(OrthoProjVecNorm(dx,ux));
	arma::colvec Q = WtProjQuant(dx, ux, w, dk, a, b, l);

	double norml = std::inner_product(Q.begin(), Q.end(), Q.begin(), 0.0);	

	double alpha = (1+sqrt(Normu))/2;

	double beta = alpha * Normu/norml;

	double Depth = exp(-beta);

	return Depth;
	
}


/********************************************************************************************************
 *<WtProjQuantMod> Generates a weighted projection quantile profile along a polar coordinate vector.    *
 ********************************************************************************************************/
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]


arma::colvec WtProjQuantMod(NumericMatrix dx, NumericVector ux, float k){



	arma::colvec u(ux.begin(), ux.size(), false);
	arma::mat d(dx.begin(), dx.nrow(), dx.ncol(), false);

	int alpha = k*dx.nrow();

	double Normu = std::inner_product(u.begin(),u.end(),u.begin(),0.0);
	arma::colvec U = u/sqrt(Normu);

	arma::colvec P = d*U;

	arma::mat Pu = P*arma::strans(U);

	arma::mat P_Ortho = d-Pu;

	arma::mat P_Ortho_Sq = arma::pow(P_Ortho,2);

	arma::colvec P_Ortho_Norm = P_Ortho_Sq*arma::ones(dx.ncol());

	arma::uvec P_Ortho_Sorted = arma::sort_index(P_Ortho_Norm,"ascend");

	arma::uvec P_Ortho_Sub_Vec = P_Ortho_Sorted.subvec(0,alpha);

	arma::mat d_submat = d.rows(P_Ortho_Sub_Vec);

	arma::colvec Proj_Quant = pq(Rcpp::wrap(d_submat),ux);
	
	return Proj_Quant;

}


/********************************************************************************************************
 *<WtProjQuantProfileMod> Generates a weighted projection quantile profile along a polar coordinate vector.*
 ********************************************************************************************************/

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::NumericVector WtProjQuantProfileMod(NumericMatrix dx, NumericVector ux, float k, int n){

	arma::mat um = genU(n, ux);

	arma::Mat<double> Quantiles(n,ux.size());

	for(int i=0; i<n; i++){

		Rcpp::NumericVector u = Rcpp::wrap(arma::strans(um.row(i)));
		Rcpp::NumericVector Quant = Rcpp::wrap(arma::strans(WtProjQuantMod(dx, u, k)));
		Quantiles.row(i) = Rcpp::as<arma::rowvec>(Quant);

	}

	return Rcpp::wrap(Quantiles);


}




















