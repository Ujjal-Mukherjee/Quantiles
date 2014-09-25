library(inline)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)

gdf <- "
 #include <algorithm>
 #include <vector>
 #include <stdexcept>
 #include <cmath>
 #include <iostream>
 #include <R.h>
 #include <Rmath.h>
 #include <Rcpp.h>

 Rcpp::NumericVector xa(a);
 Rcpp::NumericVector xb(b);
 int n_xa = xa.size(), n_xb = xb.size();

 Rcpp::NumericVector xab(n_xa + n_xb - 1);
 for (int i = 0; i < n_xa; i++)
 for (int j = 0; j < n_xb; j++)
 xab[i + j] += xa[i] * xb[j];
 return xab;
 "
sink("source.txt")
fun <- cxxfunction(signature(a="numeric", b="numeric"),
 gdf, plugin="Rcpp", verbose=TRUE)
sink();

s<-scan("source.txt", what=character(), sep="\n");
s<-s[6:49];
for(i in 1:length(s)){
	if(i==1)s1<-strsplit(s[i], split=" :")[[1]][2];
	if(i>1)s1<-c(s1,strsplit(s[i], split=" :")[[1]][2]);
}

write.table(s1, "source1.cpp", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\n");


g1<-"double innerprod(std::vector<double> xa){
  int n = xa.size();
  double result=0;  // create output vector
  double ss = 0;

  for(int i = 0; i < n; ++i) {
    ss =ss + xa[i]*xa[i];
  }

  result = ss;
  return result;
}"


g2<-"Rcpp::NumericVector xs(s);
	std::vector<double> x = xs;
	double y = innerprod(xs);
	double z = sqrt(y);
	return wrap(z);
"

a2<-cxxfunction(signature(s="numeric"),incl=g1,
 body=g2, plugin="Rcpp", verbose=TRUE);

##############


a1<-"#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double F1(NumericVector a) {
  int n = a.size();
  double result=0;  // create output vector
  double ss = 0;

  for(int i = 0; i < n; ++i) {
    ss += pow(a[i],2);
  }

  result = ss;
  return result;
}"


a2<-"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

 Rcpp::NumericVector a(x);

  double result=0;

  // grab the R function F1
  Function F1( \"F1\" ) ; 
  result = as<double>( F1(a) );

  return result;

"


a3<-cxxfunction(signature(x="numeric"),incl=a1,
 body=a2, plugin="Rcpp", verbose=TRUE);

#####################################################################


	arma::Mat<double> Quantiles(dx.nrow(),ux.size());

	for(int i=0; i<n; i++){


		ux[0] = um(i,0);
		ux[1] = um(i,1);

		arma::colvec temp = pq(dx, ux);
		Quantiles(i,0) = temp[0];
		Quantiles(i,1) = temp[1];

	}

######################################################################


src <-
"       arma::mat x = arma::randu(5, 5);
       arma::mat y = x;
       arma::uvec indice;
       indice << 0 << 3;
       int nRm = indice.n_elem;
       for(int j = 0; j < nRm; j++)
               y.shed_row( indice(j) - j );
       return List::create(Named(\"x\")=x, Named(\"y\")=y, Named(\"ind\")=indice, Named(\"nRm\")=nRm);"
colX <- cxxfunction(body=src,plugin="RcppArmadillo")
