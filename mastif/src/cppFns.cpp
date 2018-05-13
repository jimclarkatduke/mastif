#include <RcppArmadillo.h>
#include <Rmath.h>
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat kernYrRcpp(arma::mat dmat, arma::mat fec, arma::uvec years,
                     arma::uvec seedyear, arma::uvec treeyear,
                     arma::uvec seedrow, arma::uvec treecol){
  int ny = years.size();
  int nr = seedyear.size();
  int nf = fec.n_cols;
  int sindex, tindex, dsindex, dtindex;
  arma::mat lambda(nr,nf); lambda.fill(0);
  
  for(int j = 0; j < ny; j++){
    
    uvec ws = find(seedyear == years(j));
    uvec wt = find(treeyear == years(j));
    if(ws.size() == 0) continue;
    
    uvec ds = seedrow.elem( ws ) - 1;
    uvec dt = treecol.elem( wt ) - 1;
    
    for(int l = 0; l < nf; l++){
 
      for(unsigned int i = 0; i < ws.n_elem; i++){
      
        sindex = ws(i);
        dsindex = ds(i);
        double lsum = 0.0;
      
        for(unsigned int k = 0; k < wt.n_elem; k++){
        
          tindex  = wt(k);
          dtindex = dt(k);
      
          lsum = lsum + dmat(dsindex,dtindex)*fec(tindex,l);
        }
        lambda(sindex,l) = lsum;
      }
    }
  }
  return lambda;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List byRcpp(const int nr, const arma::mat frommat,
                  arma::mat totmat, arma::mat summat, 
                  arma::mat minmat, arma::mat maxmat){
  int i, j;
  double s;
  
  for(int k = 0; k < nr; k++){
    
    i = frommat(k,0) - 1;
    j = frommat(k,1) - 1;
    s = frommat(k,2);
    totmat(i,j) = totmat(i,j) + 1;
    summat(i,j) = summat(i,j) + s;
    
    if(s > maxmat(i,j))
      maxmat(i,j) = s;   
    
    if(s < minmat(i,j))
      minmat(i,j) = s;
  }
  
  return Rcpp::List::create(Rcpp::Named("total")=totmat,
                            Rcpp::Named("sum")=summat,
                            Rcpp::Named("min")=minmat,
                            Rcpp::Named("max")=maxmat);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double tnormRcpp(double lo, double hi, double mu, double sig){
  
  double q1, q2, z;
  
  q1 = Rf_pnorm5(lo,mu,sig,1,0);
  q2 = Rf_pnorm5(hi,mu,sig,1,0);
  z = Rf_runif(q1,q2);
  z = Rf_qnorm5(z, mu, sig, 1, 0);
  
  if(z > hi){
    z = lo;
  }
  
  if(z < lo){
    z = hi;
  }
  
  return(z);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat trMVNmatrixRcpp(arma::mat avec, arma::mat muvec, 
                          arma::mat smat, arma::mat lo,
                          arma::mat hi, arma::uvec whichSample, 
                          arma::uvec idxALL){
  int cindex;
  arma::rowvec av;
  arma::rowvec mv;
  arma::vec mAs(2);
  int nm = smat.n_rows;
  int nr = muvec.n_rows;
  arma::rowvec p1(nm-1);
  arma::mat sin(nm-1, nm-1);
  arma::uvec cid(1);
  arma::uvec idx;
  arma::mat m1(1,1);
  arma::mat s1(1,1);
  double tiny = min(smat.diag())*.0001;
  int nk = whichSample.n_elem;
  
  arma::mat A(nr, nm); A.fill(NA_REAL);
  arma::umat idxALLm(nm-1, nm);
  
  for(int j=0; j < nm; j++)
    
    idxALLm.col(j) = idxALL.elem( find(idxALL != j) );
  
  for(int i = 0; i < nr ; i++){
    
    for(int k = 0; k < nk; k++){
      
      cindex = whichSample[k]-1;
      
      av = avec.row(i);
      mv = muvec.row(i);
      
      cid(0) = cindex;
      idx = idxALLm.col(cindex);
      sin = arma::inv_sympd(smat.submat(idx, idx));
      p1 = trans(smat.submat(idx, cid)) * sin;
      
      m1 = mv[cindex] + dot(p1, (av.elem(idx) - mv.elem(idx)));
      s1 = smat(cindex,cindex) - dot(p1, smat.submat(cid, idx)) ;
      
      mAs[0] = m1(0,0);
      mAs[1] = s1(0,0);
      if(mAs[1] < 0) mAs[1] = tiny;  
      
      double sss = pow(mAs[1],.5);
      
      avec(i,cindex) = tnormRcpp(lo(i,cindex), hi(i,cindex), mAs[0], sss);
      A(i,cindex) = avec(i,cindex);
      
    }
  }
  return A;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat betaRcpp(int n, arma::mat X, arma::vec y, 
                   double sigma,  arma::mat AI){
  int ncols = AI.n_cols;
  
  arma::mat XX = X.t() * X;
  arma::colvec v = X.t()/sigma * y;
  arma::mat IXX = XX/sigma + AI;
  arma::mat V = inv_sympd(IXX);
  arma::colvec mu = V * v;
  
  arma::mat z = randn(n, ncols);
  
  return arma::repmat(mu, 1, n).t() + z * chol(V);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat randEffectRcpp(arma::uvec gindex, arma::uvec groups,
                         arma::mat X, arma::colvec y,
                         double sigma, arma::mat AI) {
  
  int ngroup = groups.n_elem;
  int q = X.n_cols;
  mat Z(ngroup,q); Z.fill(0);
  
  for(int j = 0; j < ngroup; j++){
    
    uvec ws = find(gindex == groups(j));
    int nj = ws.size();
    if(nj < 3) continue;
    
    mat tempW = X.rows(ws);
    vec tempY = y.elem(ws);
    
    Z.row(j) = betaRcpp(1, tempW, tempY, sigma, AI);
  }
  return Z;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat solveRcpp(arma::mat A) {
  arma::mat AA(A);
  arma::mat Ainv = arma::inv_sympd(AA);
  return Ainv;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat rmvnormRcpp(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  bool success = false;
  arma::mat S = sigma;
  
  arma::mat Y = randn(n, ncols);
  
  success = chol(S, sigma);
  if(success == false){
    sigma += eye(ncols,ncols) * 1e-5;
  }
  success = chol(S, sigma);
  if(success == false){
  //    throw std::range_error("sigma not positive definite");
      return arma::repmat(mu*0, 1, n).t();
  }
  return arma::repmat(mu, 1, n).t() + Y * chol(sigma);
}

