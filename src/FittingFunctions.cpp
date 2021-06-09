// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>
#include <iostream>
#include <vector>
// Returns sign of a vector
vec sgn(vec val) {
  vec out(val.size());
  for (int i=0; i<val.size(); i++)
  {
    if (val(i)>0) out(i)=1;
    else if  (val(i)<0) out(i)=-1;
    else out(i)=0;
  }
  return out;
}

// Make one step for the t-group of a gene (soft threshold)
vec glmnetSimple(mat X,vec Y,double lam1)
{
  int ng = Y.size();
  vec beta0(ng); beta0 = -solve( X, Y );
  vec sign0(ng); sign0 = sgn(beta0);
  vec beta(ng); beta = -solve( X, Y + sign0*lam1);
  vec sign(ng); sign = sgn(beta);
  for (int g=0;g<ng;g++) if (sign0(g)!=sign(g)) return(beta*0);
  return(beta);
}



// Calculates first (b) and second (a) derivatives of the function to minimize 
void GetHessian(mat x, vec beta, vec p, vec y, double lam2, vec w, vec& b, mat& a) 
{
  int ng = beta.size();
  int ns = y.size();
  a.fill(0.0);
  b.fill(0.0);
  
  // Calc with lambda1=0
  for (int g1=0; g1<ng; g1++)
  {
    if (g1>0) {b(g1) += 2*lam2*(beta(g1)-beta(g1-1)); a(g1,g1) += 2*lam2;}
    if (g1<(ng-1)) {b(g1) += 2*lam2*(beta(g1)-beta(g1+1)); a(g1,g1) += 2*lam2;}
    for (int s=0; s<ns; s++)
    {
      b(g1) += -x(s,g1)*(y(s)*(1.0-p(s))-(1-y(s))*p(s))*w(s);
      a(g1,g1) += x(s,g1)*x(s,g1)*p(s)*(1.0-p(s))*w(s);    
    }
    
    for (int g2=g1+1; g2<ng; g2++)
    {
      if (g2==g1-1) a(g1,g2) += -2*lam2;
      if (g2==g1+1) a(g1,g2) += -2*lam2;
      for (int s=0; s<ns; s++)
      {
        a(g1,g2) += x(s,g1)*x(s,g2)*p(s)*(1.0-p(s))*w(s);  
      }
    }
    for (int g2=0; g2<g1; g2++) a(g1,g2) = a(g2,g1);
  }
  return;
}

// Calculates intercept in iterative way 
double CalculateDeltaIntercept(vec y, vec p, vec w)
{
  int ns = y.size();
  double DeltaBeta0=1.0;
  double DeltaBeta0Prev=50.0;
  double Sum1=0;
  double Sum2=0;
  while(abs(DeltaBeta0-DeltaBeta0Prev)/(1e-5+abs(DeltaBeta0)+abs(DeltaBeta0Prev))>1.0e-5)
  {
    DeltaBeta0Prev=DeltaBeta0;
    for (int s=0; s<ns; s++) {Sum1+=w(s)*y(s); Sum2+=w(s)/(1.0/p(s)-1.0+DeltaBeta0); }
    DeltaBeta0 = Sum1/Sum2;
  }
  return(log(DeltaBeta0));
}



vec GroupRound(mat x0, vec y, vec tV, double lam1, double lam2, vec beta, double& Intercept, vec w, vec IndFor0,vec IndTFor0, vec& M, double& LLmin)
{
  int m = beta.size();
  int nt = tV.size();
  int ns = y.size();
  int m0 = m/nt;
  uvec dGg(nt);
  mat x0G(ns,nt);
  vec b(nt);
  mat a(nt,nt);
  vec p(ns); 
  vec betaOld(nt);
  vec betaNew(nt);
  double LL;
  
  p = 1.0/(1.0+exp(-M)); for (int s=0;s<ns;s++) {if (p(s) > 1.0-1.0e-5) p(s)=1.0-1.0e-5; else if (p(s)<1.0e-5) p(s)=1.0e-5;}
  double F=0; for (int s=0;s<ns;s++) F += -y(s)*log(p(s))*w(s) - (1-y(s))*log(1.0-p(s))*w(s);
  double dIntercept = CalculateDeltaIntercept(y, p, w);
  Intercept = Intercept + dIntercept;
  M = M + dIntercept;
  p = 1.0/(1.0+exp(-M)); for (int s=0;s<ns;s++) {if (p(s) > 1.0-1e-5) p(s)=1.0-1e-5; else if (p(s)<1.0e-5) p(s)=1.0e-5;}
  for (int s=0;s<ns;s++) LLmin += -y(s)*log(p(s))*w(s) - (1-y(s))*log(1.0-p(s))*w(s);
  LLmin += - F;
  // return(beta);
  for (int g=0;g<m0;g++)
  {
    for (int it=0;it<nt;it++) {dGg(it) = it*m0 + g;}
    betaOld = beta(dGg);
    x0G.fill(0);
    
    for (int s=0;s<ns;s++) {x0G(s,IndTFor0(s)) = x0(IndFor0(s),g);}
    GetHessian(x0G, betaOld, p, y, lam2,w,b,a);
    betaNew = glmnetSimple(a,b - a * betaOld,lam1);
    if (any(betaNew!=betaOld))
    {
      vec Mnew = M + x0G * (betaNew-betaOld);
      vec pnew = 1.0/(1.0+exp(-Mnew)); for (int s=0;s<ns;s++) {if (pnew(s) > 1.0-1e-5) pnew(s)=1.0-1e-5; else if (pnew(s)<1.0e-5) pnew(s)=1.0e-5;}
      beta(dGg) = betaNew;
      LL=0;
      for (int s=0;s<ns;s++) LL += -y(s)*log(pnew(s))*w(s) - (1-y(s))*log(1.0-pnew(s))*w(s);
      for (int it=0;it<nt;it++) LL += lam1*abs(betaNew(it))-lam1*abs(betaOld(it));
      for (int it=1;it<nt;it++) LL += lam2*(betaNew(it)-betaNew(it-1))*(betaNew(it)-betaNew(it-1))-lam2*(betaOld(it)-betaOld(it-1))*(betaOld(it)-betaOld(it-1));
      for (int it=0;it<(nt-1);it++) LL += lam2*(betaNew(it)-betaNew(it+1))*(betaNew(it)-betaNew(it+1))-lam2*(betaOld(it)-betaOld(it+1))*(betaOld(it)-betaOld(it+1));
      if (LL<LLmin)
      {
        LLmin = LL;
        p = pnew;
        M = Mnew;
      }
      else
      {
        beta(dGg) = betaOld;
      }
    }
  }
  return(beta);
}

vec SingleGeneRound(mat x0, vec y, vec tV, double lam1, double lam2, vec beta, double& Intercept, vec w, vec IndFor0,vec IndTFor0, vec& M, double& LLmin)
{
  int m = beta.size();
  int nt = tV.size();
  int ns = y.size();
  int m0 = m/nt;
  vec b(nt);
  mat a(nt,nt);
  vec p(ns); 
  uvec dGg(nt);
  mat x0G(ns,nt);
  double LL;
  p = 1.0/(1.0+exp(-M)); for (int s=0;s<ns;s++) {if (p(s) > 1.0-1e-5) p(s)=1.0-1e-5; else if (p(s)<1.0e-5) p(s)=1.0e-5;}
  double F=0; for (int s=0;s<ns;s++) F += -y(s)*log(p(s))*w(s) - (1-y(s))*log(1.0-p(s))*w(s);
  double dIntercept = CalculateDeltaIntercept(y, p, w);
  Intercept = Intercept + dIntercept;
  M = M + dIntercept;
  p = 1.0/(1.0+exp(-M)); for (int s=0;s<ns;s++) {if (p(s) > 1.0-1e-5) p(s)=1.0-1e-5; else if (p(s)<1.0e-5) p(s)=1.0e-5;}
  for (int s=0;s<ns;s++) LLmin += -y(s)*log(p(s))*w(s) - (1-y(s))*log(1.0-p(s))*w(s);
  LLmin += - F;
  for (int g=0;g<m0;g++)
  {
    for (int it=0;it<nt;it++) {dGg(it) = it*m0 + g;}
    x0G.fill(0);
    for (int s=0;s<ns;s++) {x0G(s,IndTFor0(s)) = x0(IndFor0(s),g);}
    GetHessian(x0G, beta(dGg), p, y, lam2,w,b,a);

    for (int it=0;it<nt;it++)
    {
      double betaOld = beta(g+m0*it);
      double betaTry = (b(it) - a(it,it) * betaOld)/a(it,it);
      int sigma0 = (betaTry > 0 ? 1 : (betaTry < 0 ? -1 : 0 ) );
      double betaNew = (b(it) - a(it,it) * betaOld + sigma0*lam1)/a(it,it);
      int sigma = (betaNew > 0 ? 1 : (betaNew < 0 ? -1 : 0 ) );
      if (sigma0!=sigma) betaNew=0;
      
      if (betaNew!=betaOld)
      {
        vec Mnew = M;
        for (int s=0;s<ns;s++) Mnew(s) += x0G(s,it) * (betaNew-betaOld);
        vec pnew = 1.0/(1.0+exp(-Mnew)); for (int s=0;s<ns;s++) {if (pnew(s) > 1.0-1e-5) pnew(s)=1.0-1e-5; else if (pnew(s)<1.0e-5) pnew(s)=1.0e-5;}
        beta(g+it*m0) = betaNew;
        LL=0;
        for (int s=0;s<ns;s++) LL += -y(s)*log(pnew(s))*w(s) - (1-y(s))*log(1.0-pnew(s))*w(s);
        LL += lam1*(abs(betaNew)-abs(betaOld));
        if (it>0) LL+= lam2*(betaNew-beta(g+(it-1)*m0))*(betaNew-beta(g+(it-1)*m0))-lam2*(betaOld-beta(g+(it-1)*m0))*(betaOld-beta(g+(it-1)*m0));
        if (it<(nt-1)) LL+= lam2*(betaNew-beta(g+(it+1)*m0))*(betaNew-beta(g+(it+1)*m0))-lam2*(betaOld-beta(g+(it+1)*m0))*(betaOld-beta(g+(it+1)*m0));
        if (LL<LLmin)
        {
          LLmin = LL;
          p = pnew;
          M = Mnew;
        }
        else
        {
          beta(g+it*m0) = betaOld;
        }
      }
    }
  }
  return(beta);
}


// [[Rcpp::export]]
vec _TimePointsPenalized_FitRound(mat x0, vec y, vec tV, double lam1, double lam2, vec beta, double Intercept, vec w, vec IndFor0,vec IndTFor0)
{
  IndFor0 = IndFor0-1;
  IndTFor0 = IndTFor0-1;
  int nt = tV.size();
  w=w/accu(w);
  int m = beta.size();
  int ns = y.size();
  int m0 = m/nt;
  vec M(ns); M.fill(Intercept);
  for (int s=0;s<ns;s++) for (int g=0;g<m0;g++) {M(s) = M(s) + x0(IndFor0(s),g) * beta(g+IndTFor0(s)*m0);}
  vec p(ns); p = 1.0/(1.0+exp(-M)); for (int s=0;s<ns;s++) {if (p(s) > 1.0-1e-5) p(s)=1.0-1e-5; else if (p(s)<1.0e-5) p(s)=1.0e-5;}

  double LL=0;
  for (int s=0;s<ns;s++) LL += -y(s)*log(p(s))*w(s) - (1-y(s))*log(1.0-p(s))*w(s);
  for (int g=1;g<m0;g++) for (int it=0;it<nt;it++) LL += lam1*abs(beta(g+it*m0));
  for (int g=1;g<m0;g++) for (int it=0;it<nt-1;it++) LL+= lam2*(beta(g+(it+1)*m0)-beta(g+it*m0))*(beta(g+(it+1)*m0)-beta(g+it*m0));
  double LLprev = 0;
  vec betaPrev = -beta;
  while (abs(LL-LLprev)/sqrt(LLprev*LLprev+LL*LL)>1.0e-5 | any(sgn(beta)!=sgn(betaPrev)))
  {
    LLprev = LL;
    betaPrev = beta;
    beta = SingleGeneRound(x0, y, tV, lam1, lam2, beta, Intercept, w, IndFor0,IndTFor0, M, LL);
    beta = GroupRound(x0, y, tV, lam1, lam2, beta, Intercept, w, IndFor0,IndTFor0, M, LL);
  }
  return(beta);
}




