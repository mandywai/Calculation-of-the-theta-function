# Calculation-of-the-theta-function
Following the result by Reineke in Poisson automorphisms and quiver moduli (arXiv: 0804.3214), there is the MAPLE program to compute theta fuction as in definition 4.1 in the paper.

pmumut:= proc(q,d1,d2,m,inf)
  #This is the p_mu(t) in Prop4.1 (3)
  
  local i,pmu,pd;
    pmu:=1;
    
    for i from 1 to inf do
    
      pd:=simplify(pdq(m,i*d1,i*d2,q));
      
      pmu:=pmu + pd*T^i;
      
    od;
    
  pmu;
  
end;

mutwistedpt:= proc(q,d1,d2,m,e1,e2,inf)
  # This is the twisted pmut by the twisted in P.4 of the paper,{e,d}=m*(-e2d1+d2e1) this is for m arrow
  local i,pmu,pd;
    pmu:=1;
    for i from 1 to inf do
      pd:=simplify(pdq(m,i*d1,i*d2,q));
      pmu:=pmu + q^(m*i*(-e2*d1+d2*e1))*pd*T^i; 
    od;
  pmu;
end;


muqt:=proc(q,d1,d2,m,e1,e2,inf)
  local pt,c,inv,k,p,twp,twpq,mult,deg,n,ans,i,r,cc;
    pt:=pmumut(q,d1,d2,m,inf); inv:=1; p:=1; twpq:=1; ans:=1;
     for k from 1 to inf do
        c:=coeftayl(pt,T=0,k);
        p:=p+ c*q^((2-m)*k*(k-1)/2) *T^k;
     od;
    inv:=convert(taylor(1/p,T,inf+1),polynom);
    twp:=mutwistedpt(q,d1,d2,m,e1,e2,inf);
      for r from 1 to inf do
        twpq:=twpq +simplify(coeftayl(twp,T=0,r)*q^((2-m)*r*(r-1)/2))*T^r;
      od;
    mult:=convert(taylor(twpq*inv,T,inf),polynom);
    deg:=degree(mult,T);
      for n from 1 to deg do 
        cc:=coeftayl(mult,T=0,n);
        ans:=ans + simplify(cc*q^((m-2)*n*(n-1)/2))*t^n;
      od;
   print(subs(q=1,ans));
   ans;  
end;
