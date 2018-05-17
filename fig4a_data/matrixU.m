function U = matrixU(const,theta,alpha,mu,b)

  a0 = exp(theta(1));
  k = exp(theta(2));
  nstar = const(1);
  tstar = const(2);
  H = const(3);
  
  
  Nc = fix(b*k*tstar + 20*sqrt(b*k*tstar));
  ncolumn = (Nc+1)*nstar;
  nrow = ncolumn;
 
  U = zeros(nrow,1);

  x = nstar-1;

  for a=0:Nc
       i = (a+1) + x*(Nc+1);
       f_wr = alpha*a^H/(a^H+a0^H);
       U(i,1) = f_wr;
    
  end



end
