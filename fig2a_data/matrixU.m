function U = matrixU(a0,k,alpha,nstar,tstar)


 Nc = fix(k*tstar + 10*sqrt(k*tstar));
 ncolumn = (Nc+1)*nstar;
 nrow = ncolumn;
 
 U = zeros(nrow,1);

for i=1:nrow
    
    if i > nrow -(Nc+1)
       a = (i-1) - fix((i-1)/(Nc+1))*(Nc+1);
       f = alpha*a^3/(a^3+a0^3);
       U(i,1) = f;
     end
    
end

end
