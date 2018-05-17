function A = matrixA(a0,k,alpha,nstar,tstar)

      Nc = fix(k*tstar + 10*sqrt(k*tstar));
      ncolumn = (Nc+1)*nstar;
      nrow = ncolumn;
      A = zeros(nrow,ncolumn);

      
       for i=1:nrow
           
                  for j=1:ncolumn
     
                    if i == j
         
                      a = (j-1) - fix((j-1)/(Nc+1))*(Nc+1);
                      f = alpha*a^3/(a^3+a0^3); 
                      A(i,i) = -f-k;
     
                        if i < ncolumn
                          A(i+1,i) = k;
                        end
     
                        if i<=  ncolumn- (Nc+1)
                           A(i+Nc+1,j) = f;
                        end
           
                        if i== ncolumn
                           A(i,i) = -f;
                        end
           
                    end
        
                  end
       end
           
end
