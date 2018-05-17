function A = matrixA(Nr,r0,k,alpha,nstar)

ncolumn = (Nr+1)*nstar;
nrow = ncolumn;
A = zeros(nrow,ncolumn);

for i=1:nrow         
    
 for j=1:ncolumn
     if i == j
         
        r = (j-1) - fix((j-1)/(Nr+1))*(Nr+1)
        rp = r + 1;
     
          if rp == Nr+1
            rp = 0;
          end
        f = alpha*r0/(r+r0); 
         A(i,i) = -f-r*k;
     
           if i < ncolumn
             A(i,i+1) = k*rp;
           end
     
           if i<=  ncolumn- (Nr+1)
             A(i+Nr+1,j) = f;
           end
     
     end
        
 end
           
end
