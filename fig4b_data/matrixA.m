function A = matrixA(const,alpha,beta,theta,L,mu)


    Nr = const(1);
    nstar = const(2);
    tstar = const(3);
    H = const(4);
    J = const(5);
    
    
    nrow = (Nr+1)*nstar;
    ncol = nrow;
    
A = zeros(nrow,ncol);

r0 = exp(theta(1));
%r0 = 5/(exp(-theta(1))+1);
k  = exp(theta(2));
%k = 5/(exp(-theta(2))+1);

for x_=0:nstar-1;
    for r_=0:Nr;
        
            for x=0:nstar-1;
               for r=0:Nr;
                 
                     
                   i = (r+1) + x*(Nr+1);
                   j = (r_+1) + x_*(Nr+1);
                    
                  
                   
                   if r_ == r  && x_ == x-1
                       f_rx = alpha*r0^H/(r^H+r0^H);
                       A(i,j) = f_rx;
                       
                   elseif r_ == r+1 && x_ ==x
                       A(i,j) = k*r_;
                       
                   %elseif r_ == r  && x_ == x+1
                   %    A(i,j) = mu*x_;
                   elseif r_ == r-1 && x_ == x
                       A(i,j) = mu;
                       
                   elseif r_ == r  && x_ == x
                      
                       f_xr =alpha*r0^H/(r^H+r0^H);
                       A(i,j) = -(k*r_+mu+f_xr) ; % all term
                   
                   end
        
           
                end
             end
            
    end
end
           
    
    