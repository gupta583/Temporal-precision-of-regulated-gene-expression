function A = matrixA(const,theta,alpha,mu,b)

    
      a0 = exp(theta(1));
      k = exp(theta(2));
      

      nstar = const(1);
      tstar = const(2);
      H = const(3);
      
    


      Nc = fix(b*k*tstar + 20*sqrt(b*k*tstar));
      ncolumn = (Nc+1)*nstar;
      nrow = ncolumn;
      A = zeros(nrow,ncolumn);

  for x_=0:nstar-1;
    for a_=0:Nc;
        
            for x=0:nstar-1;
               for a=0:Nc;
                 
                     
                   i = (a+1) + x*(Nc+1);
                   j = (a_+1) + x_*(Nc+1);
      
                   if a_ == a  && x_ == x-1
                       f_ax = alpha*a^H/(a^H+a0^H);
                       A(i,j) = f_ax;
                       
                   %elseif a_ == a-1 && x_ ==x
                   %    A(i,j) = k;
        
                   elseif a_ == a+1 && x_ == x
                       A(i,j) = mu*a_;
                       
                   elseif a_ == a  && x_ == x
                      
                       f_xa =alpha*a^H/(a^H+a0^H);
                       k_b = k;
                       A(i,j) = -(k_b + mu*a_ + f_xa) ; % all term
                   end
                   
                   for delta=1:a
                       if a_ == a-delta && x_ ==x
                           i = (a+1) + x*(Nc+1);
                           j = (a_+1) + x_*(Nc+1);
                           %geo_pro = (1/(b+1))*(b/(b+1))^delta;
                           geo_pro = (1/b)*((b-1)/b)^(delta-1);
                           A(i,j) = k*geo_pro;
                       end
                   end
                   
                   
               end
            end
    end
  end
  
       
      
           
end
