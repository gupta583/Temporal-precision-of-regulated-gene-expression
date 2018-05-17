function U = matrixU(const,alpha,beta,theta,L,mu)


    Nr = const(1);
    nstar = const(2);
    tstar = const(3);
    H = const(4);
    J = const(5);



nrow = (Nr+1)*nstar;
ncol = nrow;

U = zeros(nrow,1);


r0 = exp(theta(1));
%r0 = 5/(exp(-theta(1))+1);
k  = exp(theta(2));
%k = 5/(exp(-theta(2))+1);

x = nstar-1;
for r=0:Nr
       i = (r+1) + x*(Nr+1);
       f_wr = alpha*r0^H/(r^H+r0^H);
       U(i,1) = f_wr;
    
end

end