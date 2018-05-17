function U = matrixU(Nr,r0,alpha,nstar)

ncolumn = (Nr+1)*nstar;
nrow = ncolumn;
U = zeros(nrow,1);

for i=1:nrow
    
    if i > nrow -(Nr+1)
        r = (i-1) - fix((i-1)/(Nr+1))*(Nr+1);
        f = alpha*r0^3/(r^3+r0^3);
        U(i,1) = f;
    end
end

end
