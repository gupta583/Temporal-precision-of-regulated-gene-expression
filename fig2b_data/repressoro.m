clear all;

nstar = 15;
tstar = 20;
Nr = 365;

r0s = logspace(log10(45.0),log10(50),20);
ks = logspace(log10(0.09),log10(0.5),20);

ncolumn = (Nr+1)*nstar;
nrow = ncolumn;

P0 = zeros(nrow,1);
P0(Nr+1,1) = 1;
matrixA(2,1,0,1,2)
pause


for rr=1:length(r0s)
    r0 = r0s(rr)
    for kk=1:length(ks)
        k = ks(kk)
        matrixA(Nr,r0,k,alpha,nstar)
        pause;
        l = @(alpha) (1*transpose(matrixU(Nr,r0,alpha,nstar))*...
                  (inv(matrixA(Nr,r0,k,alpha,nstar)))^2*P0 - tstar)^2;
              a1 = fminbnd(l,0,50000);
              
 mean = 1*transpose(matrixU(Nr,r0,a1,nstar))*...
          (inv(matrixA(Nr,r0,k,a1,nstar)))^2*P0;
 mean2= -1*2*transpose(matrixU(Nr,r0,a1,nstar))*...
          (inv(matrixA(Nr,r0,k,a1,nstar)))^3*P0;
 variance(rr,kk) = mean2 - mean^2;
 bestalpha(rr,kk) = a1;
 bestmean(rr,kk)  = mean;
        
    end
end

imagesc(log10(r0s),log10(ks),variance')
xlabel('log_{10}r_0')
ylabel('log_{10}k')
colorbar;
%print(gcf,'-depsc', 'mat9-var.eps')
%dlmwrite('varianceh1300point.dat',variance)

