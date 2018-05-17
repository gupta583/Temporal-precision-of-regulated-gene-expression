clear all;

nstar = 15;
tstar = 20;

a0s = logspace(log10(290),log10(292.5),2);     % 1 to 15
ks = logspace(log10(0.2),log10(15.0),2);   %intitial 2 --> 150    0.2 to 1

for aa=1:length(a0s)
    a0 = a0s(aa);
    aa
    for kk=1:length(ks)
      k = ks(kk);
      
      
      Nc = fix(k*tstar + 10*sqrt(k*tstar));
      ncolumn = (Nc+1)*nstar;
      nrow = ncolumn;
      
      P0 = zeros(nrow,1);
      P0(1,1) = 1;
      
      
       l = @(alpha) (1*transpose(matrixU(a0,k,alpha,nstar,tstar))*...
                  (inv(matrixA(a0,k,alpha,nstar,tstar)))^2*P0 - tstar)^2;
       a1 = fminbnd(l,0,500000);
       
       mean = 1*transpose(matrixU(a0,k,a1,nstar,tstar))*...
              (inv(matrixA(a0,k,a1,nstar,tstar)))^2*P0;
      mean2= -1*2*transpose(matrixU(a0,k,a1,nstar,tstar))*...
              (inv(matrixA(a0,k,a1,nstar,tstar)))^3*P0;
          
      variance(aa,kk) = mean2 - mean^2;
      bestalpha(aa,kk) = a1;
      bestmean(aa,kk)  = mean;
        
    end
end



imagesc(log10(a0s),log10(ks),variance')
xlabel('log_{10}a0')
ylabel('log_{10}k')
colorbar;
%print(gcf,'-depsc', 'mat10var.eps')
%dlmwrite('actvarianceh1.dat',variance)




