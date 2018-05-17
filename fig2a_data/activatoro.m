clear all;

nstar = 15;
tstar = 20;

a0s = logspace(log10(1),log10(15),150);     % K, Half maxima
ks = logspace(log10(0.2),log10(10),150);    % Production rate, k

for aa=1:length(a0s)
    a0 = a0s(aa);
    aa
    for kk=1:length(ks)
      k = ks(kk);
      
      
      Nc = fix(k*tstar + 10*sqrt(k*tstar));     %cutoff for activator mean + 10*standard deviation
      ncolumn = (Nc+1)*nstar;
      nrow = ncolumn;
      
      P0 = zeros(nrow,1);
      P0(1,1) = 1;
      
      % fix the target condition by finding correct alpha
       l = @(alpha) (1*transpose(matrixU(a0,k,alpha,nstar,tstar))*...
                  (inv(matrixA(a0,k,alpha,nstar,tstar)))^2*P0 - tstar)^2;    
       a1 = fminbnd(l,0,500000);          
       
       mean = 1*transpose(matrixU(a0,k,a1,nstar,tstar))*...          % mean first passage time
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
dlmwrite('fig2a.dat',variance)




