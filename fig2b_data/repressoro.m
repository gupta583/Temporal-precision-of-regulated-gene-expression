clear all;

nstar = 15;
tstar = 20;
Nr = 15;

r0s = logspace(log10(1),log10(15),300);    %Half maxima K
ks = logspace(log10(0.01),log10(5),300);    % Decay rate, \mu

ncolumn = (Nr+1)*nstar;
nrow = ncolumn;

P0 = zeros(nrow,1);
P0(Nr+1,1) = 1;


for rr=1:length(r0s)
    r0 = r0s(rr)
    for kk=1:length(ks)
        k = ks(kk)
      
        l = @(alpha) (1*transpose(matrixU(Nr,r0,alpha,nstar))*...
                  (inv(matrixA(Nr,r0,k,alpha,nstar)))^2*P0 - tstar)^2;
              a1 = fminbnd(l,0,50000);    % set the target condition at t = t_* | x = x_*
              
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
dlmwrite('fig2b.dat',variance)

