clear all;

nstar = 15;
tstar = 20;
Nr = 15;
Nc = 45;
H = 3;
J = 1;

r0s = logspace(log10(10),log10(10),1);
ks = logspace(log10(0.01),log10(0.5),15);

ncolumn = (Nc+1)*nstar;
nrow = ncolumn;

P0 = zeros(nrow,1);
r =Nr; x = 0;


int_cond = 1 +r + x*(Nr+1);
P0(int_cond,1) = 1;


%P0 = zeros(nrow,1);
%P0(Nr+1,1) = 1;

beta= 0;
%mu = 0.5;
L = 10000;
r0 = Nr/exp(2);

%mu is production rate in code actualy in plot decay rate 



const = [Nc nstar tstar H J];

for rr=1:length(r0s)
    r0 = r0s(rr);
    for kk=1:length(ks)
        k = ks(kk);
        mu = 0;%0.5*r0*k;       
        theta = [log(r0) log(k)];
        l = @(alpha) (1*transpose(matrixU(const,alpha,beta,theta,L,mu))*...
                  (inv(matrixA(const,alpha,beta,theta,L,mu)))^2*P0 - tstar)^2;
              a1 = fminbnd(l,0,50000);
              
 mean = 1*transpose(matrixU(const,a1,beta,theta,L,mu))*...
          (inv(matrixA(const,a1,beta,theta,L,mu)))^2*P0;
 mean2= -1*2*transpose(matrixU(const,a1,beta,theta,L,mu))*...
          (inv(matrixA(const,a1,beta,theta,L,mu)))^3*P0;
 variance(rr,kk) = mean2 - mean^2;
 bestalpha(rr,kk) = a1;
 bestmean(rr,kk)  = mean;
        
    end
end



figure(2);

plot(log10(ks*tstar), variance(1,:)*nstar/tstar^2,'-or','linewidth',2)
xlabel('Decay rate, $\mut_*$','Interpreter','latex')
ylabel('$\sigma_t^2x_*/t_* ^2$','Interpreter','latex')

hold on;

%print(gcf,'-depsc', 'mat9-var.eps')
%dlmwrite('varianceh1300point.dat',variance)

xlabel('Decay rate, $\mu t_*$','Interpreter','latex')
ylabel('Variance, $\sigma_t^2x_*/t_* ^2$','Interpreter','latex')
%legend({'$k = 0$','$k =\mu K/2$','$k=\mu K $'},'Interpreter','latex')
set(gca,'fontsize',20)
