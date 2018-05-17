clear all;

nstar = 15;
tstar = 20;
H = 3;


a0s = logspace(log10(10),log10(10),1);
ks = logspace(log10(0.2),log10(1),20);    


const = [nstar tstar H];
    
    bs = 5:1:5;
   
for bb = 1:length(bs)
        b = bs(bb)
for aa=1:length(a0s)
    a0 = a0s(aa);
    for kk=1:length(ks)
      k = ks(kk)
      mu = 0;%2*k/a0;
      
      
      Nc = fix(b*k*tstar + 20*sqrt(b*k*tstar));
      ncolumn = (Nc+1)*nstar;
      nrow = ncolumn;
      
      P0 = zeros(nrow,1);
      P0(1,1) = 1;
      
      theta = [log(a0) log(k)];
       l = @(alpha) (1*transpose(matrixU(const,theta,alpha,mu,b))*...
                  (inv(matrixA(const,theta,alpha,mu,b)))^2*P0 - tstar)^2;
       a1 = fminbnd(l,0,50);
       
       mean = 1*transpose(matrixU(const,theta,a1,mu,b))*...
              (inv(matrixA(const,theta,a1,mu,b)))^2*P0;
      mean2= -1*2*transpose(matrixU(const,theta,a1,mu,b))*...
              (inv(matrixA(const,theta,a1,mu,b)))^3*P0;
          
      variance(aa,kk) = mean2 - mean^2;
      bestalpha(aa,kk) = a1;
      bestmean(aa,kk)  = mean;
        
    end
end
var_burst(bb) = variance(1,1)
    end
    


%figure(2);
plot(bs,var_burst,'linewidth',2)
%plot(ks*tstar, variance(1,:)*nstar/tstar^2,'-b','linewidth',2)
xlabel('$kt_*$','Interpreter','latex')
%xlabel('Burst size $b$','Interpreter','latex')
ylabel('Variance, $\sigma_t^2x_*/t_* ^2$','Interpreter','latex')
set(gca,'fontsize',20)
hold on
%hold on;


%xlabel('Production rate, $k t_*$','Interpreter','latex')
%ylabel('Variance, $\sigma_t^2x_*/t_* ^2$','Interpreter','latex')
%legend({'$k = 0$','$k =\mu K/2$','$k=\mu K $'},'Interpreter','latex')
%set(gca,'fontsize',20)
%legend({'$\mu = 0$','$\mu =k/K$','$\mu=2k/K $'},'Interpreter','latex')

