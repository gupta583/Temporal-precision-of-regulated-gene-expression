clear all;

tstar = 20;
nstar = 15;
tmax = 70;
H = 3;
Nr = 15;
traj_time= 50;
dt = 0.20;

const = [nstar tstar tmax H Nr];

ensemble = 10000;


r0s = logspace(log10(10),log10(10),1);
ks = logspace(log10(0.01),log10(0.5),15);



%for i=1:ensemble
%            [x,r,t] = stochatic(const,variable,alpha);
%            x_(1:length(x),i) = x;
%            r_(1:length(x),i) = r;
%            t_(1:length(x),i) = t;     
%         end



for rr=1:length(r0s)
    r0 = r0s(rr);
    for kk = 1:length(ks)
        k = ks(kk);
        variable = [r0,k];

        l = @(alpha) (meanfpt(stochatic(const,variable,alpha,ensemble),nstar,ensemble) ...
                       -tstar)^2;
        a1 = fminbnd(l,0,5)
        mod = stochatic(const,variable,a1,ensemble);
        [fpt V_fpt]= meanfpt(mod,nstar,ensemble);
        variance(rr,kk) = V_fpt
        bestalpha(rr,kk) = a1;
        bestmean(rr,kk)  = fpt;
     end
end

[sz1 sz2 sz3] = size(mod);
nr_noise(1:sz2)= mod(1,1:sz2,2);
%mod = stochatic(const,variable,alpha,ensemble);

%x_(:,:) = mod(:,:,1);
%r_(:,:) = mod(:,:,2);
%t_(:,:) = mod(:,:,3);



plot(ks*tstar, variance(1,:)*nstar/tstar^2,'-g','linewidth',2)
xlabel('Decay rate, $\mu t_*$','Interpreter','latex')
ylabel('Variance, $\sigma_t^2x_*/t_* ^2$','Interpreter','latex')
set(gca,'fontsize',20)



   

