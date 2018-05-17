clear all;

tstar = 20;
nstar = 15;
tmax = 50;
H = 1;
Nr = 365;
traj_time= 50;
dt = 0.25;
sigma = 0.1373;
td = 0.7328;

const = [nstar tstar tmax H Nr sigma];

ensemble = 5000;
r0 = Nr/exp(2);
k = 0.107587;
alpha = 1.8324;


%r0s = logspace(log10(5),log10(5),1);
%ks = logspace(log10(0.01),log10(0.5),20);

figure(1);
clf;




for hh=3:3
        
        H= hh;
        clear mod x_ r_ t_
        const = [nstar tstar tmax H Nr sigma td];
        variable = [r0,k];

        l = @(alpha) (meanfpt(stochatic(const,variable,alpha,ensemble),nstar,ensemble) ...
                       -tstar)^2;
        [a1 sl] = fminbnd(l,1,10);
        mod = stochatic(const,variable,a1,ensemble);
        
        [fpt V_fpt]= meanfpt(mod,nstar,ensemble);
        
        t0 = (td-sigma)*tstar;
        
        %mod = stochatic(const,variable,alpha,ensemble);
        
       x_(:,:) = mod(:,:,1);
       r_(:,:) = mod(:,:,2);
       t_(:,:) = mod(:,:,3);
        
        [x,r,t] =meantraj(x_,r_,t_,t0,dt,ensemble);
        
        
        
        if length(x) ==2
        s1 = (x(2) - x(1))/dt;
        end
        
        clear x r t;
        
        
        t0 = (td+sigma)*tstar;
        [x,r,t] =meantraj(x_,r_,t_,t0,dt,ensemble);
        if length(x) ==2
        s2 = (x(2) - x(1))/dt;
        end
        
        s(hh) = s2-s1
        
%        r1=Nr*exp(-k*(td-sigma)*tstar);
%        r2=Nr*exp(-k*(td+sigma)*tstar);
        
%        s_ana(hh) = a1*r0^H*(r1^H-r2^H)/(r0^H+r1^H)/(r0^H+r2^H);
         
        variance(1,hh) = V_fpt;
        bestalpha(1,hh) = a1;
        bestmean(1,hh)  = fpt;
        
        clear x r t;
        
        [x,r,t] =meantraj(x_,r_,t_,0,25,ensemble);
        
        
        plot(t,x,'linewidth',2)
        hold on
        plot(t,r/10,'linewidth',2)
        
end

%plot(tstar*[td td], [0 15],'--k','linewidth',2)
%plot(0.5*tstar*[td td], [0 15],'--k','linewidth',2)
%plot(tstar*[td-sigma td-sigma], [0 15],'--k','linewidth',2)

%plot(tstar*[td+sigma td+sigma], [0 15],'--k','linewidth',2)
%mod = stochatic(const,variable,alpha,ensemble);

%x_(:,:) = mod(:,:,1);
%r_(:,:) = mod(:,:,2);
%t_(:,:) = mod(:,:,3);


 %[fpt V_fpt]= meanfpt(stochatic(const,variable,alpha,ensemble),nstar,ensemble)   



%[fpt V_fpt] = meanfpt(x_,r_,t_,nstar,ensemble)
     



%[x,r,t] =meantraj(x_,r_,t_,traj_time,dt,ensemble);

%plot(t,x,'g',t,r,'r')

%p = polyfit(t,x,1)

%f = fit(t',r','exp1')

   

