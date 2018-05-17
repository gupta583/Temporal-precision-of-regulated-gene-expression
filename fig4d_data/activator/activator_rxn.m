clear all;

tstar = 20;
nstar = 15;
tmax = 50;
H = 1;
Na = 0;
traj_time= 50;
dt = 0.25;
t0 = traj_time;
sigma = 0.1373;
td = 0.7328;

const = [nstar tstar tmax H Na sigma];

ensemble = 10000;
a0 = 292.5;
k = 15;
alpha = 1.8324;


figure(1);
clf;
H =3;

for hh=1:9
        
        td = hh*0.1;
        clear mod x_ a_ t_
        const = [nstar tstar tmax H Na sigma td];
        variable = [a0,k];

        l = @(alpha) (meanfpt_act(stochatic_act(const,variable,alpha,ensemble),nstar,ensemble) ...
                       -tstar)^2;
       
        [a1_ sl] = fminbnd(l,0,10);
        mod = stochatic_act(const,variable,a1_,ensemble);
        
        [fpt V_fpt]= meanfpt_act(mod,nstar,ensemble);
        
        t0 = (td-sigma)*tstar;
        
        x_(:,:) = mod(:,:,1);
        a_(:,:) = mod(:,:,2);
        t_(:,:) = mod(:,:,3);
        
        
        
        [x,a,t] =meantraj_act(x_,a_,t_,t0,dt,ensemble);
        
        
        
        if length(x) ==2
        s1 = (x(2) - x(1))/dt;
        end
        
        clear x a t;
        
        
        t0 = (td+sigma)*tstar;
        [x,a,t] =meantraj_act(x_,a_,t_,t0,dt,ensemble);
        if length(x) ==2
        s2 = (x(2) - x(1))/dt;
        end
        
        s(hh) = s2-s1
         
        variance(1,hh) = V_fpt;
        bestalpha(1,hh) = a1_; 
        bestmean(1,hh)  = fpt;
        
        clear x a t;
        
        [x,a,t] =meantraj_act(x_,a_,t_,0,25,ensemble);
        
        
        plot(t,x,'linewidth',2)
        hold on
        plot(t,a/10,'linewidth',2)
        
 %       a1 = k*(td-sigma)*tstar;
 %       a2 = k*(td+sigma)*tstar;
 %       s2 = a1_*a2^H/(a0^H+a2^H);
 %       s1 = a1_*a1^H/(a0^H+a1^H);
        
 %       s_ana(hh)= s2-s1
        
        
        
end




plot(tstar*[td td], [0 15],'--k','linewidth',2)
plot(0.5*tstar*[td td], [0 15],'--k','linewidth',2)
plot(tstar*[td-sigma td-sigma], [0 15],'--k','linewidth',2)

plot(tstar*[td+sigma td+sigma], [0 15],'--k','linewidth',2)


   

