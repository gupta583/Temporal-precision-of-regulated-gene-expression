clear all;

nstar = 15;
tstar = 20;
td = 0.7328;        % Average division time: from experimental data
sigmad = 0.1373;     % variance in division time from experimental data
a0s = logspace(log10(5),log10(5),1);    %Half maxima for activator 
ks_a = logspace(log10(0.2),log10(1),15); %production rate activator
ks_a1 = logspace(log10(0.2),log10(1),20); %production rate activator
r0s = logspace(log10(10),log10(10),1);  % Half maxima, K repressor
mus = logspace(log10(0.01),log10(0.5),15); % Decay rate 
mus_v = logspace(log10(0.01),log10(0.5),20); % Decay rate 
%figure properties
lw = 2;  % linewidth
fs = 15;  %fontsize

figure(11);
clf;

subplot(2,2,1);
%activator burst and decay. K = 10; H = 1;
a0 = a0s(1);
var_0 = dlmread('act_K10_mu0_h3'); % \mu = 0 genrate data from folder fig4a_data/ by setting mu = 0
var_1 = dlmread('act_K10_mu1_h3'); % genrate data from folder fig4a_data/ by setting \mu = k/2K 
var_2 = dlmread('act_K10_mu2_h3'); % genrate data from folder fig4a_data/ by setting \mu = k/K
var_3 = dlmread('act_K10_mu0_b3'); % genrate data from folder fig4a_data/ by setting \mu = 0 burst b = 3
var_5 = dlmread('act_K10_mu0_b5') % genrate data from folder fig4a_data/ by setting \mu = 0 burst b = 5

hold on; box on;
plot(ks_a*tstar, var_0*nstar/tstar^2,'-g','linewidth',lw)
plot(ks_a*tstar, var_1*nstar/tstar^2,'-b','linewidth',lw)
plot(ks_a*tstar, var_2*nstar/tstar^2,'-r','linewidth',lw)

plot(ks_a1*tstar, var_3*nstar/tstar^2,'--c','linewidth',lw)
plot(ks_a1*tstar, var_5*nstar/tstar^2,'--m','linewidth',lw)
plot([3 20], [1 1],'--k','linewidth',lw)
xlim([3 20])
ylim([0 5])
legend({'$\mu = 0$','$\mu = k/2K$','$\mu = k/K$',...
    '$\mu = 0$, $b =3$','$\mu = 0$, $b =5$'},'Interpreter','latex');

xlabel('Production rate, $kt_*$','Interpreter','latex')
ylabel('Variance, $\sigma_t^2x_*/t_* ^2$','Interpreter','latex')
set(gca,'fontsize',fs)
title('Activator')

subplot(2,2,2);
r0 = r0s(1);
box on; hold on;
title('Repressor')
var_0 = dlmread('rep_K10_mu0_h3');%genrate data from folder fig4b_data/ by setting  \mu = 0 
var_1 = dlmread('rep_K10_mu1_h3'); %k = muK/2
var_2 = dlmread('rep_K10_mu2_h3');%k = \muK
var_3 = dlmread('rep_mu0_p'); % k = 0; poisson genrate data from fig4b_data/withnoise

plot(mus*tstar, var_0*nstar/tstar^2,'-g','linewidth',lw)
plot(mus*tstar, var_1*nstar/tstar^2,'-b','linewidth',lw)
plot(mus*tstar, var_2*nstar/tstar^2,'-r','linewidth',lw)
plot(mus*tstar, var_3*nstar/tstar^2,'--g','linewidth',lw)
plot([0 10], [1 1],'--k','linewidth',lw)
ylim([0.5 2.2])
legend({'$k = 0$','$k = \mu K/2$','$k = \mu K$'...
    ,'$k = 0$, Poisson'},'Interpreter','latex');
xlabel('Decay rate, $\mu t_*$','Interpreter','latex')
ylabel('Variance, $\sigma_t^2x_*/t_* ^2$','Interpreter','latex')
set(gca,'fontsize',fs)




subplot(2,2,3)
% genrate trajectory from fig4c_data/activator and /repressor

traj = dlmread('act_traj'); 
t_a = traj(:,1);
x_a = traj(:,2);
a   = traj(:,3);

traj_ = dlmread('rep_traj');
t_r = traj_(:,1);
x_r = traj_(:,2);
r   = traj_(:,3);
plot([0 0.1],[0 0],'-k','linewidth',lw)
hold on
plot([0 0],[0 0.1],'--k','linewidth',lw)
plot(t_a/tstar,x_a,'-b','linewidth',lw);
hold on
plot(t_r/tstar,r/10,'--r','linewidth',lw);
plot(t_r/tstar,x_r,'--b','linewidth',lw);
hold on
plot(t_a/tstar,a/10,'-g','linewidth',lw);

plot([td td],[0 40],':k','linewidth',lw)
legend({'Activator','Repressor'},'interpreter','latex')
xlim([0 1.25])
ylim([0 35])
xlabel('Time, $t/t_*$','Interpreter','latex')
ylabel('Molecule number','Interpreter','latex')
set(gca,'fontsize',fs)
text(0.45,20,'$\bar{a}/10$','Interpreter','latex','fontsize',fs,'color','g')
text(1,11,'$\bar{x}$','Interpreter','latex','fontsize',fs,'color','b')
text(0.75,30,'$\bar{t}_d/t_*$','Interpreter','latex','fontsize',fs,'color','black')
text(0.1,15,'$\bar{r}/10$','Interpreter','latex','fontsize',fs,'color','r')


subplot(2,2,4);

%generate data from fig4d_data/activator and repressor
box on;
hold on;

tds = 0.1:0.1:0.9;
rep_var = dlmread('rep_var_td_h3');  %s simple division
act_var = dlmread('act_var_td_h3_1');
rep_var_nod = 2.36;
act_var_nod = 3.89;

plot(tds,act_var*nstar/tstar^2,'-g','linewidth',lw)
hold on;
plot(tds,rep_var*nstar/tstar^2,'-r','linewidth',lw)
plot(-[0.1 0.9],[act_var_nod act_var_nod]*nstar/tstar^2,'--k','linewidth',lw)
plot([0.1 0.9],[act_var_nod act_var_nod]*nstar/tstar^2,'--g','linewidth',lw)
plot([0.1 0.9],[rep_var_nod rep_var_nod]*nstar/tstar^2,'--r','linewidth',lw)
patch([td-sigmad td+sigmad td+sigmad td-sigmad],...
    [0 0 0.6 0.6],[0.75 0.75 0.75],'edgecolor','none')
alpha(0.3)
xlabel('Time, $t_d/t_*$','Interpreter','latex')
ylabel('Variance, $\sigma_t^2x_*/t_* ^2$','Interpreter','latex')
set(gca,'fontsize',fs)
xlim([0.1 0.9])
text(0.72,0.2,'Exp. $t_d/t_*$','fontsize',fs,'rotation',90,'Interpreter','latex')
legend({'Activator','Repressor','No division'},'interpreter','latex','location','NW')




