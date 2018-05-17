clear all;

figure(1);
clf;

nstar = 15;
tstar = 20;

t=0:0.1:tstar;


Nr = 15;
alpha1 = 1.6922;
k1 = 0.1376;
r0 = 2.6371;

% analytic mean trajectory of x and r.
meanx1 = alpha1*log((r0*exp(k1*t)+Nr)/(r0+Nr))/k1;
meanr = Nr*exp(-k1*t);
t0r = log(Nr/r0)/k1;
N0r = Nr*exp(-k1*t0r)
%_________---______________________________________

subplot(2,2,4)
plot(t/tstar,meanx1,'b',t/tstar,meanr,'r','linewidth',2)
ylim([0 20])
xlabel('Time, $t/t_*$','Interpreter','latex')
ylabel('Molecule number','Interpreter','latex')
hold on;
set(gca,'fontsize',17)
plot([0 t0r]/tstar,[r0 r0],'--r',...
    [t0r t0r]/tstar, [N0r 0],'--r',...
    [0 5]/tstar,[Nr Nr],':k','linewidth',2)

text((tstar-2)/tstar, 3,'$\bar{r}$','Interpreter','latex','color','r'...
    ,'fontsize',15)
text((tstar-3)/tstar, Nr-2,'$\bar{x}$','Interpreter','latex','color','b'...
    ,'fontsize',15)
text(1/tstar, r0+1,'$K$','Interpreter','latex'...
    ,'fontsize',15)
text(0.65*tstar/tstar, r0-1.6,'$t_0/t_*$','Interpreter','latex'...
    ,'fontsize',15)
text(3/tstar,Nr+1,'$N$','Interpreter','latex','fontsize',15)
title('Repressor')


% mean trajectory for product and regulator(activator)
alpha2 = 2.1052;
k2 = 1.0;
a0 = 15;
t0a = a0/k2;
meanx2 = alpha2*(k2*t-a0*log((k2*t+a0)/a0))/k2;
meana  = k2*t;

%-------------------------------------


subplot(2,2,3)
plot(t/tstar,meanx2,'b',t/tstar,meana,'g','linewidth',2)
xlabel('Time, $t/t_*$','Interpreter','latex')
title('Activator')
hold on;
plot([0 a0/k2]/tstar,[a0 a0],'--g','linewidth',2)
plot([t0a t0a]/tstar, [0 k2*t0a],'g--','linewidth',2)
text((tstar-3)/tstar, Nr,'$\bar{a}$','Interpreter','latex','color','g'...
    ,'fontsize',15)
text((tstar-3)/tstar, Nr-4,'$\bar{x}$','Interpreter','latex','color','b'...
    ,'fontsize',15)
text(1/tstar, a0+1,'$K$','Interpreter','latex'...
    ,'fontsize',15)
text(0.59*tstar/tstar, a0-10,'$t_0/t_*$','Interpreter','latex'...
    ,'fontsize',15)
ylabel('Molecule number','Interpreter','latex')
set(gca,'fontsize',17)

%------------------
subplot(2,2,1)
image( imread('1ac.png') );
set(gca,'XtickLabel',[],'YtickLabel',[]);
set(gca,'xtick',[])
set(gca,'ytick',[])

%-------------------
clear all;
subplot(2,2,2)

Nr = 15;
nstar = 15;
tstar = 20;
k = 0.1391;

M = dlmread('1a.dat') %  Variance obtainted from fig_1_b_variance by creating ensemble from fig_1_b_ensemble.sh
M1 = dlmread('1.dat'); % sample trajectory obtained from fig_1_b_trajectory
M2 = dlmread('2.dat'); % sample trajectory obtained from fig_1_b_trajectory
M3 = dlmread('3.dat'); % sample trajectory obtained from fig_1_b_trajectory

j =1
for i=1:length(M1)
    
   if M1(i,2) < 16
       NEWM1(j,1) = M1(i,1)
       NEWM1(j,2) = M1(i,2)
       j = j+1
   end
end

j =1
for i=1:length(M2)
    
   if M2(i,2) < 16
       NEWM2(j,1) = M2(i,1)
       NEWM2(j,2) = M2(i,2)
       j = j+1
   end
end
j =1
for i=1:length(M3)
    
   if M3(i,2) < 16
       NEWM3(j,1) = M3(i,1)
       NEWM3(j,2) = M3(i,2)
       j = j+1
   end
end

stairs(NEWM1(:,1)/tstar,NEWM1(:,2),'linewidth',2)
hold on;
stairs(NEWM2(:,1)/tstar,NEWM2(:,2),'linewidth',2)
stairs(NEWM3(:,1)/tstar,NEWM3(:,2),'linewidth',2)
plot((0.5+M(:,1))/tstar,15+50*M(:,2),'linewidth',2)
text(41/tstar, nstar,'$x_*$','Interpreter','latex'...
    ,'fontsize',17)
xlim([0 40]/tstar)
xlabel('Time, $t/t_*$','Interpreter','latex')

r0 = 2.6371;
t0 = log(Nr/r0)/k;

plot([0 40]/tstar,[nstar nstar],':k','linewidth',2)
plot([tstar tstar]/tstar,[0 30],':k','linewidth',2)


text((tstar+0.1)/tstar,nstar+2, '$\rightarrow$','Interpreter','latex',...
    'fontsize',17)
text((tstar+0.5)/tstar,nstar+8, '$\sigma_t^2$','Interpreter','latex',...
    'fontsize',17)
ylabel('Molecule number, $x$','Interpreter','latex')
set(gca,'fontsize',17)






