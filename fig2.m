clear all;
figure(1);
clf;
repvariance = dlmread('repvarianceh1r01.dat');

nstar = 15;
tstar = 20;
Nr = 15;

r0s = logspace(log10(1),log10(15),300); % r0 row
repks = logspace(log10(0.01),log10(5),300); %k's 


[M,I] = min(repvariance(:));
[I_row,I_col] = ind2sub(size(repvariance),I);

xn = r0s(I_row);    %global minima
yn = repks(I_col);  %global minima

xa = Nr/exp(2);     %h --> inf analytical minima
ya = (nstar*exp(2)/(2*Nr*tstar)) + (2/tstar);

ah = subplot(1,2,2);
currentPos = get(ah,'Position');
set(ah,'Position',currentPos + [0.05,0,0,0])
get(ah,'Position')
%imagesc(log10(r0s),log10(repks*tstar),(repvariance*nstar/tstar^2)');
contourf(log10(r0s),log10(repks*tstar),repvariance'*nstar/tstar^2)
hold on;
%plot(log10(xn),log10(yn),'.','color','magenta','MarkerSize',25);
plot(log10(xa),log10(ya*tstar),'.','color','white','MarkerSize',30);
axis('square')
%text(log10(0.5),log10(7*tstar),'B','fontsize',14)
xlabel('Half-maximal number, $K$','Interpreter','latex');
ylabel('Degradation rate, $\mu t_*$','Interpreter','latex');
title('Repressor')
c = colorbar;
%c.Label.String = '\sigma_t^2';
ylabel(c,'Timing variance, $\sigma_t^2x_*/t_*^2$','fontsize',11,...
       'Interpreter','latex','fontsize',15)

set(gca,'XTickLabel',{ '10^0', '10^1' });
set(gca,'XTick',[log10(10^0) log10(10^1) ]);


set(gca,'YTickLabel',{'10^{0}', '10^{1}', '10^2'  });
set(gca,'YTick',[log10(10^0) log10(10^1) log10(10^2)   ],...
    'ydir','normal');

set(gca,'fontsize',15)
hold off;

subplot(1,2,1);

a0s = logspace(log10(1),log10(15),150);
actks = logspace(log10(0.2),log10(1.0),150);

actvariance = dlmread('actvarianceh1a01.dat');

axa = 1:0.1:15;
aya = (axa + nstar/2)/tstar;



%imagesc(log10(a0s),log10(actks*tstar),(actvariance*nstar/tstar^2)');
contourf(log10(a0s),log10(actks*tstar),actvariance'*nstar/tstar^2)
hold on;

plot(log10(axa),log10(aya*tstar),'--','color','white','linewidth',2);
%xlim(log10([.1 15]))
ylim(log10([4 20]))
axis('square')

%text(log10(0.5),log10(1.13*tstar),'A','fontsize',14)
xlabel('Half-maximal number, $K$','Interpreter','latex');
ylabel('Production rate, $kt_*$','Interpreter','latex');
title('Activator')
c = colorbar;
%c.Label.String = '\sigma_t^2';
ylabel(c,'Timing variance, $\sigma_t^2x_*/t_*^2$','fontsize',11,...
 'Interpreter','latex','fontsize',15)

%set(gca,'XTickLabel',{  '10^0', '10^1' });
%set(gca,'XTick',[ log10(10^0) log10(10^1) ]);


%set(gca,'YTickLabel',{'10^{-2}', '10^{-1}', '10^0'  });
%set(gca,'YTick',[log10(10^-2) log10(10^-1) log10(10^0)   ]);


set(gca,'ydir','normal',...
    'xtick',[0 1],'xticklabel',{'10^0','10^1'},...
    'ytick',log10(.2*tstar:.2*tstar:1*tstar),...
    'yticklabel',.2*tstar:.2*tstar:1*tstar)
set(gca,'fontsize',15)
print(gcf,'-dpdf', 'review1')




