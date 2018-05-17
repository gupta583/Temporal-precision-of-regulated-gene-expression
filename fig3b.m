clear all;

M = dlmread('exp_data.csv');
t = M(:,1);
x = M(:,2);
area = 0;

Nr = 15;
%nstar = 15;

delxmax = 12;
delxmin = 3;

delxs = delxmin:1:delxmax;

oo = 1;
clf;
hold on


for i=1:length(delxs)

    delx = delxs(i);
    nstar = 1;
    
    
    o = 1;
    ll =1;
    while nstar < 25
        nstar = (2*(o-1)+1)*delx*0.5;
        o = o + 1;
        if nstar > 10
         nstars(ll) = nstar;
         ll = ll+1;
        end
    end
    
    
     for j=1:length(x)
         newx(j) = (2*fix(x(j)/delx)+1)*delx*0.5;
     end
     
     
   
     for j=1:length(nstars)-1
         nstar = nstars(j);
         
         l =1;
         for k=1:length(newx)
              if nstar ==newx(k)
                  tstars(l) = t(k);
                  l = l + 1;
              end
         end
         
        
         tstar(oo) = mean(tstars(1:l-1));
         variance(oo) = var(tstars(1:l-1));
         
         variance(oo) = nstar*variance(oo)/tstar(oo)^2;
        
         var1 = variance(oo);
         ll = 1;
         area = 0;
         
         while x(ll) < nstar
             area = area + (x(ll)+x(ll+1))*(t(ll+1)-t(ll))*0.5;
             ll = ll + 1;
         end
         
         area = area + (nstar-x(ll-1))*(tstar(oo)-t(ll-1))*0.5;
         
         area = 2*area/(nstar*tstar(oo));
         
         rho(oo) = area;
         
         
          oo = oo + 1;
          
 %         plot(area,var1,'.','MarkerSize',2*nstar,...
 %             'color',[delx/13+0.06 rand rand])
         
         
%          plot(area,var1,'.','MarkerSize',3*nstar,...
%              'color',[0 0 (delx-3)/10])
          
          
       hold on
       
     [area var1]
         
     end
     
      clear nstars
      clear tstar
     
     

   
end

mean_rho = mean(rho)
std_rho= std(rho)
varnce = mean(variance)
std_variance= std(variance)

h3 = plot(mean_rho,varnce,'ob','markersize',14,'linewidth',3)

abar_nstar = [0.5 2 20];
rbar_nstar = [0.5 2 20];
rhostar = 0:0.01:1.0;



green_shade = [0.4 0.7 1.0];

for i=1:length(abar_nstar)
     cons = (exp(2)-1)/rbar_nstar(i)/8;
     sminstar_act = rhostar.^2 + (1-rhostar)*0.5/abar_nstar(i);
     sminstar_rep = rhostar.^2 + cons*(1-rhostar).^3;
     h1(i)=plot(rhostar,sminstar_act,'color',[0 green_shade(i) 0],'linewidth',2);
     h2(i)=plot(rhostar,sminstar_rep,'--','color',[green_shade(i) 0 0],'linewidth',2);
     
     
     
     lstr1{i} = ['$\langle a \rangle/{x_*} =' num2str(abar_nstar(i)) '$'];
     lstr2{i} = ['$\langle r \rangle/{x_*} =' num2str(rbar_nstar(i)) '$'];
    % lstr3{i} = ['$n_* =' num2str(5*(i+1)) '$' ];
end


%     h3(1)= plot(-1,-1,'.k','Markersize',30);
%     lstr3{1} = ['$x_* = 10$' ];
     
%      h3(2)= plot(-1,-1,'.k','Markersize',60);
%         lstr3{2} = ['$x_* = 20$' ];
         
lstr3 = ['Experimental']         


legend([h3 h1 h2],[lstr3 lstr1 lstr2],'Interpreter','latex',...
    'Position',[0.71 0.305 0.01 0.2])

axis('square')
%h2,lstr2)
%c = colorbar;
%ylabel(c,'Bin size, $\Delta x$','Interpreter','latex')
%caxis([delxmin delxmax])

%cvec = 0:0.01:100;
%map = zeros(length(cvec),3);
%for i =1:length(cvec)
%map(i,3) = i/length(cvec);
%end
%colormap(map)

plot(1,1,'s','color','black','markersize',20,'linewidth',2)
text(0.44,1,'No regulation $\rightarrow$','Interpreter','latex','fontsize',23)
plot([mean_rho-std_rho mean_rho+std_rho],[varnce varnce],'-b','linewidth',3)
plot([mean_rho mean_rho],[varnce-std_variance varnce+std_variance],'-b','linewidth',3)

ylim([0 1.2])
box on;
%cons = nstar*(exp(2)-1)*0.25/Nr


%sminstar_rep = 4*rhostar.^2 + 1.5973*(1-2*rhostar).^2;
%plot(rhostar,sminstar_act1...
%     ,rhostar,sminstar_act2...
%     ,rhostar,sminstar_act3,'color',[0 1 0]...
%    ,rhostar,sminstar_rep);

xlabel('Linearity, $\rho$','Interpreter','latex')
ylabel('Timing variance, $ \sigma_t^2x_*/t_*^2$'...
    ,'Interpreter','latex')
set(gca,'Fontsize',23)
%a = get(gca,'XTickLabel');
%set(gca,'fontsize',24)


%print(gcf,'-dpdf', 'fig3b')


 
