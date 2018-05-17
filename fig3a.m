clear all;

nstar1 = 15;
tstar1 = 20;
M = xlsread('exp_mig.xlsx');
shift = 3.2;
%figure(1);
%clf;

%plot(M(1:98,1),M(1:98,2),'or',M(99:123,1),M(99:123,2)...
%    ,'og',M(124:126,1),M(124:126,2),'oy',...
%    M(127:175,1),M(127:175,2),'ob',...
%    M(176:179,1),M(176:179,2),'ok',...
%    M(180:205,1),M(180:205,2),'om');
%legend({'QR','QR.p','QR.p about to divide',...
%    'QR.p data Ni Ji','QR.p dividing','QR.pa'},'Interpreter','latex')

mrnab = [M(99:123,2)'];
timeb = [M(99:123,1)'];

mrnaa = [M(180:205,2)'];
timea = [M(180:205,1)'];

mrnad = [M(176:179,2)' M(124:126,2)'];
timed = [M(176:179,1)' M(124:126,1)'];

timeb = shift-timeb;
timea = shift-timea;
timed = shift-timed;

x(:,2) = [mrnab mrnad mrnaa];
x(:,1) = [timeb timed timea];


c = 1;
for i=1:length(x(:,1))
    
    if x(i,2) <= 25 && x(i,2) >= 10
        b(c,1) = x(i,1);
        b(c,2) = x(i,2);
        c = c+1;
    end
    
end
tstar = mean(b(:,1));


M(:,1) = (shift - M(:,1))*tstar1/tstar;
x(:,1) = [timeb timed timea]*tstar1/tstar;
timea  = timea*tstar1/tstar;
timeb  = timeb*tstar1/tstar;
timed  = timed*tstar1/tstar;

% check tstar
c = 1;
for i=1:length(x(:,1))
    
    if x(i,2) <= 25 && x(i,2) >= 10
        b(c,1) = x(i,1);
        b(c,2) = x(i,2);
        c = c+1;
    end
    
end
tstar = mean(b(:,1));
sigmat = std(b(:,1));

td = mean(timed);
sigmad = std(timed);



figure(3);
clf;
hold on
box on;
plot(timeb,mrnab,'oc','markerfacecolor','c','markersize',8)
plot(timed,mrnad,'ok','markerfacecolor','k','markersize',8)
plot(timea,mrnaa,'ob','markerfacecolor','b','markersize',8)

legend({'Before cell division','During cell division','After cell division'},...
    'Interpreter','latex','location','northwest','fontsize',18)
plot([td td],[0 35],'--k','linewidth',2)
plot([td-sigmad td-sigmad],[0 35],':k','linewidth',2)
text(td-sigmad,10,'$\longleftarrow$','Interpreter','latex'...
    ,'fontsize',25)
text(td-sigmad/2-0.2,11.5,'$\sigma _d$','Interpreter','latex',...
    'fontsize',20)


text(td-sigmad/2,15,'Cell division','Interpreter','latex','Rotation',90,...
       'fontsize',20)

 % migration stop
 
 %plot([td+sigmad-0.35 td+sigmad-0.35],[0 35],'--m',... by eye
 %   [tstar+sigmat-0.6 tstar+sigmat-0.6], [0 35],'--m')
temp = (tstar+sigmat-0.6 +td+sigmad-0.35)/2;
temp1 = tstar+sigmat-0.6;
temp2 = td+sigmad-0.35;

patch([temp1 temp1 temp2 temp2],...
    [0 35 35 0],[255/255 102/255 255/255],'EdgeColor','none')
alpha(0.3)
plot([td+sigmad td+sigmad],[0 35],':k','linewidth',2)
set(gca,'xticklabel',[])
plot([temp temp],[0 0.5],'-m','linewidth',0.8)
text(temp-0.5,-1.3,'$t_*$','Interpreter','latex','color','magenta','fontsize',22)

text(20,2,'Migration stops','Interpreter','latex','Rotation',90,...
    'color','magenta','fontsize',20)

xlabel('Time, $t$ (a.u.)','Interpreter','latex','fontsize',22)
ylabel('Number of mig-1 mRNA \ molecules, $x$','Interpreter','latex','fontsize',22)
vec_pos = get(get(gca, 'XLabel'), 'Position');
set(get(gca,'xlabel'),'Position',vec_pos+[-2 -2 0])

text(td-0.5,-1.5,'$\bar{t}_d$','Interpreter','latex','fontsize',22)

axis square

yyb(:,1) = timeb;
yyb(:,2) = mrnab;

yyb(:,1) = yyb(:,1)-td+sigmad;

[~,idx1] = sort(abs(yyb(:,1)));
yyb(:,1) = yyb(idx1,1) + td-sigmad;
yyb(:,2) = yyb(idx1,2);
c =1;

for N=2:length(yyb)
         l = 0;
         while l <= N-1
          nxb(l+1) = yyb(length(yyb)-l,1)+td-sigmad;
          nyb(l+1) = yyb(length(yyb)-l,2);
          l = l + 1;
         end
         
    
     P1 = polyfit(nxb,nyb,1);
     slopeyb(c) = P1(1);
     slopexb(c) = N;
     c = c +1;
     
     
      
     
     
end
    
    
x2 = yyb(8,1):0.01:yyb(4,1);
%plot(x2,0.8523*x2 -7.426,'-c','linewidth',2)

%text(10,7,'slope $=1.0 \pm 0.1$','color','green','Interpreter','latex')

yya(:,1) = timea;
yya(:,2) = mrnaa;

yya(:,1) = yya(:,1)-td-sigmad;

[~,idx] = sort(abs(yya(:,1)));
yya(:,1) = yya(idx,1) + td+sigmad;
yya(:,2) = yya(idx,2);


c =1;
for N=3:length(yya)
   
         l = 0;
         while l <= N-1
          nxa(l+1) = yya(l+1,1);
          nya(l+1) = yya(l+1,2);
          l = l + 1;
         end
         
     
    
     P2 = polyfit(nxa,nya,1)
     slopeya(c) = P2(1);
     slopexa(c) = N;
     c = c +1;
    
      
     
end
  
x3 = yya(5,1):0.01:yya(8,1);
%plot(x3,2.1861*x3-21.5252,'-b','linewidth',2)

axis('square')
set(gca,'Fontsize',20)
%print(gcf,'-dpdf', 'fig3a')



%patch([tstar-sigmat tstar-sigmat tstar+sigmat tstar+sigmat ],...
%    [0 35 35 0],[0.6 0.8 1],'EdgeColor','none')

%text(18,36,'$\xleftrightarrow^{Migration stops}$','Interpreter','latex','fontsize',30)
%text(20,32,'slope $=1.9 \pm 0.3$','color','magenta','Interpreter','latex')
%text(5,17.5,'$\updownarrow$','interpreter','latex','fontsize',45)
%text(6.5,17.5,'$\Delta x$','interpreter','latex','fontsize',22)

%print(gcf,'-dpdf', 'n12')



mean(slopeyb(4:14))
mean(slopeya(4:14))


sqrt(std(slopeyb(4:14))^2+std(slopeya(4:14))^2)
mean(slopeya(4:14))-mean(slopeyb(4:14))










