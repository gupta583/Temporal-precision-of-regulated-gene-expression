clear all;

nstar1 = 15;
tstar1 = 20;
M = xlsread('exp_mig.xlsx');  %experimental data
shift = 3.2;


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

set(gca,'Fontsize',20)











