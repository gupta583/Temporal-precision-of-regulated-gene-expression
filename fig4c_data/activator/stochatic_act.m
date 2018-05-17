function mod=stochatic_act(const,variable,alpha,ensemble)


 
nstar = const(1); tstar = const(2); tmax = const(3); H = const(4);
td = const(7); sigma = const(6);
Na = const(5); a0 = variable(1); k = variable(2);

clear mod;
for ee=1:ensemble

        
        division_switch = 1;
        tsplit = normrnd(td*tstar,sigma*tstar);
     
        a0 = variable(1);
        k = variable(2);
        i = 1;
        t_(i) = 0; x_(i) = 0; a_(i) = Na;
        
        
        while t_(i) < tmax
             rn = rand(2,1);
             x_prop = alpha*a_(i)^H/(a_(i)^H+a0^H);
             a_prop = k;
             tot_prop = x_prop + a_prop;
             tau = log(1/rn(1))/tot_prop;
             
             
              
             
              if division_switch ==1 && t_(i) >= tsplit && tsplit > 0
                division_switch = 0;
                
                count = 0;
                
                for c=1:x_(i)
                    if rand(1) <= 0.5;
                        count = count +1;
                    end
                end
                x_(i) = count;
                 
                count = 0;
                for c=1:a_(i)
                    if rand(1) <= 0.5;
                        count = count +1;
                    end
                end
                 a_(i) = count;
                 
                 t_(i) = tsplit;
                 a0 = 0.5*a0;
             end
             
             if rn(2) < x_prop/tot_prop
                x_(i+1) = x_(i) +1;
                a_(i+1) = a_(i);
             else
                x_(i+1) = x_(i);
                a_(i+1) = a_(i)+1;
             end
                t_(i+1) = t_(i) + tau;
                i = i +1;
        end
     
        mod(1:length(x_),ee,1) = x_;
        mod(1:length(a_),ee,2) = a_;
        mod(1:length(t_),ee,3) = t_;

        clear x_ a_ t_;
end


end