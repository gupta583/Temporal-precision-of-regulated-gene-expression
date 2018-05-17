function mod=stochatic(const,variable,alpha,ensemble,v_par)


nstar = const(1); tstar = const(2); tmax = const(3); H = const(4);
td = const(7); sigma = const(6);
Nr = const(5); r0 = variable(1); k = variable(2);


clear mod;
for ee=1:ensemble
        
        division_switch = 1;
        r0 = variable(1);
        tsplit = normrnd(td*tstar,sigma*tstar);
        
        i = 1;
        t_(i) = 0; x_(i) = 0; r_(i) = Nr; 

        while t_(i) < tmax
             rn = rand(2);
             x_prop = alpha*r0^H/(r_(i)^H+r0^H);
             r_prop = k*r_(i);
             tot_prop = x_prop + r_prop;
             tau = log(1/rn(1))/tot_prop;
             
             
             
             
             
             if division_switch ==1 && t_(i) >= tsplit && tsplit > 0
                division_switch = 0;
                
                count = 0;
                for c=1:x_(i)
                   
                    if rand(1) <= 0.5
                        count = count +1;
                    end
                end
                x_(i) = count;
                 
                count = 0;
                for c=1:r_(i)
                    if rand(1) <= 0.5
                        count = count +1;
                    end
                end
                 r_(i) = count;
                 
                 t_(i) = tsplit;
                 r0 = 0.5*r0;
             end
                 
             if rn(2) < x_prop/tot_prop
                x_(i+1) = x_(i) +1;
                r_(i+1) = r_(i);
             else
                x_(i+1) = x_(i);
                r_(i+1) = r_(i)-1;
             end
                t_(i+1) = t_(i) + tau;
                
                i = i +1;
        end
     
        mod(1:length(x_),ee,1) = x_;
        mod(1:length(r_),ee,2) = r_;
        mod(1:length(t_),ee,3) = t_;
       
       
        clear x_ r_ t_;
end


end