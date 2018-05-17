function mod=stochatic(const,variable,alpha,ensemble)


 
nstar = const(1); tstar = const(2); tmax = const(3); H = const(4);
Nr = const(5); r0 = variable(1); k = variable(2);

clear mod;
for ee=1:ensemble


        i = 1;
        t_(i) = 0; x_(i) = 0; r_(i) = poissrnd(Nr);
        while t_(i) < tmax
             rn = rand(2);
             x_prop = alpha*r0^H/(r_(i)^H+r0^H);
             r_prop = k*r_(i);
             tot_prop = x_prop + r_prop;
             tau = log(1/rn(1))/tot_prop;

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