function [X_mean,R_mean,T_mean] =meantraj(x_,r_,t_,t0,dt,ensemble)

  dt_ = 0.25;

   ll=1;
   
   for t0=t0:dt_:t0+dt
       
       X = [];
       R = [];
       for j=1:ensemble
           i =1;
           while t_(i,j) <=t0
               i = i+1;
           end
           if i ==1
               X = [X; x_(i,j)];
               R = [R; r_(i,j)];
           else
               X = [X; x_(i-1,j)];
               R = [R; r_(i-1,j)];
           end

       end

     
       x_mean = 0;
       r_mean = 0;
      
       x_edges = 0:1:max(X) +1;
       [N_x,edges1_x] = histcounts(X,x_edges);
       
       r_edges = 0:1:max(R) +1;
       [N_r,edges1_r] = histcounts(R,r_edges);

          for k =1:length(N_x)
             x_mean = x_mean + edges1_x(k)*N_x(k)/length(X);
          end

          for k =1:length(N_r)
             r_mean = r_mean + edges1_r(k)*N_r(k)/length(R);
          end

       X_mean(ll) = x_mean;
       R_mean(ll) = r_mean;
       T_mean(ll) = t0;
      
       ll = ll+1;

    end

end