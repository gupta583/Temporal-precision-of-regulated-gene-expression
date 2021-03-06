function [X_mean,A_mean,T_mean] =meantraj_act(x_,a_,t_,t0,dt,ensemble)

   ll=1;
   dt_ = 0.25;
   
      
   for t0=t0:dt_:t0+dt
       
       X = [];
       A = [];
       for j=1:ensemble
           i =1;
           while t_(i,j) <=t0
               i = i+1;
           end
           if i ==1
               X = [X; x_(i,j)];
               A = [A; a_(i,j)];
           else
               X = [X; x_(i-1,j)];
               A = [A; a_(i-1,j)];
           end

       end

     
       x_mean = 0;
       a_mean = 0;
      
       x_edges = 0:1:max(X) +1;
       [N_x,edges1_x] = histcounts(X,x_edges);
       
       a_edges = 0:1:max(A) +1;
       [N_a,edges1_a] = histcounts(A,a_edges);

          for k =1:length(N_x)
             x_mean = x_mean + edges1_x(k)*N_x(k)/length(X);
          end

          for k =1:length(N_a)
             a_mean = a_mean + edges1_a(k)*N_a(k)/length(A);
          end

       X_mean(ll) = x_mean;
       A_mean(ll) = a_mean;
       T_mean(ll) = t0;
      
       ll = ll+1;
       clear N_x N_a x_edges edges1_x a_edges edges1_a

    end

end