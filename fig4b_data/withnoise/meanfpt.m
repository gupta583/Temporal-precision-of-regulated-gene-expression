function [fpt, V_fpt]  = meanfpt(mod,nstar,ensemble)


      x_(:,:) = mod(:,:,1);
      r_(:,:) = mod(:,:,2);
      t_(:,:) = mod(:,:,3);
        

       T = [];
       fpt = 0;
       V_fpt = 0;
       for j=1:ensemble
           for i =1:length(x_(:,j))
               if x_(i,j) == nstar
                  T = [T; t_(i,j)];
               end
           end
       end

      [N_T,edges_T] = histcounts(T,3*length(T)+1);

      for k =1:length(N_T)
             fpt = fpt + edges_T(k)*N_T(k)/length(T);
             V_fpt = V_fpt + edges_T(k)^2*N_T(k)/length(T);
      end

     V_fpt = V_fpt - fpt^2;
     
end
