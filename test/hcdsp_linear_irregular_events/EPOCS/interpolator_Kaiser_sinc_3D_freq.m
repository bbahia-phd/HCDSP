function [d_rec] = interpolator_Kaiser_sinc_3D_freq(U,x_irreg,x_reg,N,a,type)
%  type  == 0: from regular to irregular
%  type  == 1: from irregular to regular


dh1 = x_reg(2,1,1) - x_reg(1,1,1);
dh2 = x_reg(1,2,2) - x_reg(1,1,2);

[nh1,nh2,~] = size(x_reg);

nk = size(x_irreg,1);


if type ==0  % from regular to irregular
   
    % %     input is a 3D data (nt,nx,ny); 
    % %     Output is 2D data (nt,nk)     
    
    d_rec = zeros(nk,1);

    
    for k = 1:nk
        
        ia = floor((x_irreg(k,1)-x_reg(1,1,1))/dh1) + 1;        
        ja = floor((x_irreg(k,2)-x_reg(1,1,2))/dh2) + 1; 
        
% %   define the boundary for ia and ja
        min_ia = max(ia-N,1);
        max_ia = min(ia+N,nh1);
        
        min_ja = max(ja-N,1);
        max_ja = min(ja+N,nh2);
       
        for ia = min_ia:1:max_ia 
            for ja = min_ja:1:max_ja
                
                t = (x_irreg(k,1) - x_reg(ia,ja,1))/dh1;
                u = (x_irreg(k,2) - x_reg(ia,ja,2))/dh2;
                
                Wt = Kaiser_sinc(t,N+1,a);
                Wu = Kaiser_sinc(u,N+1,a);
                d_rec(k) = d_rec(k)+ Wt*Wu*U(ia,ja);

            end
        end
        
    end
    
else % %==1: from irregular to regular
    
    % %     input is a 1D data (nk,1) for each freq.; 
    % %     Output is 2D data (nh1,nh2) for each freq.;
    
    d_rec = zeros(nh1,nh2);

    for k = 1:nk
        ia = floor((x_irreg(k,1)-x_reg(1,1,1))/dh1) + 1;        
        ja = floor((x_irreg(k,2)-x_reg(1,1,2))/dh2) + 1; 
        
        min_ia = max(ia-N,1);
        max_ia = min(ia+N,nh1);
        
        min_ja = max(ja-N,1);
        max_ja = min(ja+N,nh2);
        
        
        for ia = min_ia:1:max_ia 
            for ja = min_ja:1:max_ja
                
                t = (x_irreg(k,1) - x_reg(ia,ja,1))/dh1;
                u = (x_irreg(k,2) - x_reg(ia,ja,2))/dh2;
                
                Wt = Kaiser_sinc(t,N+1,a);
                Wu = Kaiser_sinc(u,N+1,a);
                
                d_rec(ia,ja) = d_rec(ia,ja) +  Wt*Wu*U(k);
              
            end
        end
        

        
    
    end

end
