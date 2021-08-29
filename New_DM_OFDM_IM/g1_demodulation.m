function g1_alinan = g1_demodulation(rx_signal,g,k,n)
     
     g1_alinan = zeros(g,k);
     for a = 1:g
        temp_alici = rx_signal(a,:);
        for b = 1:n
            if  real(temp_alici(1,b))==0.7 || real(temp_alici(1,b))==-1.3
                temp_alici(1,b) = 0;
            else 
                temp_alici(1,b) = 1;
            end
        end   

        if (temp_alici(1,1) == 1 && temp_alici(1,2) == 1)  
            g1_alinan(a,:) = [0,0];
        elseif (temp_alici(1,2) == 1 && temp_alici(1,4) == 1)  
            g1_alinan(a,:) = [0,1];
        elseif (temp_alici(1,1) == 1 && temp_alici(1,3) == 1)  
            g1_alinan(a,:) = [1,0];
        elseif (temp_alici(1,3) == 1 && temp_alici(1,4) == 1)  
            g1_alinan(a,:) = [1,1];
        end
    end



end