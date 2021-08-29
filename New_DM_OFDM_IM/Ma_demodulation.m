function Ma_alinan = Ma_demodulation(rx_signal,g,n)
    
    Ma_alinan = zeros(g,n);
    c = 1;
    for a=1:g
        temp_alici = rx_signal(a,:);
        for b = 1:n
            if  real(temp_alici(1,b))==0.7 || real(temp_alici(1,b))==-1.3
                temp_alici(1,b) = 0;
            else 
                temp_alici(1,b) = temp_alici(1,b);
            end
        end 
        for b = 1:n       
            if abs(temp_alici(1,b)) > 0
               if (temp_alici(1,b)) == (-0.7+1.3i)
                   Ma_alinan(a,c) = 0;
                   Ma_alinan(a,c+1) = 0;
               elseif (temp_alici(1,b)) == (-0.7-0.7i)
                   Ma_alinan(a,c) = 0;
                   Ma_alinan(a,c+1) = 1;
               elseif (temp_alici(1,b)) == (1.3+1.3i)
                   Ma_alinan(a,c) = 1;
                   Ma_alinan(a,c+1) = 0;
               elseif (temp_alici(1,b)) == (1.3-0.7i)
                   Ma_alinan(a,c) = 1;
                   Ma_alinan(a,c+1) = 1;    
               end
               c = c+2;
               if c==5
                  c=1;
               end
            end
        end
    end

end