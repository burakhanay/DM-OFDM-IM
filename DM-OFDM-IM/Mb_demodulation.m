function Mb_alinan = Mb_demodulation(rx_signal,g,n)
    
    Mb_alinan = zeros(g,n);
    c = 1;
    for a=1:g
        temp_alici = rx_signal(a,:);
        for b = 1:n
            if  abs(temp_alici(1,b))<=2
                temp_alici(1,b) = 0;
            else 
                temp_alici(1,b) = temp_alici(1,b);
            end
        end     
        for b = 1:n       
            if abs(temp_alici(1,b)) > 2
               if single(real(temp_alici(1,b))) == single(2.7321)
                   Mb_alinan(a,c) = 0;
                   Mb_alinan(a,c+1) = 0;
               elseif single(imag(temp_alici(1,b))) == single(2.7321)
                   Mb_alinan(a,c) = 0;
                   Mb_alinan(a,c+1) = 1;
               elseif single(imag(temp_alici(1,b))) == single(-2.7321)
                   Mb_alinan(a,c) = 1;
                   Mb_alinan(a,c+1) = 0;
               elseif single(real(temp_alici(1,b))) == single(-2.7321)
                   Mb_alinan(a,c) = 1;
                   Mb_alinan(a,c+1) = 1;    
               end
                c = c+2;
              if c==5
                  c=1;
              end
            end
        end
    end
    
end