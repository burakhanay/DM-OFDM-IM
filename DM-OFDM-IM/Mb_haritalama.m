function Mb_harita = Mb_haritalama(Mb_grubu,temp_ilk2,temp_son2,g)
    Mb_harita = zeros(g,2);
    for a=1:g
                temp_ilk2 = [Mb_grubu(a,1) Mb_grubu(a,2)];
                temp_son2 = [Mb_grubu(a,3) Mb_grubu(a,4)];
                if isequal(temp_ilk2,[0,0])
                    Mb_harita(a,1) = single(2.7321);
                elseif isequal(temp_ilk2,[0,1])
                    Mb_harita(a,1) = single(2.7321i);
                elseif isequal(temp_ilk2,[1,0])
                    Mb_harita(a,1) = single(-2.7321i); 
                else
                    Mb_harita(a,1) = single(-2.7321);   
                end
                if isequal(temp_son2,[0,0])
                    Mb_harita(a,2) = single(2.7321);
                elseif isequal(temp_son2,[0,1])
                    Mb_harita(a,2) = single(2.7321i);
                elseif isequal(temp_son2,[1,0])
                    Mb_harita(a,2) = single(-2.7321i); 
                else
                    Mb_harita(a,2) = single(-2.7321);   
                end
    end
    
end