function Ma_harita = Ma_haritalama(Ma_grubu,temp_ilk2,temp_son2,g)
    Ma_harita = zeros(g,2);    
    for a=1:g
                
                temp_ilk2 = [Ma_grubu(a,1) Ma_grubu(a,2)];
                temp_son2 = [Ma_grubu(a,3) Ma_grubu(a,4)];
                if isequal(temp_ilk2,[0,0])
                    Ma_harita(a,1) = (-0.7+1.3i);
                elseif isequal(temp_ilk2,[0,1])
                    Ma_harita(a,1) = (-0.7-0.7i);
                elseif isequal(temp_ilk2,[1,0])
                    Ma_harita(a,1) = (1.3+1.3i); 
                else
                    Ma_harita(a,1) = (1.3-0.7i);   
                end
                if isequal(temp_son2,[0,0])
                    Ma_harita(a,2) = (-0.7+1.3i);
                elseif isequal(temp_son2,[0,1])
                    Ma_harita(a,2) = (-0.7-0.7i);
                elseif isequal(temp_son2,[1,0])
                    Ma_harita(a,2) = (1.3+1.3i); 
                else
                    Ma_harita(a,2) = (1.3-0.7i);   
                end         
    end
end
