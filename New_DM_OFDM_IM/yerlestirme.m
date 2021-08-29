function s = yerlestirme(Ma_harita,Mb_harita,g1_grubu,g,n)
   
    s = zeros(g,n);
    for a=1:g
        if isequal(g1_grubu(a,:),[0,0])
            s(a,1) = Ma_harita(a,1);
            s(a,2) = Ma_harita(a,2);
            s(a,3) = Mb_harita(a,1);
            s(a,4) = Mb_harita(a,2);
        elseif isequal(g1_grubu(a,:),[0,1])
            s(a,1) = Mb_harita(a,1);
            s(a,2) = Ma_harita(a,1);
            s(a,3) = Mb_harita(a,2);
            s(a,4) = Ma_harita(a,2);
        elseif isequal(g1_grubu(a,:),[1,0])
            s(a,1) = Ma_harita(a,1);
            s(a,2) = Mb_harita(a,1);
            s(a,3) = Ma_harita(a,2);
            s(a,4) = Mb_harita(a,2);
        else
            s(a,1) = Mb_harita(a,1);
            s(a,2) = Mb_harita(a,2);
            s(a,3) = Ma_harita(a,1);
            s(a,4) = Ma_harita(a,2);
        end
    end

end