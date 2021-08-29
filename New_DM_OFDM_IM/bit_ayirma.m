function [g1_grubu,Ma_grubu,Mb_grubu] = bit_ayirma(bilgi_grup,n,g1,g2,nBitPerSymbol,g)
    
    g1_grubu = zeros(g,g1);
    Ma_grubu = zeros(g,(g2/2));
    Mb_grubu = zeros(g,(g2/2));
    for a=1:g
        for b=1:g1
            g1_grubu(a,b) = bilgi_grup(a,b);
        end
    end
    for a=1:g
        for b=(g1+1):(g1+n)
            Ma_grubu(a,b-2) = bilgi_grup(a,b);
        end
    end
    for a=1:g
        for b = (g1+n+1):nBitPerSymbol
            Mb_grubu(a,b-6) = bilgi_grup(a,b);
        end
    end

end