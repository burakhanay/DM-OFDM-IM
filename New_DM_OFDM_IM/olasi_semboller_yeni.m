% olasi sembolleri oluþturma
A_cons=[-0.7+1.3i -0.7-0.7i 1.3+1.3i 1.3-0.7i];
B_cons=[-1.3+0.7i -1.3-1.3i 0.7+0.7i 0.7-1.3i];

     
aa=1;
olasi_semboller00=zeros(length(A_cons)^length(A_cons),length(A_cons));
for ii=1:4
    for jj=1:4
        for kk=1:4
            for ll=1:4
                olasi_semboller00(aa,:)=[A_cons(kk) A_cons(ll) B_cons(jj) B_cons(ii)];                
                aa=aa+1;
            end
        end
    end
end

 
olasi_semboller01=[olasi_semboller00(:,3) olasi_semboller00(:,1) olasi_semboller00(:,4) olasi_semboller00(:,2)];
olasi_semboller10=[olasi_semboller00(:,1) olasi_semboller00(:,3) olasi_semboller00(:,2) olasi_semboller00(:,4)];
olasi_semboller11=[olasi_semboller00(:,3) olasi_semboller00(:,4) olasi_semboller00(:,1) olasi_semboller00(:,2)];
olasi_semboller=[olasi_semboller00;olasi_semboller01;olasi_semboller10;olasi_semboller11];
crit=olasi_semboller;
save('DM_OFDM_ML_tablosu','crit');