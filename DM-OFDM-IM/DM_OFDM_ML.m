clear all;
close all;
clc;
x=[4096,4];
h=1;
k=1;
crit = zeros(1536,4);
for a=1:4096
    x(a,4)=h;
    h=h+1;
    k=k+1;
    if k==9
        h=1;
        k=1;
    end   
end
h=1;
k=1;
for a=1:4096
    x(a,3)=h;
    k=k+1;
    if rem(k,8)==1
        h=h+1;
    end   
    if h==9
        h=1;
        k=1;  
    end
end
h=1;
k=1;
for a=1:4096
    x(a,2)=h;
    k=k+1;
    if rem(k,64)==1
        h=h+1;
    end   
    if h==9
        h=1;
        k=1;  
    end
end

h=1;
k=1;
for a=1:4096
    x(a,1)=h;
    k=k+1;
    if rem(k,512)==1
        h=h+1;
    end   
    if h==9
        h=1;
        k=1;  
    end
end

for a=1:4096
    for b=1:4 
        if x(a,b)==1
            x(a,b)=single(1-1i);
        elseif x(a,b)==2
            x(a,b)=single(1+1i);
        elseif x(a,b)==3
            x(a,b)=single(-1-1i);
        elseif x(a,b)==4
            x(a,b)=single(-1+1i);   
        elseif x(a,b)==5
            x(a,b)=single(2.7321);
        elseif x(a,b)==6
            x(a,b)=single(2.7321i);
        elseif x(a,b)==7
            x(a,b)=single(-2.7321i);
        else
            x(a,b)=single(-2.7321);
        end
    end
end
temp=0;
kk=1;
for a=1:4096
    temp=0;
    for b=1:4
        if abs(x(a,b)) > 2
            temp = temp+1;
        end
    end
    if temp==2
        crit(kk,:)=x(a,:);
        kk=kk+1;
    end
end

save('DM_OFDM_ML_tablosu','crit');
