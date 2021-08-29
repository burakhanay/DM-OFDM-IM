clear all;
close all;
clc;
x=[4096,4];
h=1;
k=1;
crit = [];
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
            x(a,b)=(1.3+1.3i);
        elseif x(a,b)==2
            x(a,b)=(1.3-0.7i);
        elseif x(a,b)==3
            x(a,b)=(-0.7+1.3i);
        elseif x(a,b)==4
            x(a,b)=(-0.7-0.7i);   
        elseif x(a,b)==5
            x(a,b)=(0.7+0.7i);
        elseif x(a,b)==6
            x(a,b)=(0.7-1.3i);
        elseif x(a,b)==7
            x(a,b)=(-1.3+0.7i);
        else
            x(a,b)=(-1.3-1.3i);
        end
    end
end
temp_a=0;
temp_b=0;
kk=1;
for a=1:4096
    temp_a=0;
    temp_b=0;
    for b=1:4
        if real(x(a,b))==1.3 || real(x(a,b))==-0.7
            temp_a = temp_a+1;
        end
    end
    for b=1:4
        if real(x(a,b))==-1.3 || real(x(a,b))==0.7
            temp_b = temp_b+1;
        end
    end
    if temp_b==2
        crit(kk,:)=x(a,:);
        kk=kk+1;
    end
    if temp_a==2
        crit(kk,:)=x(a,:);
        kk=kk+1;
    end
 end

save('DM_OFDM_ML_tablosu','crit');
