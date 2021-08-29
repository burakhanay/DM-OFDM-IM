clear all;
close all;
clc;

load('DM_OFDM_ML_tablosu.mat');
olabilir=[1,length(crit)];

cyclic_prefix = 16;
nFFTSize = 128;
g = 32;           %grup sayisi
n = nFFTSize/g;   
k = n/2;
g1 = floor(log2( faktoriyel(n)/(faktoriyel(n-k)*faktoriyel(k)) ));
Ma = n;
Mb = n;
g2 = k*log2(Ma) + ((n-k)*log2(Mb));
nBitPerSymbol = g1+g2;
m = nBitPerSymbol*g;    %bilgi biti sayisi

g1_grubu = zeros(g,g1);
Ma_grubu = zeros(g,(g2/2));
Mb_grubu = zeros(g,(g2/2));

temp_ilk2 = zeros(1,2);
temp_son2 = zeros(1,2);
Ma_harita = zeros(g,2);
Mb_harita = zeros(g,2);

s = zeros(g,n);

temp_alici = zeros(1,n);
g1_alinan = zeros(g,k);
Ma_alinan = zeros(g,n);
Mb_alinan = zeros(g,n);
for kk=1:10
    bilgi = randi([0 1],1,m); %0 ve 1 lerden olusan bilgi biti olusturuldu
    bilgi_grup = reshape(bilgi.',nBitPerSymbol,g).'; %bilgi biti gruplara bolundu
    [g1_grubu,Ma_grubu,Mb_grubu] = bit_ayirma(bilgi_grup,n,g1,g2,nBitPerSymbol,g);
    %1 ve 2. bitler g1//3-6. bitler Ma_grubuna//7-10. bitler Mb_grubuna atandi
    
    Ma_harita = Ma_haritalama(Ma_grubu,temp_ilk2,temp_son2,g);  
    Mb_harita = Mb_haritalama(Mb_grubu,temp_ilk2,temp_son2,g);
    s = yerlestirme(Ma_harita,Mb_harita,g1_grubu,g,n);

    s_1D = reshape(s.',1,[]);
    outputIFFT = ifft(fftshift(s_1D),nFFTSize);%carpimmi hocaya sor
    outputIFFT = (nFFTSize/sqrt(n)).*outputIFFT;
    outputIFFT = [outputIFFT(:,(nFFTSize-cyclic_prefix+1):(nFFTSize)) outputIFFT]; %cyclic-prefix eklendi

    %%%%%%%%%%%%%
    %%%%ALICI%%%%
    %%%%%%%%%%%%%

    rx_signal = outputIFFT(:,(cyclic_prefix+1):(nFFTSize+cyclic_prefix)); %cyclic-prefix kaldirildi
    rx_signal = (sqrt(n)/nFFTSize)*fftshift(fft(rx_signal,nFFTSize));
    rx_signal = reshape(rx_signal,n,g).';

    for a=1:32
        temp_alici=rx_signal(a,:);
        for j=1:length(crit)
            olabilir(j) = sum(abs(temp_alici-crit(j,:)).^2,2);
        end
        [PP,I] = min(olabilir);
        h=I(1,1);
        rx_signal(a,:)=crit(h,:);
    end

             %%g1 bitleri elde edilecek
    g1_alinan = g1_demodulation(rx_signal,g,k,n);%%g1 bitleri elde edildi
    Ma_alinan = Ma_demodulation(rx_signal,g,n);%%Ma bitleri elde edildi
    Mb_alinan = Mb_demodulation(rx_signal,g,n);%%Mb bitleri elde edildi

    demod = [g1_alinan Ma_alinan Mb_alinan];
    demod_1D = reshape(demod.',1,[]);
    hata = size(find(bilgi-demod_1D),2)

end
