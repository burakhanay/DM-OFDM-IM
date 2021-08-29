clear all;
close all;


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
nSym = 10^1;

g1_grubu = zeros(g,g1);
Ma_grubu = zeros(g,(g2/2));
Mb_grubu = zeros(g,(g2/2));

temp_ilk2 = zeros(1,2);
temp_son2 = zeros(1,2);
Ma_harita = zeros(nSym,2);
Mb_harita = zeros(nSym,2);

xF = zeros(nSym,nFFTSize);
xT = zeros(nSym,nFFTSize);
yF = zeros(nSym,m);
yF_1D = zeros(1,nSym*m);
demod = zeros(nSym,m);

temp_alici_2 = zeros(1,n);
g1_alinan = zeros(g,k);
Ma_alinan = zeros(g,n);
Mb_alinan = zeros(g,n);


Eb = 0:2:12;  %SNR degerleri
nTap = 10;
EsN0dB  = Eb + 10*log10(2.5) + 10*log10(128/144);
nErr = zeros(1,length(Eb));
enerji_esitleyen = zeros(nSym,1);

for cc=4:length(Eb)
    
        sigma = sqrt(1/(10^(EsN0dB(cc)/10))); 
        nt = sigma*sqrt(1/2)*(randn(1,(nSym*(nFFTSize+cyclic_prefix))) + 1i*randn(1,(nSym*(nFFTSize+cyclic_prefix))));     
       
        bilgi = randi([0 1],1,(nSym*m)); %0 ve 1 lerden olusan bilgi biti olusturuldu
        bilgi_grup = reshape(bilgi.',m,nSym).'; %bilgi biti gruplara bolundu
        for aa=1:nSym
            temp_verici1 = reshape(bilgi_grup(aa,:)',nBitPerSymbol,g).'; %bilgi biti gruplara bolundu
            [g1_grubu,Ma_grubu,Mb_grubu] = bit_ayirma(temp_verici1,n,g1,g2,nBitPerSymbol,g);
            %1 ve 2. bitler g1//3-6. bitler Ma_grubuna//7-10. bitler Mb_grubuna atandi
            Ma_harita = Ma_haritalama(Ma_grubu,temp_ilk2,temp_son2,g);  
            Mb_harita = Mb_haritalama(Mb_grubu,temp_ilk2,temp_son2,g);
            temp_verici_2 = yerlestirme(Ma_harita,Mb_harita,g1_grubu,g,n);
            temp_verici_2 = reshape(temp_verici_2.',1,[]);
            enerji_esitleyen(aa,1) = 1/sqrt(mean(abs(temp_verici_2).^2));
            xF(aa,:) = enerji_esitleyen(aa,1)*temp_verici_2;
            xT(aa,:) = (nFFTSize/sqrt(nFFTSize))*ifft(fftshift(xF(aa,:)),nFFTSize);           
        end
       
        xT_cyclic = [xT(:,(nFFTSize-cyclic_prefix+1):(nFFTSize)) xT];  %cyclic-prefix eklendi 
        xT_1D = reshape(xT_cyclic.',1,[]);
        
        %noise = addGaussianNoise(xT_1D,EbN0dB);
        %outputIFFT = sqrt(144/128)*xT_1D + noise;
        %outputIFFT = sqrt(144/128)*xT_1D + 10^(-EsN0dB(EbN0dB)/20)*nt;
        outputIFFT = sqrt(144/128)*xT_1D + nt;
        
        %%%%%%%%%%%%%
        %%%%ALICI%%%%
        %%%%%%%%%%%%%

        rx_signal = reshape(outputIFFT,(nFFTSize+cyclic_prefix),nSym).';
        rx_signal = rx_signal(:,(cyclic_prefix+1):(nFFTSize+cyclic_prefix)); %cyclic-prefix kaldirildi
        for aa=1:nSym   
            rx_signal(aa,:) = (sqrt(nFFTSize)/nFFTSize)*fftshift(fft(rx_signal(aa,:),nFFTSize));
            rx_signal(aa,:) = rx_signal(aa,:) / enerji_esitleyen(aa,1);
            temp_alici_1 = (reshape(rx_signal(aa,:),n,g).');
            for a=1:g
                temp_alici_2=temp_alici_1(a,:);
                for j=1:length(crit)
                    olabilir(j) = sum(abs((temp_alici_2)-(crit(j,:))).^2,2);
                    if olabilir(j)==0
                        break;
                    end
                end
                [PP,I] = min(olabilir);
                h=I(1,1);
                temp_alici_1(a,:)=crit(h,:);
            end

            g1_alinan = g1_demodulation(temp_alici_1,g,k,n);%%g1 bitleri elde edildi
            Ma_alinan = Ma_demodulation(temp_alici_1,g,n);%%Ma bitleri elde edildi
            Mb_alinan = Mb_demodulation(temp_alici_1,g,n);%%Mb bitleri elde edildi

            demod = [g1_alinan Ma_alinan Mb_alinan];
            yF(aa,:) = reshape(demod.',1,[]);
        end
        
        yF_1D = reshape(yF.',1,[]);
        nErr(cc) = size(find(bilgi-yF_1D),2);
        
end


simBer = nErr/(nSym*m);
semilogy(Eb,simBer,'mx-','LineWidth',2);
axis([0 30 10^-6 1])
grid on
legend(' AWGN');
xlabel(' Eb/No, dB ')
ylabel(' Bit Error Rate ')
title(' BER for new constellation DM-OFDM in AWGN channel ')

