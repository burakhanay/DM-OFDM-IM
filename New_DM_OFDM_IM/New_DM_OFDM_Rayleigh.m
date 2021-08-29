clear all;
close all;
clc;

load('DM_OFDM_ML_tablosu.mat');
olabilir=zeros(1,length(crit));

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
nSym = 10^2;
nTap = 10;

g1_grubu = zeros(g,g1);
Ma_grubu = zeros(g,(g2/2));
Mb_grubu = zeros(g,(g2/2));

temp_tx_ilk2 = zeros(1,2);
temp_tx_son2 = zeros(1,2);
Ma_harita = zeros(nSym,2);
Mb_harita = zeros(nSym,2);

xF = zeros(nSym,nFFTSize);
xT = zeros(nSym,nFFTSize);
xhT = zeros(nSym,(nFFTSize+cyclic_prefix+nTap-1));
yF = zeros(nSym,m);
yF_1D = zeros(1,nSym*m);
demod = zeros(nSym,m);

temp_alici = zeros(g,n);
g1_alinan = zeros(g,k);
Ma_alinan = zeros(g,n);
Mb_alinan = zeros(g,n);

SE = m/(cyclic_prefix + nFFTSize);
Eb = 0:5:30;  %SNR degerleri
EsN0dB = Eb + 10*log10(nFFTSize/(nFFTSize+cyclic_prefix))+ 10*log10(SE);
%EsN0dB = Eb + 10*log10(SE);

nErr = zeros(1,length(Eb));
simBer = zeros(1,length(Eb));
%M_ary=4;
%olcek = sqrt(1/((2/3)*(M_ary-1)));

for cc=1:length(Eb)
    bilgi = randi([0 1],1,(nSym*m)); %0 ve 1 lerden olusan bilgi biti olusturuldu
    %eb_enerjisi=(norm(bilgi)^2)/length(bilgi);%anlamlý bitin enerjisi hesaplandý
    bilgi_grup = reshape(bilgi.',m,nSym).'; %bilgi biti gruplara bolundu
    for aa=1:nSym
        temp_verici1 = reshape(bilgi_grup(aa,:)',nBitPerSymbol,g).'; %bilgi biti gruplara bolundu
        [g1_grubu,Ma_grubu,Mb_grubu] = bit_ayirma(temp_verici1,n,g1,g2,nBitPerSymbol,g);
        %1 ve 2. bitler g1//3-6. bitler Ma_grubuna//7-10. bitler Mb_grubuna atandi
        Ma_harita = Ma_haritalama(Ma_grubu,temp_tx_ilk2,temp_tx_son2,g);
        Mb_harita = Mb_haritalama(Mb_grubu,temp_tx_ilk2,temp_tx_son2,g);
        temp_verici_2 = yerlestirme(Ma_harita,Mb_harita,g1_grubu,g,n);
        temp_verici_2 = reshape(temp_verici_2.',1,[]);
        xF(aa,:) = temp_verici_2;
        xT(aa,:) = (nFFTSize/sqrt(nFFTSize))*ifft(fftshift(xF(aa,:)),nFFTSize);
        %(norm(xT(aa,:))^2)/length(xT(aa,:))
        olcek(aa,1) = sqrt(mean(abs(xT(aa,:)).^2));
        xT(aa,:) = xT(aa,:)./olcek(aa,1);
        %xT_enerji(aa,1)=(norm(xT(aa,:))^2)/length(xT(aa,:))
    end
    xT_cyclic = [xT(:,(nFFTSize-cyclic_prefix+1):(nFFTSize)) xT];%cyclic-prefix eklendi
    
%     for aa=1:nSym   
%     %(norm(xF(aa,:))^2)/320%constellation enerji
%     es_enerjisi(1,aa) = (norm(xT_cyclic(aa,:)).^2)/length(xT_cyclic(aa,:));
%     end
%     ortalama_es=mean(es_enerjisi);%cyclic eklenmiþ zaman uzayýndaki sinyalin enerjisi
%     kod_oran  = eb_enerjisi/ortalama_es
%     eb_es = (10^(Eb(1,1)/10)) / (10^(EsN0dB(1,1)/10))  
    %%%%%Rayleigh Kanalý%%%%%%
    [ht,hF]=Ray_kanal(nFFTSize,nSym,nTap);
    %%%%%Rayleigh Kanalý%%%%%%
    
    for aa = 1:nSym%ht ile konvole edildi
        xhT(aa,:) = conv(xT_cyclic(aa,:),ht(aa,:));
    end
    xT_1D = reshape(xhT.',1,[]);%1x(nSym*m) uzunluðunda diziye donusturuldu
    
    %%%%noise%%%%
    snr=EsN0dB(cc);
    boyut = nSym*(cyclic_prefix+nFFTSize);
    sigma = sqrt((boyut/(nSym*m))/(10^(snr/10)));% Eb/N0
    % total tx power/total # of inf. bits   Eb=numel(xt)/numel(bitler)
    %boyut = length(xT_1D);
    nt = sqrt(1/2)*(sigma).*(randn(1,length(xT_1D)) + 1i*randn(1,length(xT_1D)));
    %(norm(nt).^2)/length(nt)
    %%%%noise%%%%
    outputIFFT = sqrt((nFFTSize+cyclic_prefix)/nFFTSize)*xT_1D + nt;
    
    %%%%%%%%%%%%%
    %%%%ALICI%%%%
    %%%%%%%%%%%%%
    
    rx_signal = reshape(outputIFFT,(nFFTSize+cyclic_prefix+nTap-1),nSym).';
    rx_signal = rx_signal(:,(cyclic_prefix+1:nFFTSize+cyclic_prefix)); %cyclic-prefix kaldirildi
    for aa = 1:nSym
        rx_signal(aa,:) = (sqrt(nFFTSize)/nFFTSize)*fftshift(fft(rx_signal(aa,:),nFFTSize));
        rx_signal(aa,:) = rx_signal(aa,:).*olcek(aa,1);
        %rx_signal(aa,:) = rx_signal(aa,:) ./ hF(aa,:);
    end
    for aa=1:nSym %rx_signal'in her bir satýrý tek tek ele alýnacaktýr
        temp_alici = reshape(rx_signal(aa,:),n,g).';%alýnan dizi 32x4'e cevrilip
        %temp_alici'ya atandý,burda her bir 4'lü karþýlaþtýrýlacaktýr
        hF_bloklari=reshape(hF(aa,:).',n,g).';
        %ilgili hF satýrý 32x4'lük matrise dönüþtürüldü,her 4'lü crit
        %ile çarpýlýp karþýlaþtýrma yapýlcaktýr.
        for a=1:g
            for j=1:length(crit)
                olabilir(j) = sum(abs(temp_alici(a,:)-hF_bloklari(a,:).*crit(j,:)).^2,2);
                if olabilir(j)==0%eðer %100 benzerlik varsa diðer bloklar kontrol edilmeyecek
                    break;
                end
            end
            [PP,I] = min(olabilir);
            h=I(1,1);
            temp_alici(a,:)=crit(h,:);
        end
              
        g1_alinan = g1_demodulation(temp_alici,g,k,n);%%g1 bitleri elde edildi
        Ma_alinan = Ma_demodulation(temp_alici,g,n);%%Ma bitleri elde edildi
        Mb_alinan = Mb_demodulation(temp_alici,g,n);%%Mb bitleri elde edildi
        
        demod = [g1_alinan Ma_alinan Mb_alinan];
        yF(aa,:) = reshape(demod.',1,[]);
        
    end
    yF_1D = reshape(yF.',1,[]);
    nErr(cc) = size(find(bilgi-yF_1D),2);
    simBer(cc) = nErr(cc)/(nSym*m);
end


DMOFDM_IM_makale = [0.15 0.07 0.025 0.007 2*1e-3 5e-4 0.15*1e-3];
figure,semilogy(Eb,simBer,'--s'),grid,hold on
semilogy(0:5:30,DMOFDM_IM_makale,'--rx'),
axis([0 30 10^-4 1])
grid on
legend('Benzetim','Makale Sonucu [6]');
xlabel(' Eb/No, dB ')
ylabel(' Bit Error Rate ')
%title(' BER for new constellation DM-OFDM-IM ')