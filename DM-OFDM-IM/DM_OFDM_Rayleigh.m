clear all;
close all;
clc;

load('DM_OFDM_ML_tablosu.mat');
% crit_2=zeros(length(crit),4);
% for cc=1:length(crit)
%     enerji_esitleyen = 1/sqrt(mean(abs(crit(cc,:)).^2));
%     crit_2(cc,:) = enerji_esitleyen * single(crit(cc,:));
% end

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
nSym = 10^2;
nTap = 10;

g1_grubu = zeros(g,g1);
Ma_grubu = zeros(g,(g2/2));
Mb_grubu = zeros(g,(g2/2));

temp_ilk2 = zeros(1,2);
temp_son2 = zeros(1,2);
Ma_harita = zeros(nSym,2);
Mb_harita = zeros(nSym,2);

xF = zeros(nSym,nFFTSize);
xT = zeros(nSym,nFFTSize);
xhT = zeros(nSym,(nFFTSize+cyclic_prefix+nTap-1));
yF = zeros(nSym,m);
yF_1D = zeros(1,nSym*m);
demod = zeros(nSym,m);

temp_alici_2 = zeros(1,n);
g1_alinan = zeros(g,k);
Ma_alinan = zeros(g,n);
Mb_alinan = zeros(g,n);


Eb = 0:5:30;  %SNR degerleri
SE = m/(cyclic_prefix + nFFTSize);
EsN0dB = Eb + 10*log10(nFFTSize/(nFFTSize+cyclic_prefix))+ 10*log10(SE);

nErr = zeros(1,length(Eb));

nt = zeros(1,nSym*(nFFTSize+cyclic_prefix+nTap-1));
olcek = zeros(nSym,1);
for cc=1:length(Eb)   
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
            xF(aa,:) = temp_verici_2;
            %xF_enerji=(norm(xF(aa,:))^2)/length(xF(aa,:))
            xT(aa,:) = (nFFTSize/sqrt(nFFTSize))*ifft(fftshift(xF(aa,:)),nFFTSize);
            olcek(aa,1) = sqrt(mean(abs(xT(aa,:)).^2));
            xT(aa,:) = xT(aa,:)./olcek(aa,1);
            %xT_enerji=(norm(xT(aa,:))^2)/length(xT(aa,:))
        end
        xT_cyclic = [xT(:,(nFFTSize-cyclic_prefix+1):(nFFTSize)) xT];  %cyclic-prefix eklendi 
        %%%%%Rayleigh Kanalý%%%%%%
        [ht,hF]=Ray_kanal(nFFTSize,nSym,nTap);
        %%%%%Rayleigh Kanalý%%%%%%
        for aa = 1:nSym
            xhT(aa,:) = conv(xT_cyclic(aa,:),ht(aa,:));
        end
        xT_1D = reshape(xhT.',1,[]);
        %%%%noise%%%%
        snr=EsN0dB(cc);
        boyut = nSym*(cyclic_prefix+nFFTSize);
        sigma = sqrt((boyut/(nSym*m))/(10^(snr/10)));% Eb/N0
        %sigma_eski = sqrt(1/10^(snr/10));
        % total tx power/total # of inf. bits   Eb=numel(xt)/numel(bitler)
        %boyut = length(xT_1D);
        nt = sqrt(1/2)*(sigma).*(randn(1,length(xT_1D)) + 1i*randn(1,length(xT_1D)));
        %%%%noise%%%%
        outputIFFT = sqrt((nFFTSize+cyclic_prefix)/nFFTSize)*xT_1D + nt;
        %%%%%%%%%%%%%
        %%%%ALICI%%%%
        %%%%%%%%%%%%%

        rx_signal = reshape(outputIFFT,(nFFTSize+cyclic_prefix+nTap-1),nSym).';
        rx_signal = rx_signal(:,(17:144)); %cyclic-prefix kaldirildi
        for aa = 1:nSym
            rx_signal(aa,:) = (sqrt(nFFTSize)/nFFTSize)*fftshift(fft(rx_signal(aa,:),nFFTSize));   
            rx_signal(aa,:) = rx_signal(aa,:)./hF(aa,:);
            rx_signal(aa,:) = rx_signal(aa,:) .*olcek(aa,1);
        end
        
        for aa=1:nSym     
            temp_alici_1 = reshape(rx_signal(aa,:),n,g).';
            for a=1:g
                temp_alici_2=temp_alici_1(a,:);
                for j=1:length(crit)
                    olabilir(j) = sum(abs(temp_alici_2-crit(j,:)).^2,2);
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

makale_hata= [0.21 0.1 0.03 0.011 0.0033 0.001 3.6e-4];%makale hatalarý
close all; figure
semilogy(Eb,makale_hata,'--s','LineWidth',2);
hold on
semilogy(Eb,simBer,'--rx','LineWidth',2);
axis([0 30 10^-4 1])
grid on
legend('Makale Sonucu [5]', 'Benzetim');
xlabel('Eb/No, dB')
ylabel('Bit Error Rate')
%title('QPSK Kullanýlmýþ DM-OFDM-IM için BER Grafiði')