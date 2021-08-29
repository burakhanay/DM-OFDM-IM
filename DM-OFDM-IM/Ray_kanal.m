function [ht,hF] = Ray_kanal(nFFTSize,nSym,nTap)
	ht = sqrt(1/2)*sqrt(1/nTap)*(randn(nSym,nTap)+1j*randn(nSym,nTap));
    hF = zeros(nSym,nFFTSize);
    for aa = 1:nSym
        hF(aa,:) = fftshift(fft(ht(aa,:),nFFTSize,2));
    end
end