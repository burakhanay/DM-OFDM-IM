function noise = addNoise(sig,reqSNR)
    sigPower = sum(abs(sig(:)).^2)/length(sig(:));
    %sigPower = 10^(sigPower/10);
    reqSNR = 10^(reqSNR/10);
    noisePower = sigPower/reqSNR;
    noise = sqrt(noisePower/2)*(randn(stream, size(sig))+1i*randn(stream, size(sig)));
end