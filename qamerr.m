function [ error ] = qamerr( EbNo,dataMod,dataIn,char)
  
k=4;
numSamplesPerSymbol=1;

snr = EbNo + 10*log10(k) - 10*log10(numSamplesPerSymbol);
receivedSignal = awgn(dataMod,snr,'measured');

% QAM Demodulation

if strcmp(char,'Without_filter')
   
    % (A) Without filtering
    
    dataSymbolsOut = qamdemod(receivedSignal,16,0,'bin');
    dataOutMatrix = de2bi(dataSymbolsOut,k);
    dataOut = dataOutMatrix(:);
    
    [error,ber] = biterr(dataIn,dataOut);
    

    
else
    mag=abs(receivedSignal);
phase=angle(receivedSignal);

sm1=mag.*cos(phase);
sm2=mag.*sin(phase);

% the orthogonal basis
fs=100;
t=0:(1/fs):1;
cosine=(1/sqrt(51))*cos(2*pi*t);
sine=(1/sqrt(50))*sin(2*pi*t);

% QAM signal
result=sm1*cosine +sm2*sine;

% Signal without filter
result_nf=result;

% Apply median filtering here using built-in function
i=1;

while i<=7500
    result(i,:)=medfilt1(result(i,:),7);
    i=i+1;
end

% to get the coefficients of the sine and cos terms
cos_coeff=result*cosine';
sin_coeff=result*sine';

% to get the phase angle
phase_result=atan(sin_coeff./cos_coeff);

% as the phase crosses \pi, we need to limit the value to obtain same
% values
i=1;
while i<=length(phase_result)
    
        if phase(i) < -pi/2
            phase_result(i) = phase_result(i) - pi;
        end
        if phase(i) > pi/2
            phase_result(i) = pi + phase_result(i);
        end
        i=i+1;
            
end

% calculating amplitude
mag_result=(cos_coeff.^2 + sin_coeff.^2).^0.5;

% final filtered QAM signal
filtered_result=mag_result.*exp(j*phase_result);


% Demodulation of filtered signal
filtered_dataoutput=qamdemod(filtered_result,16,0,'bin');
filtered_binaryresult=de2bi(filtered_dataoutput,k);
filtered_out=filtered_binaryresult(:);
error = mdnFilter(dataIn,receivedSignal, 7);

end


    




end

