clear all
close all
clc

load('./result&Images/RushilGabriel.mat');
Fs = 15000;

for sensorID = 1%:4
    sensorID
    [ Y, f, NFFT] = signalFreqencyExtract( noiseSig{sensorID}.signal, Fs );
    figure;
    plot(f,Y);
    figure;
    [c, noiseScale, noiseBandScale] = waveletAnalysis( noiseSig{sensorID}.signal );
    noiseBandScale
    bandWidth = 50;
    for personID = 1:3
        personID
        for traceID = 1:10
            traceID
            % use first sensor to detect SE
            signal = traceSig{personID, traceID, sensorID};
            signal.signal = signalDenoise(signal.signal,50);
            noiseSig{sensorID}.signal = signalDenoise(noiseSig{sensorID}.signal,50);
            figure;
            subplot(2,1,1);
            c = waveletAnalysis( signal.signal );
            filteredSig{personID, traceID, sensorID} = waveletFiltering(c, 50:90);

            subplot(2,1,2);
            plot(signal.index,signal.signal);hold on;
            plot(signal.index,filteredSig{personID, traceID, sensorID});hold off;
        end
    end
end

% save('RushilGabriel_filtered_50_90.mat','filteredSig','noiseSig','traceSig');



