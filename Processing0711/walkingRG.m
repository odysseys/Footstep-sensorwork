clear all
close all
clc

load('../result&Images/RushilGabriel.mat');
addpath('../Functions');
Fs = 6500;

for sensorID = 1:4
    sensorID
    [ Y, f, NFFT] = signalFreqencyExtract( noiseSig{sensorID}.signal, Fs );
    figure;
    plot(f,Y);
    figure;
    [c, noiseScale, noiseBandScale] = waveletAnalysis( noiseSig{sensorID}.signal, 1 );
    noiseBandScale
    bandWidth = 50;
    for personID = 1:3
        personID
        for traceID = 1:10
            traceID
            % use first sensor to detect SE
            signal = traceSig{personID, traceID, sensorID};
            
            figure;
            subplot(2,1,1);
            c = waveletAnalysis( signal.signal, 1);
            filteredSig{personID, traceID, sensorID} = waveletFiltering(c, noiseBandScale);
            filteredSig{personID, traceID, sensorID} = filteredSig{personID, traceID, sensorID}.*10;
            filteredIdx{personID, traceID, sensorID} = signal.index;
            subplot(2,1,2);
            plot(signal.index,signal.signal);hold on;
            plot(signal.index,filteredSig{personID, traceID, sensorID});hold off;
        end
    end
end

save('RG_wavelet.mat','filteredSig','filteredIdx');

