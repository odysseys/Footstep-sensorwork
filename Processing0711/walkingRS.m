clear all
close all
clc

load('../result&Images/RushilShijia.mat');
addpath('../Functions');
Fs = 6500;


for sensorID = 1:4
    sensorID
    [ Y, f, NFFT] = signalFreqencyExtract( noiseSig{sensorID}.signal, Fs );
%     figure;
%     plot(f,Y);
%     figure;
    [c, noiseScale, noiseBandScale] = waveletAnalysis( noiseSig{sensorID}.signal, 1 );
    sensorScale{sensorID} = noiseBandScale;
    bandWidth = 50;
    for personID = 1:3
        personID
        for traceID = 1:10
            traceID
            % use first sensor to detect SE
            signal = traceSig{personID, traceID, sensorID};
            c = waveletAnalysis( signal.signal, 1);
            
%             figure;
%             subplot(2,1,1);
%             imagesc(c.cfs);hold on;
%             plot([1,size(c.cfs,2)],[sensorScale{sensorID},sensorScale{sensorID}],'r');hold off;
% 
            filteredSig{personID, traceID, sensorID} = waveletFiltering(c, noiseBandScale);
            filteredSig{personID, traceID, sensorID} = filteredSig{personID, traceID, sensorID}.*10;
            filteredIdx{personID, traceID, sensorID} = signal.index;
%             subplot(2,1,2);
%             plot(signal.index,signal.signal);hold on;
%             plot(signal.index,filteredSig{personID, traceID, sensorID});hold off;
        end
    end
end

save('RS_wavelet.mat','filteredSig','filteredIdx','sensorScale');

