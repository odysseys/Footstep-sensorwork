clear all
close all
clc

load('../result&Images/AllenGabriel_raw.mat');
addpath('../Functions');
Fs = 6500;

figure;
subplot(2,1,1);
for i = 1 : 4
    plot(sI{i}(1:10:end), sS{i}(1:10:end));hold on;
end
hold off;
subplot(2,1,2);
for i = 1 : 2
    plot(tag{i}.index,tag{i}.signal);hold on;
end
hold off;

figure;
for i = 1 : 4
    plot(sI{i}(1:10:end), sS{i}(1:10:end));hold on;
end
for i = 1 : 2
    plot(tag{i}.index,tag{i}.signal);hold on;
end
xlim([2.5*10^9,3.15*10^9,]);
hold off;

%% for each case we extract 6 traces
personTraceTimestamp{1,1} = [2.512*10^9, 2.522*10^9];
personTraceTimestamp{1,2} = [2.524*10^9, 2.533*10^9];
personTraceTimestamp{1,3} = [2.538*10^9, 2.547*10^9];
personTraceTimestamp{1,4} = [2.5485*10^9, 2.558*10^9];
personTraceTimestamp{1,5} = [2.56*10^9, 2.57*10^9];
personTraceTimestamp{1,6} = [2.572*10^9, 2.582*10^9];

personTraceTimestamp{2,1} = [2.698*10^9, 2.708*10^9];
personTraceTimestamp{2,2} = [2.723*10^9, 2.733*10^9];
personTraceTimestamp{2,3} = [2.744*10^9, 2.754*10^9];
personTraceTimestamp{2,4} = [2.766*10^9, 2.775*10^9];
personTraceTimestamp{2,5} = [2.788*10^9, 2.7985*10^9];
personTraceTimestamp{2,6} = [2.815*10^9, 2.825*10^9];

personTraceTimestamp{3,1} = [2.886*10^9, 2.896*10^9];
personTraceTimestamp{3,2} = [2.897*10^9, 2.907*10^9];
personTraceTimestamp{3,3} = [2.936*10^9, 2.946*10^9];
personTraceTimestamp{3,4} = [3.002*10^9, 3.012*10^9];
personTraceTimestamp{3,5} = [3.015*10^9, 3.03*10^9];
personTraceTimestamp{3,6} = [3.03*10^9, 3.044*10^9];

noiseTimestamp = [3.098,3.108].*10^9;

for personID = 1:3
    for traceID = 1:6
        for sensorID = 1:4
            traceSig{personID, traceID, sensorID}.signal = sS{sensorID}( ...
                            sI{sensorID}>personTraceTimestamp{personID,traceID}(1) & ...
                            sI{sensorID}<personTraceTimestamp{personID,traceID}(2));
            traceSig{personID, traceID, sensorID}.index = sI{sensorID}( ...
                            sI{sensorID}>personTraceTimestamp{personID,traceID}(1) & ...
                            sI{sensorID}<personTraceTimestamp{personID,traceID}(2));
        end
        for tagID = 1:2
            groundTruth{personID, traceID, tagID}.signal = tag{tagID}.signal( ...
                            tag{tagID}.index>personTraceTimestamp{personID,traceID}(1) & ...
                            tag{tagID}.index<personTraceTimestamp{personID,traceID}(2));
            groundTruth{personID, traceID, tagID}.index = tag{tagID}.index( ...
                            tag{tagID}.index>personTraceTimestamp{personID,traceID}(1) & ...
                            tag{tagID}.index<personTraceTimestamp{personID,traceID}(2));
        end
    end
end

for sensorID = 1:4
    noiseSig{sensorID}.index = sI{sensorID}( ...
                            sI{sensorID}>noiseTimestamp(1) & ...
                            sI{sensorID}<noiseTimestamp(2));
    noiseSig{sensorID}.signal = sS{sensorID}( ...
                            sI{sensorID}>noiseTimestamp(1) & ...
                            sI{sensorID}<noiseTimestamp(2));
end
     
save('../result&Images/AllenGabriel.mat','traceSig','noiseSig','groundTruth');

% for sensorID = 1:4
%     sensorID
%     [ Y, f, NFFT] = signalFreqencyExtract( noiseSig{sensorID}.signal, Fs );
%     figure;
%     plot(f,Y);
%     figure;
%     [c, noiseScale, noiseBandScale] = waveletAnalysis( noiseSig{sensorID}.signal, 1 );
%     noiseBandScale
%     bandWidth = 50;
%     for personID = 1:3
%         personID
%         for traceID = 1:10
%             traceID
%             % use first sensor to detect SE
%             signal = traceSig{personID, traceID, sensorID};
%             
%             figure;
%             subplot(2,1,1);
%             c = waveletAnalysis( signal.signal, 1);
%             filteredSig{personID, traceID, sensorID} = waveletFiltering(c, noiseBandScale);
%             filteredSig{personID, traceID, sensorID} = filteredSig{personID, traceID, sensorID}.*10;
%             filteredIdx{personID, traceID, sensorID} = signal.index;
%             subplot(2,1,2);
%             plot(signal.index,signal.signal);hold on;
%             plot(signal.index,filteredSig{personID, traceID, sensorID});hold off;
%         end
%     end
% end
% 
% save('AG_wavelet.mat','filteredSig');

