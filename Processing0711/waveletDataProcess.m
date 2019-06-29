clear all
close all
clc

% load('../result&Images/AllenGabriel.mat');
% load('AG_wavelet.mat');

% load('../result&Images/RushilGabriel.mat');
% load('RG_wavelet.mat');

load('../result&Images/RushilShijia.mat');
load('RS_wavelet.mat');


personNum = size(traceSig,1);
traceNum = size(traceSig,2);
sensorNum = size(traceSig,3);
tagNum = 2;

addpath('../Functions');

for personID = 1:personNum
    personID
    for traceID = 1:traceNum
        traceID
        windowEnergy = [];
        windowEnergyRelative = [];
        windowAmplitude = [];
        windowShift = [];
        
        timestampStep = 10^6/10;
        for sensorID = 1 : sensorNum
            idx{sensorID} = filteredIdx{personID, traceID, sensorID};
%             idx{sensorID} = idx{sensorID}./2^32; % for RG
            idx{sensorID} = idx{sensorID}-idx{sensorID}(1);
            if sensorID == 1
                timestampRange = idx{1}(end) - idx{1}(1);
            end
            
            sig{sensorID} = filteredSig{personID, traceID, sensorID};
        end
            
        for timestampStart = 1 : timestampStep/2 : timestampRange
            windowStart = timestampStart - 1;
            windowStop = windowStart + timestampStep;
            minLen = -1;
                
            for sensorID = 1:sensorNum
                timestampRange = idx{sensorID}(end) - idx{sensorID}(1);
                
                windowSig{sensorID} = sig{sensorID}(idx{sensorID} >= windowStart & idx{sensorID} < windowStop);
                if sensorID == 1
                    minLen = length(windowSig{sensorID});
                elseif minLen > length(windowSig{sensorID})
                    minLen = length(windowSig{sensorID});
                end
                
            end
            
            if minLen > 0
                for sensorID = 1:sensorNum
                    windowSig{sensorID} = windowSig{sensorID}(1:minLen);
                end
            end
            tempA = []; tempE = []; tempS =[];
            for sensorID = 1:sensorNum
                tempA = [tempA, mean(windowSig{sensorID})];
                tempE = [tempE, sum(windowSig{sensorID}.^2)];
                [v,s] = max(abs(xcorr(windowSig{1},windowSig{sensorID})));
                s = s - length(windowSig{1});
                tempS = [tempS, s];
            end
            windowAmplitude = [windowAmplitude; tempA];
            windowEnergy = [windowEnergy; tempE];
            windowEnergyRelative = [windowEnergyRelative; tempE-tempE(1)];
            windowShift = [windowShift; tempS];
        end
        
        figure;
        subplot(4,2,1);
        for sensorID = 1 : sensorNum
            plot(traceSig{personID, traceID, sensorID}.index, filteredSig{personID, traceID, sensorID}); hold on;
        end
        for tagID = 1 : tagNum
            if ~isempty(groundTruth{personID, traceID, tagID}.index)
                plot(groundTruth{personID, traceID, tagID}.index,groundTruth{personID, traceID, tagID}.signal); hold on;
            end
        end
        hold off;
        axis tight;
        subplot(4,2,3);
        plot(windowEnergy);
        axis tight;
        subplot(4,2,5);
        for sensorID = 1 : sensorNum
            plot(traceSig{personID, traceID, sensorID}.index, traceSig{personID, traceID, sensorID}.signal);hold on;
        end
        hold off;
        axis tight;
%         waveletAnalysis(traceSig{personID, traceID, sensorID}.signal,1);
%         [S,F,T,P] = spectrogram(traceSig{personID, traceID, sensorID}.signal,64,32,1024,6500);
%         imagesc(P);
%         ylim([0,50]);
        
        subplot(4,2,2);
        spectrogram(traceSig{personID, traceID, 1}.signal,64,32,1024,6500,'yaxis');
        subplot(4,2,4);
        spectrogram(traceSig{personID, traceID, 2}.signal,64,32,1024,6500,'yaxis');
        subplot(4,2,6);
        spectrogram(traceSig{personID, traceID, 3}.signal,64,32,1024,6500,'yaxis');
        subplot(4,2,8);
        spectrogram(traceSig{personID, traceID, 4}.signal,64,32,1024,6500,'yaxis');
    end
end