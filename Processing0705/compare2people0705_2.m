% test out identification

clear all
close all
clc

load('RushilGabriel.mat');
Fs = 15000;
cutoffFrequency = 200;
sensorID = 1;
idAccuracyRate = [];
addpath('./libsvm-master/matlab/');
addpath('../Functions');

for filterRange = 1:25:200
    filterRange
    [ Y, f, NFFT] = signalFreqencyExtract( noiseSig{sensorID}.signal, Fs );
    figure;
    plot(f,Y);
    figure;
    bandWidth = 40;
    [c, noiseScale, noiseBandScale] = waveletAnalysis( noiseSig{sensorID}.signal, bandWidth);
    
    %% define the results format
    trainingFeatures = [];
    trainingLabels = [];
    testingFeatures = [];
    testingLabels = [];
    
    %%
    for personID = [1,3]
        personID
        % training
        for traceID = 1:5
            traceID
            % use first sensor to detect SE
            signal = traceSig{personID, traceID, sensorID};
            signal.signal = signalDenoise(signal.signal,50);
            noiseSig{sensorID}.signal = signalDenoise(noiseSig{sensorID}.signal,bandWidth);
            % wavelet filtering
            c = waveletAnalysis( signal.signal, waveletAnalysis );
            filteredSig{personID, traceID, sensorID} = waveletFiltering(c, filterRange:filterRange+40);
            
            [ stepEventsSig, stepEventsIdx, stepEventsVal, ...
            stepStartIdxArray, stepStopIdxArray, ... 
            windowEnergyArray, noiseMu, noiseSigma, ...
            noiseRange ] = SEDetection( filteredSig{personID, traceID, sensorID}, noiseSig{sensorID}.signal,16,Fs );
            
            [~, maxSigIdx] = sort(abs(stepEventsVal), 'descend');
            maxSigIdx = maxSigIdx(1:5);
            maxSigIdx = sort(maxSigIdx,'ascend');
            stepEventsIdx = stepEventsIdx(maxSigIdx);
            stepEventsVal = stepEventsVal(maxSigIdx);
            stepStartIdxArray = stepStartIdxArray(maxSigIdx);
            stepStopIdxArray = stepStopIdxArray(maxSigIdx);
            stepNum = length(stepEventsIdx);
            for stepIdx = 1 : stepNum
                stepSig = filteredSig{personID, traceID, sensorID}(stepStartIdxArray(stepIdx):stepStopIdxArray(stepIdx));
                [ Y, f, NFFT] = signalFreqencyExtract( stepSig, Fs );
                Y = Y(f<=cutoffFrequency);
                f = f(f<=cutoffFrequency);
                Y = signalNormalization(Y);
                trainingFeatures = [trainingFeatures; Y];
                trainingLabels = [trainingLabels; personID];
            end
        end
    end
    
    accBase = 0;
    gc = 0;
    for gi = [0.1,0.5,1,3,5,7]
        tempStruct = svmtrain(trainingLabels, trainingFeatures, ['-s 0 -t 2 -b 1 -g ' num2str(gi) ' -c 100']);
        [predicted_label, accuracy, decision_values] = svmpredict(trainingLabels, trainingFeatures, tempStruct,'-b 1');
        if accuracy(1) > accBase
            svmstruct = tempStruct;
            accBase = accuracy(1);
            gc = gi;
        end
    end
        
    for personID = [1,3]
        personID
        for traceID = 6:10
            traceID
            % use first sensor to detect SE
            signal = traceSig{personID, traceID, sensorID};
            signal.signal = signalDenoise(signal.signal,50);
            noiseSig{sensorID}.signal = signalDenoise(noiseSig{sensorID}.signal,50);
            % wavelet filtering
            c = waveletAnalysis( signal.signal, waveletAnalysis );
            filteredSig{personID, traceID, sensorID} = waveletFiltering(c, filterRange:filterRange+40);
            
            [ stepEventsSig, stepEventsIdx, stepEventsVal, ...
            stepStartIdxArray, stepStopIdxArray, ... 
            windowEnergyArray, noiseMu, noiseSigma, ...
            noiseRange ] = SEDetection( signal.signal, noiseSig{sensorID}.signal,16,Fs );
            
            [~, maxSigIdx] = sort(abs(stepEventsVal), 'descend');
            maxSigIdx = maxSigIdx(1:5);
            maxSigIdx = sort(maxSigIdx,'ascend');
            stepEventsIdx = stepEventsIdx(maxSigIdx);
            stepEventsVal = stepEventsVal(maxSigIdx);
            stepStartIdxArray = stepStartIdxArray(maxSigIdx);
            stepStopIdxArray = stepStopIdxArray(maxSigIdx);
            
            stepNum = length(stepEventsIdx);
            for stepIdx = 1 : stepNum
                stepSig = filteredSig{personID, traceID, sensorID}(stepStartIdxArray(stepIdx):stepStopIdxArray(stepIdx));
                [ Y, f, NFFT] = signalFreqencyExtract( stepSig, Fs );
                Y = Y(f<=cutoffFrequency);
                f = f(f<=cutoffFrequency);
                Y = signalNormalization(Y);
                testingFeatures = [testingFeatures; Y];
                testingLabels = [testingLabels; personID];
            end
            
        end
    end
    [predicted_label, accuracy, decision_values] = svmpredict(testingLabels, testingFeatures, svmstruct,'-b 1');
    idAccuracyRate = [idAccuracyRate; filterRange, accuracy(1)];
end
figure;
plot(idAccuracyRate(:,1),idAccuracyRate(:,2));
save('RushilGabriel_idresutls.mat','idAccuracyRate');
