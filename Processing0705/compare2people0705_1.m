% compare localization

clear all
close all
clc

load('RushilGabriel_filtered_50_90.mat');
Fs = 6500;
SensorLocations = [0,0;0,3;5,0;3.5,0];
options = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','none');

sensorID = 1;
sensorNum = 4;
for personID = [1,3]
    personID
    for traceID = 1:10
        traceID
        
        signal.index = traceSig{personID, traceID, sensorID}.index;
        signal.signal = signalDenoise(filteredSig{personID, traceID, sensorID},50);
        noiseSig{sensorID}.signal = signalDenoise(noiseSig{sensorID}.signal,50);           

        % extract steps and localize
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
        
        figure;
        plot(signal.index,signal.signal);hold on;
        scatter(signal.index(stepEventsIdx), signal.signal(stepEventsIdx));hold on;
        stepNum = length(stepEventsIdx);
        for stepIdx = 1 : stepNum
            tempIdx = stepStartIdxArray(stepIdx);
            plot([signal.index(tempIdx),signal.index(tempIdx)],[-550, 500],'r');
            tempIdx = stepStopIdxArray(stepIdx);
            plot([signal.index(tempIdx),signal.index(tempIdx)],[-550, 500],'g');
        end
        hold off;
        
        % extract steps within one trace from 4 sensors
        figure;
        peakTimestamp = zeros(stepNum,sensorNum);
        firstPeakTDoA = zeros(stepNum,sensorNum);
        xcorrTDoA = zeros(stepNum,sensorNum);
        calcLoc = zeros(stepNum,4); % 1,2 for first peak, 3,4 for xcorr
            
        for stepIdx = 1 : stepNum
            subplot(1,stepNum,stepIdx);
            startTimestamp = signal.index(stepStartIdxArray(stepIdx));
            stopTimestamp = signal.index(stepStopIdxArray(stepIdx));
            
            for sensorIdx = 1 : sensorNum
                tempSigIdx = traceSig{personID, traceID, sensorIdx}.index;
                tempSigSig = filteredSig{personID, traceID, sensorIdx};
                tempIdx = find(tempSigIdx > startTimestamp & ...
                    tempSigIdx < stopTimestamp);
                tempSig = tempSigSig(tempIdx); 
                plot(tempIdx, tempSig); hold on;
                [peaks, locs,w,p] = findpeaks(tempSig,'MinPeakDistance',Fs/500,'MinPeakHeight',max(tempSig)/2,'Annotate','extents');
                scatter(tempIdx(locs(1)),peaks(1),'rv');
                peakTimestamp(stepIdx, sensorIdx) = tempIdx(locs(1));
                stepSigSet{personID, traceID, stepIdx, sensorIdx} = [tempIdx';tempSig];
                
                if sensorIdx == 1
                    referenceStep = stepSignalNormalization( tempSig );
                else
                    compareStep = stepSignalNormalization( tempSig );
                    [v,sh] = max(abs(xcorr(referenceStep, compareStep)));
                    xcorrTDoA(stepIdx, sensorIdx) = sh -length(referenceStep);
                end
            end
            hold off;
            
            %% obtaining TDoA
            % use first peak
            firstPeakTDoA(stepIdx,:) = peakTimestamp(stepIdx,:) - peakTimestamp(stepIdx,1);
            
            for v = [0.05]
                x0 = [0;0];
                pc = SensorLocations(1,:);
                pi = SensorLocations(2,:);
                pj = SensorLocations(3,:);
                v12 = v;
                v13 = v*1.5;

                tic = -firstPeakTDoA(stepIdx,2);
                tjc = -firstPeakTDoA(stepIdx,3);
                [x_firstPeak, fval, exitflag] = fsolve(@(x) localizationEquations( x, pi, pj, pc, tic, tjc, v12, v13 ),x0,options);
                if validateLocation( SensorLocations, x_firstPeak ) == 1
                    calcLoc(stepIdx,1) = x_firstPeak(1);
                    calcLoc(stepIdx,2) = x_firstPeak(2);
                end
                
                
                tic = xcorrTDoA(stepIdx,2);
                tjc = xcorrTDoA(stepIdx,3);
                [x_xcorr, fval, exitflag] = fsolve(@(x) localizationEquations( x, pi, pj, pc, tic, tjc, v12, v13 ),x0,options);
                if validateLocation( SensorLocations, x_xcorr ) == 1
                    calcLoc(stepIdx,3) = x_xcorr(1);
                    calcLoc(stepIdx,4) = x_xcorr(2);
                end
                
            end
        end
        figure;
        subplot(2,1,1);
        plot(firstPeakTDoA);
        subplot(2,1,2);
        plot(xcorrTDoA);
        
        %% trilateration with these two
        % first peak
        figure;
        subplot(2,1,1);
        scatter(SensorLocations(:,1),SensorLocations(:,2),'ro');hold on;
        scatter(calcLoc(:,1),calcLoc(:,2),'b*'); hold off;
        subplot(2,1,2);
        scatter(SensorLocations(:,1),SensorLocations(:,2),'ro');hold on;
        scatter(calcLoc(:,3),calcLoc(:,4),'b*'); hold off;
    end
end

