Fs = 6500;
TR{1} = [2.158*10^9, 2.174*10^9];
TR{2} = [2.183*10^9, 2.198*10^9];
TR{3} = [2.202*10^9, 2.216*10^9];

SensorLocations = [0,0;0,3;5,0;5,3];

for trIdx = 1%:3
    traceRange = TR{trIdx};
    %% specific range
    for i = 1 : 4
        tempIdx = find(sI{i} > traceRange(1) & sI{i} < traceRange(2));

        anchorI{i} = sI{i}(tempIdx);
        anchorS{i} = sS{i}(tempIdx);
        noiceS{i} = sS{i}(sI{i} > 2.42*10^9 & sI{i} < 2.48*10^9);
        radioS{i,1} = d{i,1}(d{i,1}(:,3) > traceRange(1) & d{i,1}(:,3) < traceRange(2) ,2);
        radioI{i,1} = d{i,1}(d{i,1}(:,3) > traceRange(1) & d{i,1}(:,3) < traceRange(2) ,3);
        radioS{i,2} = d{i,2}(d{i,2}(:,3) > traceRange(1) & d{i,2}(:,3) < traceRange(2) ,2);
        radioI{i,2} = d{i,2}(d{i,2}(:,3) > traceRange(1) & d{i,2}(:,3) < traceRange(2) ,3);
    end
    for i = 1 : 2
        tagS{i} = tag{i}.signal(tag{i}.index > traceRange(1) & tag{i}.index < traceRange(2));
        tagI{i} = tag{i}.index(tag{i}.index > traceRange(1) & tag{i}.index < traceRange(2));
    end
    for i = 1 : 4
        [cN, ~, noiseBandScale] = waveletAnalysis( noiceS{i}, 1 );
        noiseFiltered{i} = waveletFiltering(cN, noiseBandScale);
        c = waveletAnalysis( anchorS{i}, noiseBandScale );
        anchorS{i} = waveletFiltering(c, noiseBandScale);
    end

    figure;
    subplot(3,1,1);
    for i = 1 : 4
        plot(anchorI{i}, anchorS{i}.*25); hold on;
    end
    for i = 1 : 2
        plot(tagI{i},tagS{i}); hold on;
    end
    subplot(3,1,2);
    for i = 1 : 4
        plot(radioI{i,1}, radioS{i,1} ); hold on;
    end
    subplot(3,1,3);
    for i = 1 : 4
        plot(radioI{i,2}, radioS{i,2} ); hold on;
    end    

    %% calculate location of person 1 and person 2 based on radio signal
    for personID = 1:2
        for sensorID = [1,2,4]
            ToF{personID, sensorID} = mean(radioS{sensorID,personID}) -  1.2;
        end
    end

    windowEnergy = [];
    windowEnergyRelative = [];
    windowAmplitude = [];
    windowShift = [];
    windowIdx = [];

    timestampRange = anchorI{1}(end) - anchorI{1}(1);
    timestampStep = 10^6/40;

    for timestampStart = 1 : timestampStep/2 : timestampRange
        windowStart = anchorI{sensorID}(1)+timestampStart - 1;
        windowStop = windowStart + timestampStep;
        minLen = -1;
        for sensorID = 1:4
            windowSig{sensorID} = anchorS{sensorID}(anchorI{sensorID} >= windowStart & anchorI{sensorID} < windowStop);
            if sensorID == 1
                minLen = length(windowSig{sensorID});
            elseif minLen > length(windowSig{sensorID})
                minLen = length(windowSig{sensorID});
            end
        end
        if minLen > 0
            for sensorID = 1:4
                windowSig{sensorID} = windowSig{sensorID}(1:minLen);
            end
        end
        tempA = []; tempE = []; tempS =[];
        for sensorID = 1:4
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
        windowIdx = [windowIdx; windowStart];
    end

    figure;
    subplot(2,1,1);
    plot(windowEnergy);
    subplot(2,1,2);
    for sensorID = 1:4
        plot(anchorI{sensorID}, anchorS{sensorID});hold on;
        plot(windowIdx, windowEnergy./1000);
    end


    %% step extraction based on Gaussian noise model
    [ stepEventsSig, stepEventsIdx, stepEventsVal, ...
                stepStartIdxArray, stepStopIdxArray, ... 
                windowEnergyArray, noiseMu, noiseSigma, noiseRange ] = SEDetection( anchorS{1}', noiseFiltered{1} );
    figure;
    subplot(2,1,1);
    plot(anchorI{1},anchorS{1});hold on;
    for stepIdx = 1:length(stepStartIdxArray)
        plot([anchorI{1}(stepStartIdxArray(stepIdx)),anchorI{1}(stepStartIdxArray(stepIdx))],[-15,15],'r');
        plot([anchorI{1}(stepStopIdxArray(stepIdx)),anchorI{1}(stepStopIdxArray(stepIdx))],[-15,15],'g');
    end
    axis tight;
    subplot(2,1,2);
    plot(anchorI{1}(windowEnergyArray(:,2)),windowEnergyArray(:,1));
    axis tight;   

    %% step extraction
    figure;
    for sensorID = 1:4
        plot(windowEnergy(:,sensorID));hold on;
        sigMean = mean(windowEnergy(:,sensorID));
        sigStd = std(windowEnergy(:,sensorID));

        [ val ,idx ] = findpeaks(windowEnergy(:,sensorID),'MinPeakDistance',Fs/40/40,...
                                                'MinPeakHeight',sigMean+sigStd,'Annotate','extents');
        scatter(idx,val,'rv');
        %% find energy peaks within the detected step range
        for stepIdx = 1:length(stepStartIdxArray)
             peaksInWindow{stepIdx, sensorID} = windowIdx(windowIdx(idx)>anchorI{1}(stepStartIdxArray(stepIdx)) &...
                                windowIdx(idx)<anchorI{1}(stepStopIdxArray(stepIdx)));
        end

    end 
    hold off;

    %% localization
    for stepIdx = 1:length(stepStartIdxArray)
        onestep = 1;
        for sensorID = 1:4
            if length(peaksInWindow{stepIdx, sensorID}) > 1
                onestep = 0;
            end
        end
        if onestep == 1
            % find the peak shift
        else
            
        end
    end


   

    
end




%% gccphat
% gccphatResult = [];
% for sensorID = 2:4
%     [v,i] = gccphat(stSig{1}(1:3227), stSig{sensorID}(1:3227));
%     gccphatResult = [gccphatResult; v,i];
% end

%% get steps and normalize windowed segment
% [ val ,idx ] = findpeaks(windowEnergy(:,1),'MinPeakDistance',Fs/40/10,...
%                                             'MinPeakHeight',sigMean,'Annotate','extents');
% figure; plot(anchorI{1},anchorS{1});hold on;
% plot(anchorI{2},anchorS{2});
% plot(anchorI{3},anchorS{3});
% plot(anchorI{4},anchorS{4});
% 
% scatter(idx*timestampStep/2, ones(length(idx),1)*15,'rv');
% for stepIdx  = 1 : length(idx)
%     plot([idx*timestampStep/2-150000,idx*timestampStep/2-150000],[-15,15],'r');
% end
% for stepIdx  = 1 : length(idx)
%     plot([idx*timestampStep/2+250000,idx*timestampStep/2+250000],[-15,15],'g');
% end
% hold off;
% 
% for stepIdx = 1 : length(idx)
%     windowStart = idx(stepIdx)*timestampStep/2-150000;
%     windowStop = idx(stepIdx)*timestampStep/2+250000;
%     sig1 = anchorS{1}(anchorI{1} >= windowStart & anchorI{1} < windowStop);
%     sig4 = anchorS{4}(anchorI{4} >= windowStart & anchorI{4} < windowStop);
%     
%     % segment each step signal into 2,4,8 parts
%     windowInSig1 = floor(length(sig1)/4);
%     windowInSig2 = floor(length(sig4)/4);
%     shifts = [];
%     for windowIdx = 1 : 8
%         range1 = (windowIdx-1)*floor(windowInSig1/4)+1:(windowIdx-1)*floor(windowInSig1/4)+windowInSig1;
%         range2 = (windowIdx-1)*floor(windowInSig2/4)+1:(windowIdx-1)*floor(windowInSig2/4)+windowInSig2;
%         windowXcorr = normXcorr(sig1(range1),sig4(range2),1);
%         [~,shift1] = max(windowXcorr);
%         shift1 = shift1-length(range1);
%         shifts = [shifts, shift1];
%     end
% end

%% dtw
% for stepIdx = 1 : length(idx)
%     windowStart = idx(stepIdx)*timestampStep/2-150000;
%     windowStop = idx(stepIdx)*timestampStep/2+250000;
%     sig1 = anchorS{1}(anchorI{1} >= windowStart & anchorI{1} < windowStop);
%     sig4 = anchorS{4}(anchorI{4} >= windowStart & anchorI{4} < windowStop);
%     
%     [Dist,D,k,w] = dtw(sig1, sig4);
%     dtwResult{stepIdx}.Dist = Dist;
%     dtwResult{stepIdx}.D = D;
%     dtwResult{stepIdx}.k = k;
%     dtwResult{stepIdx}.w = w;
%     figure;imagesc(D);
% end
% save('dtwresults.mat','dtwResult');
% load ('dtwresults.mat');
% distAll = [];
% for stepIdx = 1 : length(idx)
%     windowStart = idx(stepIdx)*timestampStep/2-150000;
%     windowStop = idx(stepIdx)*timestampStep/2+250000;
%     sig1 = anchorS{1}(anchorI{1} >= windowStart & anchorI{1} < windowStop);
%     sig4 = anchorS{4}(anchorI{4} >= windowStart & anchorI{4} < windowStop);
%     figure;
%     subplot(2,1,1); plot(sig1);hold on; plot(sig4);hold off; 
%     subplot(2,1,2); imagesc(dtwResult{stepIdx}.D);
%     distAll = [distAll, dtwResult{stepIdx}.Dist];
% end
% figure; plot(distAll);

