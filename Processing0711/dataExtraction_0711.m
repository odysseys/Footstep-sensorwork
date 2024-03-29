%%%%%%  fixed impact locations from two people %%%%%%
clear all
close all
clc
count = 0;
tagIDOffset = 10;
tagNum = 2;
anchorNum = 4;
addpath('../Functions');

load('RS_stepping.mat');

separatorIdx = find(tagSig == 10000);
status = 0;
tag{1}.signal = [];
tag{2}.signal = [];
tag{1}.index = [];
tag{2}.index = [];
timestamp = -1;
for i = 1 : tagNum
    priorTimestamp{i} = -1;
end
for i = 1 : length(separatorIdx)-1
    if separatorIdx(i+1)-separatorIdx(i) == 2 && status == 0
        tagID = tagSig(separatorIdx(i)+1);
        if tagID == 11 || tagID == 12
            status = 1;
        end
    elseif separatorIdx(i+1)-separatorIdx(i) == 2 && status == 1
        timestamp = tagSig(separatorIdx(i)+1);
%         if timestamp < 0
%             timestamp = timestamp + 2^32;
%         end
        if priorTimestamp{tagID-tagIDOffset} == -1
            priorTimestamp{tagID-tagIDOffset} = timestamp;
        else
            if priorTimestamp{tagID-tagIDOffset} - timestamp > 2^31
                timestamp = timestamp + 2^32;
            end
            priorTimestamp{tagID-tagIDOffset} = timestamp;
        end
        status = 2;
    elseif status == 2 && separatorIdx(i+1)-separatorIdx(i) > 10
        tempSig = tagSig(separatorIdx(i)+1:separatorIdx(i+1)-1)';
        tempSig = tempSig(tempSig>0);
        tag{tagID-tagIDOffset}.signal = [tag{tagID-tagIDOffset}.signal, tempSig];
        tempIdx = zeros(1, length(tempSig));
        tempIdx(end) = timestamp;
        for j = length(tempSig)-1:-1:1
            tempIdx(j) = tempIdx(j+1) - 10^5;
        end
        tag{tagID-tagIDOffset}.index = [tag{tagID-tagIDOffset}.index, tempIdx];
        
        status = 0;
    end 
end

figure;
subplot(2,1,1);
plot(tag{1}.signal(tag{1}.signal<10^9));
subplot(2,1,2);
plot(tag{2}.signal(tag{2}.signal<10^9));
% 
%% assign timestamp to each tag datapoint
for i = 1 : tagNum
    sigMean = mean(tag{i}.signal);
    tag{i}.signal = tag{i}.signal - sigMean;
    sigRange = max(abs(tag{i}.signal));
    tag{i}.signal = tag{i}.signal / sigRange * 250;
end

figure;
for i = 1 : anchorNum
    subplot(anchorNum,1,i);
    plot(s{i}(s{i}<10000));
end
% 
%% assign timestamp to each anchor datapoint
for i = 1 : anchorNum
    timestamps{i} = s{i}(s{i}>10000);
    timestampIdxs{i} = find(s{i}>10000);
    deltaTimestamp = timestamps{i}(2:end)-timestamps{i}(1:end-1);
    timeEva{i} = deltaTimestamp./(timestampIdxs{i}(2:end)-timestampIdxs{i}(1:end-1));
    checkDrop = find(-timeEva{i} > 2^31);
    if ~isempty(checkDrop)
        for tIdx = checkDrop+1: length(timestampIdxs{i})
            s{i}(timestampIdxs{i}(tIdx)) = s{i}(timestampIdxs{i}(tIdx)) + 2^32;
        end
    end
    timestamps{i} = s{i}(s{i}>10000);
    timestampIdxs{i} = find(s{i}>10000);
    
end
% 
%% alignment
% baseline
unifiedTime = round([timestamps{1}; timestamps{2}; timestamps{3}; timestamps{4}]./10^6);
unifiedTime = unique(unifiedTime);
for i = 1 : length(unifiedTime)
    if ismember(unifiedTime(i),round(timestamps{1}./10^6)) ...
        && ismember(unifiedTime(i),round(timestamps{2}./10^6)) ...
        && ismember(unifiedTime(i),round(timestamps{3}./10^6)) ...
        && ismember(unifiedTime(i),round(timestamps{4}./10^6)) 
        baseline = unifiedTime(i);
        break;
    end
end
for i = 1 : 4
    cuttingPoint{i} = timestampIdxs{i}(round(timestamps{i}./10^6) == baseline);
end
% endline
% endline = min(round([timestamps{1}(end),timestamps{2}(end),timestamps{3}(end),timestamps{4}(end)]./10^6));
% for i = 1 : 4
%     endPoint{i} = timestampIdxs{i}(round(timestamps{i}./10^6) == endline);
% end
for i = 1 : 4
    if ~isempty(cuttingPoint{i})
        s{i}(1:cuttingPoint{i}(1)-1) = [];
    end
    timestamps{i} = s{i}(s{i}>10000);
    timestampIdxs{i} = find(s{i}>10000);
%     s{i}(endPoint{1}:end) = [];
end

figure;
subplot(4,1,1);
plot(s{1}(s{1}<10000)); title('after cutting');
axis tight;
subplot(4,1,2);
plot(s{2}(s{2}<10000));
axis tight;
subplot(4,1,3);
plot(s{3}(s{3}<10000));
axis tight;
subplot(4,1,4);
plot(s{4}(s{4}<10000));
axis tight;
% 
% 
%% extract ranging signal
% extract signal between each timestamp
for sensorID = 1 : anchorNum
    sensorID
    % find all the timestamps locations
%     tIdx = (s{sensorID} == 20000);
%     tIdx2 = [0; tIdx(1:end-1)];
    
    dIdx = (s{sensorID} == 10000);
    dIdx2 = [0; dIdx(1:end-1)];

    distance{sensorID} = s{sensorID}.*(dIdx2);
    dst{sensorID} = distance{sensorID}(distance{sensorID}>0);
    distIdx{sensorID} = find(distance{sensorID}>0);
    distTime{sensorID} = zeros(1,length(distIdx{sensorID}));
    for i = 1 : length(distIdx{sensorID})
        for j = 1 : length(timestampIdxs{sensorID})-1
            if distIdx{sensorID}(i) > timestampIdxs{sensorID}(j) && distIdx{sensorID}(i) < timestampIdxs{sensorID}(j+1)
                distTime{sensorID}(i) = timestamps{sensorID}(j)+ ...
                    (timestamps{sensorID}(j+1)-timestamps{sensorID}(j))/...
                    (timestampIdxs{sensorID}(j+1)-timestampIdxs{sensorID}(j))*...
                    (distIdx{sensorID}(i) - timestampIdxs{sensorID}(j));
                break;
            end
        end
        if  distTime{sensorID}(i) == 0
%             fprintf('Meh');
            for j = 1 : length(timestampIdxs{sensorID})
                if distIdx{sensorID}(i) > timestampIdxs{sensorID}(j)
                    distTime{sensorID}(i) = timestamps{sensorID}(j)+...
                        mean(timeEva{sensorID})*(distIdx{sensorID}(i) - timestampIdxs{sensorID}(j));
                    break;
                end
            end
        end
        if  distTime{sensorID}(i) == 0
            distTime{sensorID}(i) = timestamps{sensorID}(1)+...
                mean(timeEva{sensorID})*(distIdx{sensorID}(i) - ...
                timestampIdxs{sensorID}(1));
        end
    end
    for tagID = 1 : tagNum
        d{sensorID,tagID} = [];
    end
    for i = 1 : 2 : length(dst{sensorID})
        if dst{sensorID}(i) == 11 && dst{sensorID}(i+1) < 100
            d{sensorID,1} = [d{sensorID,1}; distIdx{sensorID}(i+1) dst{sensorID}(i+1) distTime{sensorID}(i+1)];
        elseif dst{sensorID}(i) == 12 && dst{sensorID}(i+1) < 100
            d{sensorID,2} = [d{sensorID,2}; distIdx{sensorID}(i+1) dst{sensorID}(i+1) distTime{sensorID}(i+1)];
        end
    end
    
    % remove distance information
    dIdx = find(s{sensorID} == 10000);
    dIdx2 = dIdx+1;
    sigVib{sensorID} = s{sensorID};
    sigVib{sensorID}([dIdx;dIdx2]) = [];
end

figure;
plot(d{1,1}(:,3),d{1,1}(:,2));hold on;
plot(d{2,1}(:,3),d{2,1}(:,2));hold on;
plot(d{3,1}(:,3),d{3,1}(:,2));hold on;
plot(d{4,1}(:,3),d{4,1}(:,2));hold off;
figure;
plot(d{1,2}(:,3),d{1,2}(:,2));hold on;
plot(d{2,2}(:,3),d{2,2}(:,2));hold on;
plot(d{3,2}(:,3),d{3,2}(:,2));hold on;
plot(d{4,2}(:,3),d{4,2}(:,2));hold off;

% 
figure;
subplot(4,1,1);
plot(sigVib{1});
subplot(4,1,2);
plot(sigVib{2});
subplot(4,1,3);
plot(sigVib{3});
subplot(4,1,4);
plot(sigVib{4});

% fine grain alignment
timestamps{1} = sigVib{1}(sigVib{1}>10000);
timestampIdxs{1} = find(sigVib{1}>10000);
timestamps{2} = sigVib{2}(sigVib{2}>10000);
timestampIdxs{2} = find(sigVib{2}>10000);
timestamps{3} = sigVib{3}(sigVib{3}>10000);
timestampIdxs{3} = find(sigVib{3}>10000);
timestamps{4} = sigVib{4}(sigVib{4}>10000);
timestampIdxs{4} = find(sigVib{4}>10000);

for i = 1 : 4
    deltaTime = abs(timestamps{i}(2:end)-timestamps{i}(1:end-1));
    deltaIdx = abs(timestampIdxs{i}(2:end)-timestampIdxs{i}(1:end-1));
    deltaRatio = deltaTime./deltaIdx;
    badIdx = find(deltaRatio > 10^13);
    if ~isempty(badIdx)
        badIdx = badIdx(end);
        sigVib{i}(timestampIdxs{i}(badIdx)) = [];
        timestamps{i} = sigVib{i}(sigVib{i}>10000);
        timestampIdxs{i} = find(sigVib{i}>10000);
    end
end

[sS{1}, sI{1}] = timestampAlignment( sigVib{1}, timestamps{1}, timestampIdxs{1} );
[sS{2}, sI{2}] = timestampAlignment( sigVib{2}, timestamps{2}, timestampIdxs{2} );
[sS{3}, sI{3}] = timestampAlignment( sigVib{3}, timestamps{3}, timestampIdxs{3} );
[sS{4}, sI{4}] = timestampAlignment( sigVib{4}, timestamps{4}, timestampIdxs{4} );
% 
for i = 1 : 4
    sS{i} = sS{i}-mean(sS{i});
end
% 
figure;
subplot(3,1,1);
for i = 1 : 4
    plot(sI{i}(1:10:end), sS{i}(1:10:end));hold on;
end
hold off;
subplot(3,1,2);
for i = 1 : 2
    plot(tag{i}.index,tag{i}.signal);hold on;
end
hold off;
subplot(3,1,3);
plot(d{1,2}(:,3),d{1,2}(:,2));hold on;
plot(d{2,2}(:,3),d{2,2}(:,2));hold on;
plot(d{3,2}(:,3),d{3,2}(:,2));hold on;
plot(d{4,2}(:,3),d{4,2}(:,2));hold off;

figure;
for i = 1 : 4
    plot(sI{i}(1:10:end), sS{i}(1:10:end));hold on;
end
for i = 1 : 2
    plot(tag{i}.index,tag{i}.signal);hold on;
end
xlim([2.158*10^9,2.174*10^9,]);
hold off;

%% specific range
for i = 1 : 4
    tempIdx = find(sI{i} > 2.158*10^9 & sI{i} < 2.174*10^9);
    
    anchorI{i} = sI{i}(tempIdx);
    anchorS{i} = sS{i}(tempIdx);
    noiceS{i} = sS{i}(sI{i} > 2.42*10^9 & sI{i} < 2.48*10^9);
    radioS{i,1} = d{i,1}(d{i,1}(:,3) > 2.158*10^9 & d{i,1}(:,3) < 2.174*10^9 ,2);
    radioI{i,1} = d{i,1}(d{i,1}(:,3) > 2.158*10^9 & d{i,1}(:,3) < 2.174*10^9 ,3);
    radioS{i,2} = d{i,2}(d{i,2}(:,3) > 2.158*10^9 & d{i,2}(:,3) < 2.174*10^9 ,2);
    radioI{i,2} = d{i,2}(d{i,2}(:,3) > 2.158*10^9 & d{i,2}(:,3) < 2.174*10^9 ,3);
end
for i = 1 : 2
    tagS{i} = tag{i}.signal(tag{i}.index > 2.158*10^9 & tag{i}.index < 2.174*10^9);
    tagI{i} = tag{i}.index(tag{i}.index > 2.158*10^9 & tag{i}.index < 2.174*10^9);
end
for i = 1 : 4
    [~, ~, noiseBandScale] = waveletAnalysis( noiceS{i}, 1 );
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

%% step example 
figure;
subplot(3,2,1);
for i = 1% : 4
    indexR1 = find(anchorI{i} > 2.15926 * 10^9 & anchorI{i} < 2.15971 * 10^9);
    stepA1S1 = anchorS{i}(indexR1);
    stepA1I = anchorI{i}(indexR1);
    plot(stepA1I-stepA1I(1), stepA1S1); hold on;
end
% plot([2.15941*10^9, 2.15941*10^9],[-20,20],'r.-');
hold off;
axis tight;
xlabel('Timestamp (us)');
ylabel('Amplitude');
ylim([-17,15]);
title('Sensor 1 Person 1');
% legend('Sensor 1','Sensor 2','Sensor 3','Sensor 4','Person 1 Peak');

subplot(3,2,2);
for i = 4
    indexR1 = find(anchorI{i} > 2.15926 * 10^9 & anchorI{i} < 2.15971 * 10^9);
    stepA1S2 = anchorS{i}(indexR1);
    stepA1I = anchorI{i}(indexR1);
    plot(stepA1I-stepA1I(1), stepA1S2); hold on;
end
% plot([2.15941*10^9, 2.15941*10^9],[-20,20],'r.-');
hold off;
axis tight;
xlabel('Timestamp (us)');
ylabel('Amplitude');
ylim([-17,15]);
title('Sensor 2 Person 1');
% legend('Sensor 1','Sensor 2','Sensor 3','Sensor 4','Person 1 Peak');

subplot(3,2,3);
for i = 1
    indexR2 = find(anchorI{i} > 2.16721 * 10^9 & anchorI{i} < 2.16766 * 10^9);
    stepA2S1 = anchorS{i}(indexR2);
    stepA2I = anchorI{i}(indexR2);
    plot(stepA2I-stepA2I(1), stepA2S1); hold on;
end
% plot([2.16727*10^9, 2.16727*10^9],[-20,20],'b.-');
hold off;
axis tight;
xlabel('Timestamp (us)');
ylabel('Amplitude');
ylim([-17,15]);
title('Sensor 1 Person 2');
% legend('Sensor 1','Sensor 2','Sensor 3','Sensor 4','Person 2 Peak');

subplot(3,2,4);
for i = 4
    indexR2 = find(anchorI{i} > 2.16721 * 10^9 & anchorI{i} < 2.16766 * 10^9);
    stepA2S2 = anchorS{i}(indexR2);
    stepA2I = anchorI{i}(indexR2);
    plot(stepA2I-stepA2I(1), stepA2S2); hold on;
end
% plot([2.16727*10^9, 2.16727*10^9],[-20,20],'b.-');
hold off;
axis tight;
xlabel('Timestamp (us)');
ylabel('Amplitude');
ylim([-17,15]);
title('Sensor 2 Person 2');
% legend('Sensor 1','Sensor 2','Sensor 3','Sensor 4','Person 2 Peak');

subplot(3,2,5);
for i = 1
    indexR3 = find(anchorI{i} > 2.16834 * 10^9 & anchorI{i} < 2.16879 * 10^9);
    stepA3S1 = anchorS{i}(indexR3);
    stepA3I = anchorI{i}(indexR3);
    plot(stepA3I(1:length(stepA1S1))-stepA3I(1), stepA3S1(1:length(stepA1S1))); hold on;
end
together1 = stepA1S1 + stepA2S1(1:length(stepA1S1));
plot(stepA3I(1:length(stepA1S1))-stepA3I(1),together1,'r-.');
% plot([2.16849*10^9, 2.16849*10^9],[-20,20],'r.-');
% plot([2.1684*10^9, 2.1684*10^9],[-20,20],'b.-');
hold off;
axis tight;
xlabel('Timestamp (us)');
ylabel('Amplitude');
ylim([-17,15]);
title('Sensor 1 Person 1&2');
% legend('Sensor 1','Sensor 2','Sensor 3','Sensor 4','Person 1 Peak','Person 2 Peak');

subplot(3,2,6);
for i = 4
    indexR3 = find(anchorI{i} > 2.16834 * 10^9 & anchorI{i} < 2.16879 * 10^9);
    stepA3S2 = anchorS{i}(indexR3);
    stepA3I = anchorI{i}(indexR3);
    plot(stepA3I(1:length(stepA1S2))-stepA3I(1), stepA3S2(1:length(stepA1S2))); hold on;
end
together2 = stepA1S2 + stepA2S2(1:length(stepA1S2));
error2 = (abs(together2 - stepA3S2(1:length(stepA1S2)))./abs(stepA3S2(1:length(stepA1S2))));
plot(stepA3I(1:length(stepA1S2))-stepA3I(1), together2,'r-.');
% plot([2.16849*10^9, 2.16849*10^9],[-20,20],'r.-');
% plot([2.1684*10^9, 2.1684*10^9],[-20,20],'b.-');
hold off;
axis tight;
xlabel('Timestamp (us)');
ylabel('Amplitude');
ylim([-17,15]);
title('Sensor 2 Person 1&2');
% legend('Sensor 1','Sensor 2','Sensor 3','Sensor 4','Person 1 Peak','Person 2 Peak');

% subplot(4,2,7);
% for i = 1
%     indexR4 = find(anchorI{i} > 2.16905 * 10^9 & anchorI{i} < 2.16965 * 10^9);
%     stepA4S = anchorS{i}(indexR4);
%     stepA4I = anchorI{i}(indexR4);
%     plot(stepA4I, stepA4S); hold on;
% end
% plot([2.16913*10^9, 2.16913*10^9],[-20,20],'r.-');
% plot([2.16922*10^9, 2.16922*10^9],[-20,20],'b.-');
% hold off;
% axis tight;
% xlabel('Timestamp');
% % ylabel('Amplitude');legend('Sensor 1','Sensor 2','Sensor 3','Sensor 4','Person 1 Peak','Person 2 Peak');
% 
% subplot(4,2,8);
% for i = 4
%     indexR4 = find(anchorI{i} > 2.16905 * 10^9 & anchorI{i} < 2.16965 * 10^9);
%     stepA4S = anchorS{i}(indexR4);
%     stepA4I = anchorI{i}(indexR4);
%     plot(stepA4I, stepA4S); hold on;
% end
% plot([2.16913*10^9, 2.16913*10^9],[-20,20],'r.-');
% plot([2.16922*10^9, 2.16922*10^9],[-20,20],'b.-');
% hold off;
% axis tight;
% xlabel('Timestamp');
% % ylabel('Amplitude');legend('Sensor 1','Sensor 2','Sensor 3','Sensor 4','Person 1 Peak','Person 2 Peak');

    
