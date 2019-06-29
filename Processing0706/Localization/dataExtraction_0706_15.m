%  import to check
clear all
close all
clc

listing = dir('../Anchor1/copyFromPi15/Footstep/');
count = 0;
tagIDOffset = 10;

for i = 1 : length(listing)
    if sum(ismember(listing(i).name,'txtout')) == 9
        count = count + 1;
        signal1{count} = importdata(['../Anchor1/copyFromPi15/Footstep/' listing(i).name]);
    end
end

listing = dir('../Anchor2/copyFromPi15/Footstep/');
count = 0;
for i = 1 : length(listing)
    if sum(ismember(listing(i).name,'txtout')) == 9
        count = count + 1;
        signal2{count} = importdata(['../Anchor2/copyFromPi15/Footstep/' listing(i).name]);
    end
end

listing = dir('../Anchor3/copyFromPi15/Footstep/');
count = 0;
for i = 1 : length(listing)
    if sum(ismember(listing(i).name,'txtout')) == 9
        count = count + 1;
        signal3{count} = importdata(['../Anchor3/copyFromPi15/Footstep/' listing(i).name]);
    end
end

listing = dir('../Anchor4/copyFromPi15/Footstep/');
count = 0;
for i = 1 : length(listing)
    if sum(ismember(listing(i).name,'txtout')) == 9
        count = count + 1;
        signal4{count} = importdata(['../Anchor4/copyFromPi15/Footstep/' listing(i).name]);
    end
end

listing = dir('../copyFromMac15/');
count = 0;
for i = 1 : length(listing)
    if sum(ismember(listing(i).name,'txt')) == 3
        count = count + 1;
        base{count} = importdata(['../copyFromMac15/' listing(i).name]);
    end
end

s{1} = [];
for i = 1 : length(signal1)
    s{1} = [s{1}; signal1{i}];
end
s{2} = [];
for i = 1 : length(signal2)
    s{2} = [s{2}; signal2{i}];
end
s{3} = [];
for i = 1 : length(signal3)
    s{3} = [s{3}; signal3{i}];
end
s{4} = [];
for i = 1 : length(signal4)
    s{4} = [s{4}; signal4{i}];
end

tagSig = [];
for i = 1 : length(base)
    tagSig = [tagSig; base{i}];
end

separatorIdx = find(tagSig == 10000);
status = 0;
tag{1}.signal = [];
tag{2}.signal = [];
tag{1}.index = [];
tag{2}.index = [];
timestamp = -1;
for i = 1 : length(separatorIdx)-1
    if separatorIdx(i+1)-separatorIdx(i) == 2 && status == 0
        tagID = tagSig(separatorIdx(i)+1);
        if tagID == 11 || tagID == 12
            status = 1;
        end
    elseif separatorIdx(i+1)-separatorIdx(i) == 2 && status == 1
        timestamp = tagSig(separatorIdx(i)+1);
        status = 2;
    elseif status == 2 && separatorIdx(i+1)-separatorIdx(i) > 10
        tempSig = tagSig(separatorIdx(i)+1:separatorIdx(i+1)-1)';
        tempSig = tempSig(tempSig>0);
        tag{tagID-tagIDOffset}.signal = [tag{tagID-tagIDOffset}.signal, tempSig];
        tempIdx = zeros(1, length(tempSig));
        tempIdx(end) = timestamp;
        for j = length(tempSig)-1:-1:1
            tempIdx(j) = tempIdx(j+1) - 1.0320*10^5;
        end
        tag{tagID-tagIDOffset}.index = [tag{tagID-tagIDOffset}.index, tempIdx];
        
        status = 0;
    end 
end

% figure;
% subplot(2,1,1);
% plot(tag{1}.signal(tag{1}.signal<10^9));
% subplot(2,1,2);
% plot(tag{2}.signal(tag{2}.signal<10^9));

%% assign timestamp to each tag datapoint
for i = 1 : 2
    sigMean = mean(tag{i}.signal);
    tag{i}.signal = tag{i}.signal - sigMean;
    sigRange = max(abs(tag{i}.signal));
    tag{i}.signal = tag{i}.signal / sigRange * 250;
end

figure;
for i = 1 : 4
    subplot(4,1,i);
    plot(s{i}(s{i}<10000));
end

%% assign timestamp to each anchor datapoint
for i = 1 : 4
    tempSig = s{i};
    dIdx = (tempSig == 10000);
    dIdx2 = [0; dIdx(1:end-1)];
    tempSig([find(dIdx>0), find(dIdx2>0)]) = [];
    
    timestamps{i} = tempSig(tempSig>10000);
    timestampIdxs{i} = find(tempSig>10000);
    deltaTime = timestamps{i}(2:end)-timestamps{i}(1:end-1);
    deltaTimeIdx = timestampIdxs{i}(2:end)-timestampIdxs{i}(1:end-1);
    
    timeEva{i} = deltaTime./deltaTimeIdx/10^11;
    meanEva = mean(timeEva{i});
    stdEva = std(timeEva{i});
    if stdEva > 1
        checkDrop = find(abs(timeEva{i}) > meanEva + 3*stdEva);
        s{i}(timestampIdxs{i}(checkDrop-1):end) = [];
    end
    
    timestamps{i} = s{i}(s{i}>10000);
    timestampIdxs{i} = find(s{i}>10000);
    
end

%% alignment
% baseline
unifiedTime = round([timestamps{1}; timestamps{2}; timestamps{3}; timestamps{4}]./10^16);
unifiedTime = unique(unifiedTime);
for i = 1 : length(unifiedTime)
    if ismember(unifiedTime(i),round(timestamps{1}./10^16)) ...
        && ismember(unifiedTime(i),round(timestamps{2}./10^16)) ...
        && ismember(unifiedTime(i),round(timestamps{3}./10^16)) ...
        && ismember(unifiedTime(i),round(timestamps{4}./10^16)) 
        baseline = unifiedTime(i);
        break;
    end
end
for i = 1 : 4
    cuttingPoint{i} = timestampIdxs{i}(round(timestamps{i}./10^16) == baseline);
end
% endline
endline = min(round([timestamps{1}(end),timestamps{2}(end),timestamps{3}(end),timestamps{4}(end)]./10^16));
for i = 1 : 4
    endPoint{i} = timestampIdxs{i}(round(timestamps{i}./10^16) == endline);
end
for i = 1 : 4
    if ~isempty(cuttingPoint{i})
        s{i}(1:cuttingPoint{i}(1)-1) = [];
    end
%     s{i}(endPoint{1}:end) = [];
end

figure;
subplot(4,1,1);
plot(s{1}); title('after cutting');
subplot(4,1,2);
plot(s{2});
subplot(4,1,3);
plot(s{3});
subplot(4,1,4);
plot(s{4});


%% extract ranging signal
% extract signal between each timestamp
for sensorID = 1 : 4
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
                        mean(timeEva{sensorID})*10^11*(distIdx{sensorID}(i) - timestampIdxs{sensorID}(j));
                    break;
                end
            end
        end
        if  distTime{sensorID}(i) == 0
            distTime{sensorID}(i) = timestamps{sensorID}(1)+...
                mean(timeEva{sensorID})*10^11*(distIdx{sensorID}(i) - ...
                timestampIdxs{sensorID}(1));
        end
    end
    d{sensorID,1} = [];
    d{sensorID,2} = [];
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
% figure;
% plot(d{1,2}(:,3),d{1,2}(:,2));hold on;
% plot(d{2,2}(:,3),d{2,2}(:,2));hold on;
% plot(d{3,2}(:,3),d{3,2}(:,2));hold on;
% plot(d{4,2}(:,3),d{4,2}(:,2));hold off;

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
hold off;

save('AllenGabriel.mat','tag','sI','sS');