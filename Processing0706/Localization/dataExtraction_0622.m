%  import to check
clear all
close all
clc

listing = dir('../Anchor1/copyFromPi8/Footstep/');
count = 0;
for i = 1 : length(listing)
    if sum(ismember(listing(i).name,'txtout')) == 9
        count = count + 1;
        signal1{count} = importdata(['../Anchor1/copyFromPi8/Footstep/' listing(i).name]);
    end
end

listing = dir('../Anchor2/copyFromPi8/Footstep/');
count = 0;
for i = 1 : length(listing)
    if sum(ismember(listing(i).name,'txtout')) == 9
        count = count + 1;
        signal2{count} = importdata(['../Anchor2/copyFromPi8/Footstep/' listing(i).name]);
    end
end

listing = dir('../Anchor3/copyFromPi8/Footstep/');
count = 0;
for i = 1 : length(listing)
    if sum(ismember(listing(i).name,'txtout')) == 9
        count = count + 1;
        signal3{count} = importdata(['../Anchor3/copyFromPi8/Footstep/' listing(i).name]);
    end
end

listing = dir('../Anchor4/copyFromPi8/Footstep/');
count = 0;
for i = 1 : length(listing)
    if sum(ismember(listing(i).name,'txtout')) == 9
        count = count + 1;
        signal4{count} = importdata(['../Anchor4/copyFromPi8/Footstep/' listing(i).name]);
    end
end

s1 = [];
for i = 1 : length(signal1)
    s1 = [s1; signal1{i}];
end
s2 = [];
for i = 1 : length(signal2)
    s2 = [s2; signal2{i}];
end
s3 = [];
for i = 1 : length(signal3)
    s3 = [s3; signal3{i}];
end
s4 = [];
for i = 1 : length(signal4)
    s4 = [s4; signal4{i}];
end

figure;
subplot(4,1,1);
plot(s1);
subplot(4,1,2);
plot(s2);
subplot(4,1,3);
plot(s3);
subplot(4,1,4);
plot(s4);

timestamp1 = s1(s1>10000);
timestampIdx1 = find(s1>10000);
timestamp2 = s2(s2>10000);
timestampIdx2 = find(s2>10000);
timestamp3 = s3(s3>10000);
timestampIdx3 = find(s3>10000);
timestamp4 = s4(s4>10000);
timestampIdx4 = find(s4>10000);
% 
% % sensor 1
deltaTime = timestamp1(2:end)-timestamp1(1:end-1);
deltaTimeIdx = timestampIdx1(2:end)-timestampIdx1(1:end-1);
timeEva1 = deltaTime./deltaTimeIdx/10^11;
deltaVal = mean(timeEva1(timeEva1 < 7));
errorIdx = find(timeEva1 >= 7);
placeHolderLen = calculatePlaceHoldLen( deltaTime, deltaTimeIdx, deltaVal, errorIdx );
if ~isempty(placeHolderLen)
    s1 = insertTimePlaceHolder( s1, timestampIdx2(errorIdx+1)-1, placeHolderLen, 512 );
end
% 
% % sensor 2
deltaTime = timestamp2(2:end)-timestamp2(1:end-1);
deltaTimeIdx = timestampIdx2(2:end)-timestampIdx2(1:end-1);
timeEva2 = deltaTime./deltaTimeIdx/10^11;
deltaVal = mean(timeEva2(timeEva2 < 7));
errorIdx = find(timeEva2 >= 7);
placeHolderLen = calculatePlaceHoldLen( deltaTime, deltaTimeIdx, deltaVal, errorIdx );
if ~isempty(placeHolderLen)
    s2 = insertTimePlaceHolder( s2, timestampIdx2(errorIdx+1)-1, placeHolderLen, 512 );
end
% 
% % sensor 3
deltaTime = timestamp3(2:end)-timestamp3(1:end-1);
deltaTimeIdx = timestampIdx3(2:end)-timestampIdx3(1:end-1);
timeEva3 = deltaTime./deltaTimeIdx/10^11;
deltaVal = mean(timeEva3(timeEva3 < 7));
errorIdx = find(timeEva3 >= 7);
placeHolderLen = calculatePlaceHoldLen( deltaTime, deltaTimeIdx, deltaVal, errorIdx );
if ~isempty(placeHolderLen)
    s3 = insertTimePlaceHolder( s3, timestampIdx3(errorIdx+1)-1, placeHolderLen, 512 );
end
% 
% % sensor 4
deltaTime = timestamp4(2:end)-timestamp4(1:end-1);
deltaTimeIdx = timestampIdx4(2:end)-timestampIdx4(1:end-1);
timeEva4 = deltaTime./deltaTimeIdx/10^11;
deltaVal = mean(timeEva4(timeEva4 < 7));
errorIdx = find(timeEva4 >= 7);
placeHolderLen = calculatePlaceHoldLen( deltaTime, deltaTimeIdx, deltaVal, errorIdx );
if ~isempty(placeHolderLen)
    s4 = insertTimePlaceHolder( s4, timestampIdx4(errorIdx+1)-1, placeHolderLen, 512 );
end

timestamp1 = s1(s1>10000);
timestampIdx1 = find(s1>10000);
timestamp2 = s2(s2>10000);
timestampIdx2 = find(s2>10000);
timestamp3 = s3(s3>10000);
timestampIdx3 = find(s3>10000);
timestamp4 = s4(s4>10000);
timestampIdx4 = find(s4>10000);

%% alignment
unifiedTime = round([timestamp1; timestamp2; timestamp3; timestamp4]./10^16);
unifiedTime = unique(unifiedTime);
for i = 1 : length(unifiedTime)
    if ismember(unifiedTime(i),round(timestamp1./10^16)) ...
        && ismember(unifiedTime(i),round(timestamp2./10^16)) ...
        && ismember(unifiedTime(i),round(timestamp3./10^16)) ...
        && ismember(unifiedTime(i),round(timestamp4./10^16)) 
        baseline = unifiedTime(i);
        break;
    end
end

cuttingPoint1 = timestampIdx1(round(timestamp1./10^16) == baseline);
cuttingPoint2 = timestampIdx2(round(timestamp2./10^16) == baseline);
cuttingPoint3 = timestampIdx3(round(timestamp3./10^16) == baseline);
cuttingPoint4 = timestampIdx4(round(timestamp4./10^16) == baseline);

endline = min(round([timestamp1(end),timestamp2(end),timestamp3(end),timestamp4(end)]./10^16));

endPoint1 = timestampIdx1(round(timestamp1./10^16) == endline);
endPoint2 = timestampIdx2(round(timestamp2./10^16) == endline);
endPoint3 = timestampIdx3(round(timestamp3./10^16) == endline);
endPoint4 = timestampIdx4(round(timestamp4./10^16) == endline);

if ~isempty(cuttingPoint1)
    s1(1:cuttingPoint1(1)-1) = [];
end
if ~isempty(cuttingPoint2)
    s2(1:cuttingPoint2(1)-1) = [];
end
if ~isempty(cuttingPoint3)
    s3(1:cuttingPoint3(1)-1) = [];
end
if ~isempty(cuttingPoint4)
    s4(1:cuttingPoint4(1)-1) = [];
end
% 
s1(endPoint1:end) = [];
s2(endPoint2:end) = [];
s3(endPoint3:end) = [];
s4(endPoint4:end) = [];

figure;
subplot(4,1,1);
plot(s1);
subplot(4,1,2);
plot(s2);
subplot(4,1,3);
plot(s3);
subplot(4,1,4);
plot(s4);


%% alignment
% extract signal between each timestamp
dIdx = (s1 == 10000);
dIdx2 = [0; dIdx(1:end-1)];

distance1 = s1.*(dIdx2);
d1 = distance1(distance1>0);
d1Idx = find(distance1>0);
d1Timestamp = zeros(1,length(d1Idx));
for i = 1 : length(d1Idx)
    for j = 1 : length(timestampIdx1)-1
        if d1Idx(i) > timestampIdx1(j) && d1Idx(i) < timestampIdx1(j+1)
            d1Timestamp(i) = timestamp1(j)+ ...
                (timestamp1(j+1)-timestamp1(j))/(timestampIdx1(j+1)-timestampIdx1(j))*(d1Idx(i) - timestampIdx1(j));
            break;
        end
    end
    if  d1Timestamp(i) == 0
        fprintf('Meh');
        for j = 1 : length(timestampIdx1)
            if d1Idx(i) > timestampIdx1(j)
                d1Timestamp(i) = timestamp1(j)+mean(timeEva1)*10^11*(d1Idx(i) - timestampIdx1(j));
                break;
            end
        end
    end
    
    if  d1Timestamp(i) == 0
        d1Timestamp(i) = timestamp1(1)+mean(timeEva1)*10^11*(d1Idx(i) - timestampIdx1(1));
    end
end
d{1,1} = [];
d{1,2} = [];
for i = 1 : 2 : length(d1)
    if d1(i) == 11 && d1(i+1) < 100
        d{1,1} = [d{1,1}; d1Idx(i+1) d1(i+1) d1Timestamp(i+1)];
    elseif d1(i) == 12 && d1(i+1) < 100
        d{1,2} = [d{1,2}; d1Idx(i+1) d1(i+1) d1Timestamp(i+1)];
    end
end
dIdx = find(s1 == 10000);
dIdx2 = dIdx+1;
s1_vib = s1;
s1_vib([dIdx;dIdx2]) = [];

dIdx = (s2 == 10000);
dIdx2 = [0; dIdx(1:end-1)];
distance2 = s2.*(dIdx2);
d2 = distance2(distance2>0);
d2Idx = find(distance2>0);
d2Timestamp = zeros(1,length(d2Idx));
% blackList = [];
for i = 1 : length(d2Idx)
    for j = 1 : length(timestampIdx2)-1
        if d2Idx(i) > timestampIdx2(j) && d2Idx(i) < timestampIdx2(j+1)
            d2Timestamp(i) = timestamp2(j)+ ...
                (timestamp2(j+1)-timestamp2(j))/(timestampIdx2(j+1)-timestampIdx2(j))*(d2Idx(i) - timestampIdx2(j));
            break;
        end
    end
    if  d2Timestamp(i) == 0
        fprintf('Meh');
        for j = 1 : length(timestampIdx2)
            if d2Idx(i) > timestampIdx2(j)
                d2Timestamp(i) = timestamp2(j)+mean(timeEva2)*10^11*(d2Idx(i) - timestampIdx2(j));
                break;
            end
        end
    end
    if  d2Timestamp(i) == 0
        d2Timestamp(i) = timestamp2(1)+mean(timeEva2)*10^11*(d2Idx(i) - timestampIdx2(1));
    end
end
d{2,1} = [];
d{2,2} = [];
for i = 1 : 2 : length(d2)
    if d2(i) == 11 && d2(i+1) < 100
        d{2,1} = [d{2,1}; d2Idx(i+1) d2(i+1) d2Timestamp(i+1)];
    elseif d2(i) == 12 && d2(i+1) < 100
        d{2,2} = [d{2,2}; d2Idx(i+1) d2(i+1) d2Timestamp(i+1)];
    end
end

dIdx = find(s2 == 10000);
dIdx2 = dIdx+1;
s2_vib = s2;
s2_vib([dIdx;dIdx2]) = [];

dIdx = (s3 == 10000);
dIdx2 = [0; dIdx(1:end-1)];
distance3 = s3.*(dIdx2);
d3 = distance3(distance3>0);
d3Idx = find(distance3>0);
d3Timestamp = zeros(1,length(d3Idx));
for i = 1 : length(d3Idx)
    for j = 1 : length(timestampIdx3)-1
        if d3Idx(i) > timestampIdx3(j) && d3Idx(i) < timestampIdx3(j+1)
            d3Timestamp(i) = timestamp3(j)+ ...
                (timestamp3(j+1)-timestamp3(j))/(timestampIdx3(j+1)-timestampIdx3(j))*(d3Idx(i) - timestampIdx3(j));
            break;
        end
    end
    
    if  d3Timestamp(i) == 0
        fprintf('Meh');
        for j = 1 : length(timestampIdx3)
            if d3Idx(i) > timestampIdx3(j)
                d3Timestamp(i) = timestamp3(j)+mean(timeEva3)*10^11*(d3Idx(i) - timestampIdx3(j));
                break;
            end
        end
    end
    if  d3Timestamp(i) == 0
        d3Timestamp(i) = timestamp3(1)+mean(timeEva3)*10^11*(d3Idx(i) - timestampIdx3(1));
    end
end
d{3,1} = [];
d{3,2} = [];
for i = 1 : 2 : length(d3)
    if d3(i) == 11 && d3(i+1) < 100
        d{3,1} = [d{3,1}; d3Idx(i+1) d3(i+1) d3Timestamp(i+1)];
    elseif d3(i) == 12 && d3(i+1) < 100
        d{3,2} = [d{3,2}; d3Idx(i+1) d3(i+1) d3Timestamp(i+1)];
    end
end
dIdx = find(s3 == 10000);
dIdx2 = dIdx+1;
s3_vib = s3;
s3_vib([dIdx;dIdx2]) = [];

dIdx = (s4 == 10000);
dIdx2 = [0; dIdx(1:end-1)];
distance4 = s4.*(dIdx2);
d4 = distance4(distance4>0);
d4Idx = find(distance4>0);
d4Timestamp = zeros(1,length(d4Idx));
for i = 1 : length(d4Idx)
    for j = 1 : length(timestampIdx4)-1
        if d4Idx(i) > timestampIdx4(j) && d4Idx(i) < timestampIdx4(j+1)
            d4Timestamp(i) = timestamp4(j)+ ...
                (timestamp4(j+1)-timestamp4(j))/(timestampIdx4(j+1)-timestampIdx4(j))*(d4Idx(i) - timestampIdx4(j));
            break;
        end
    end
    
    if  d4Timestamp(i) == 0
        fprintf('Meh');
        for j = 1 : length(timestampIdx4)
            if d4Idx(i) > timestampIdx4(j)
                d4Timestamp(i) = timestamp4(j)+mean(timeEva4)*10^11*(d4Idx(i) - timestampIdx4(j));
                break;
            end
        end
    end
    
    if  d4Timestamp(i) == 0
        d4Timestamp(i) = timestamp4(1)+mean(timeEva4)*10^11*(d4Idx(i) - timestampIdx4(1));
    end
end
d{4,1} = [];
d{4,2} = [];
for i = 1 : 2 : length(d4)
    if d4(i) == 11 && d4(i+1) < 100
        d{4,1} = [d{4,1}; d4Idx(i+1) d4(i+1) d4Timestamp(i+1)];
    elseif d4(i) == 12 && d4(i+1) < 100
        d{4,2} = [d{4,2}; d4Idx(i+1) d4(i+1) d4Timestamp(i+1)];
    end
end
dIdx = find(s4 == 10000);
dIdx2 = dIdx+1;
s4_vib = s4;
s4_vib([dIdx;dIdx2]) = [];

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
plot(s1_vib);
subplot(4,1,2);
plot(s2_vib);
subplot(4,1,3);
plot(s3_vib);
subplot(4,1,4);
plot(s4_vib);

% fine grain alignment
timestamp1 = s1_vib(s1_vib>10000);
timestampIdx1 = find(s1_vib>10000);
timestamp2 = s2_vib(s2_vib>10000);
timestampIdx2 = find(s2_vib>10000);
timestamp3 = s3_vib(s3_vib>10000);
timestampIdx3 = find(s3_vib>10000);
timestamp4 = s4_vib(s4_vib>10000);
timestampIdx4 = find(s4_vib>10000);

sigVib{1} = s1_vib;
sigVib{2} = s2_vib;
sigVib{3} = s3_vib;
sigVib{4} = s4_vib;

timestamps{1} = timestamp1;
timestamps{2} = timestamp2;
timestamps{3} = timestamp3;
timestamps{4} = timestamp4;

timestampIdxs{1} = timestampIdx1;
timestampIdxs{2} = timestampIdx2;
timestampIdxs{3} = timestampIdx3;
timestampIdxs{4} = timestampIdx4;

% alignedSig = timestampAlignment4C( sigVib, timestamps, timestampIdxs );

[sS{1}, sI{1}] = timestampAlignment( s1_vib, timestamp1, timestampIdx1 );
[sS{2}, sI{2}] = timestampAlignment( s2_vib, timestamp2, timestampIdx2 );
[sS{3}, sI{3}] = timestampAlignment( s3_vib, timestamp3, timestampIdx3 );
[sS{4}, sI{4}] = timestampAlignment( s4_vib, timestamp4, timestampIdx4 );
% 
for i = 1 : 4
    sS{i} = sS{i}-mean(sS{i});
end
sS{3} = -sS{3};
sS{4} = -sS{4};

figure;
plot(sI{1}, sS{1});hold on;
plot(sI{2}, sS{2});hold on;
plot(sI{3}, sS{3});hold on;
plot(sI{4}, sS{4});hold off;
ylim([-700,700]);

for i = 1 : 4
    sI{i}(1:1250000) = [];
    sS{i}(1:1250000) = [];
end

figure;
plot(sI{1}, sS{1});hold on;
plot(sI{2}, sS{2});hold on;
plot(sI{3}, sS{3});hold on;
plot(sI{4}, sS{4});hold off;
ylim([-700,700]);


noiseTime = [2.43*10^18, 2.48*10^18];
for i = 1 : 4
    [~, startIdx] = min(abs(sI{i}-noiseTime(1)));
    [~, stopIdx] = min(abs(sI{i}-noiseTime(2)));
    noiseIdx{i} = sI{i}(startIdx:stopIdx);
    noiseSig{i} = sS{i}(startIdx:stopIdx);
end

Fs = 15000;            % Sampling frequency
T = 1/Fs;             % Sampling period
L = length(noiseSig{i});             % Length of signal
t = (0:L-1)*T;        % Time vector
n = 2^nextpow2(L);
Y = fft(noiseSig{i},n);
f = Fs*(0:(n/2))/n;
P = abs(Y/n);
figure;
plot(f,P(1:n/2+1))
title('Frequency Domain')
xlabel('Frequency (f)')
ylabel('|P(f)|')


locTime{1} = [2.495*10^18, 2.51*10^18];
locTime{2} = [2.54*10^18, 2.56*10^18];
locTime{3} = [2.5928*10^18, 2.605*10^18];
locTime{4} = [2.63*10^18, 2.65*10^18];
locTime{5} = [2.66*10^18, 2.68*10^18];
locTime{6} = [2.7*10^18, 2.71*10^18];
locTime{7} = [2.725*10^18, 2.74*10^18];
locTime{8} = [2.762*10^18, 2.775*10^18];
locTime{9} = [2.7968*10^18, 2.804*10^18];
locTime{10} = [2.83*10^18, 2.838*10^18];
locTime{11} = [2.876*10^18, 2.885*10^18];
locTime{12} = [2.905*10^18, 2.913*10^18];
locTime{13} = [2.942*10^18, 2.95*10^18];
locTime{14} = [2.978*10^18, 2.987*10^18];
locTime{15} = [3.005*10^18, 3.02*10^18];
locTime{16} = [3.06*10^18, 3.07*10^18];
locTime{17} = [3.106*10^18, 3.115*10^18];
locTime{18} = [3.155*10^18, 3.165*10^18];

energyRatio = zeros(18,4);

for j = 1 :18
%     figure;
    for i = 1 : 4
        [~, startIdx] = min(abs(sI{i}-locTime{j}(1)));
        [~, stopIdx] = min(abs(sI{i}-locTime{j}(2)));
        locsI{j,i} = sI{i}(startIdx:stopIdx);
        locsS{j,i} = sS{i}(startIdx:stopIdx);
%         plot(locsI{j,i},locsS{j,i}); hold on;
        sigEnergy1 = sum(locsS{j,1}.*locsS{j,1});
        sigEnergy2 = sum(locsS{j,i}.*locsS{j,i});
        energyRatio(j,i) = sigEnergy2/sigEnergy1;
    end
%     hold off;
end
distRaio = 1./energyRatio;
imagesc(energyRatio);

for j = 1 :18
%     figure;
    for i = 1 : 4
        [~, startIdx] = min(abs(d{i,1}(:,3)-locTime{j}(1)));
        [~, stopIdx] = min(abs(d{i,1}(:,3)-locTime{j}(2)));
        dI{j,i} = d{i,1}(startIdx:stopIdx,1);
        dS{j,i} = d{i,1}(startIdx:stopIdx,2);
%         plot(dI{j,i},dS{j,i}); hold on;
    end
%     hold off;
end

sensorLoc{1} = [0,0];
sensorLoc{2} = [12,0];
sensorLoc{3} = [12,10];
sensorLoc{4} = [0,10];

impulseLoc{1} = [-4,0];
impulseLoc{2} = [-2,0];
impulseLoc{3} = [2,0];
impulseLoc{4} = [4,0];
impulseLoc{5} = [6,0];
impulseLoc{6} = [8,0];
impulseLoc{7} = [10,0];
impulseLoc{8} = [12,2];
impulseLoc{9} = [12,4];
impulseLoc{10} = [12,6];
impulseLoc{11} = [12,8];
impulseLoc{12} = [10,10];
impulseLoc{13} = [8,10];
impulseLoc{14} = [6,10];
impulseLoc{15} = [4,10];
impulseLoc{16} = [2,10];
impulseLoc{17} = [-2,10];
impulseLoc{18} = [-4,10];

save('alignedSig.mat','sI','sS','locsI','locsS','sensorLoc','impulseLoc','noiseIdx','noiseSig','dI','dS','energyRatio');

% save('alignedSig.mat','sI','sS','locsI','locsS','sensorLoc','impulseLoc','noiseIdx','noiseSig','energyRatio');

[mV, mI] = min(energyRatio');
for j = 1 : 18
    figure;
    sensorSet = [1:4];
    sensorSet(sensorSet == mI(j)) = [];
    peakTime{j} = [-1,-1,-1,-1];
    sensorCombo{j} = [];
    for i = sensorSet%1 : 4
        sensorCombo{j} = [sensorCombo{j}; sensorLoc{i}];
        plot(locsI{j,i},locsS{j,i});hold on;
        if max(abs(locsS{j,i})) > 450
            firstPeak = findFirstPeak(locsS{j,i}, max(abs(locsS{j,i}))/4);
    %     elseif max(abs(locsS{j,i})) < 100
    %         firstPeak = findFirstPeak(locsS{j,i}, max(abs(locsS{j,i}))/10);
        else 
            firstPeak = findFirstPeak(locsS{j,i}, max(abs(locsS{j,i}))/10);
        end
        peakTime{j}(i) = locsI{j,i}(firstPeak);
        scatter(locsI{j,i}(firstPeak),locsS{j,i}(firstPeak));
    end
    hold off;
end

%% pick top strongest three to calculate
errorRecord = [];
for speed = 1060
    figure;
    tempError = 0;
    locCount = 0;
    for i = 1 : 18
        [locX, locY] = simpleMultilateration( peakTime{i}, sensorCombo{i}, speed );
        if locX == -1 && locY == -1
            continue;
        end
        scatter(locX, locY);hold on;
        plot([impulseLoc{i}(1),locX],[impulseLoc{i}(2),locY]);hold on;
        tempError = tempError + sqrt((impulseLoc{i}(1)-locX)^2+(impulseLoc{i}(2)-locY)^2);
        locCount = locCount + 1;
    end
    hold off;
    errorRecord = [errorRecord; speed, tempError, locCount];
end



figure;
for i = 1 : 18
    scatter(i,mean(dS{i,1}),'b');hold on;
    gt_dist = sqrt((impulseLoc{i}(1)-sensorLoc{1}(1))^2+(impulseLoc{i}(2)-sensorLoc{1}(2))^2);
    scatter(i,gt_dist/3,'r');hold on;
end
hold off;

figure;
for i = 1 : 18
    scatter(i,mean(dS{i,2}),'b');hold on;
    gt_dist = sqrt((impulseLoc{i}(1)-sensorLoc{2}(1))^2+(impulseLoc{i}(2)-sensorLoc{2}(2))^2);
    scatter(i,gt_dist/3,'r');hold on;
end
hold off;

figure;
for i = 1 : 18
    scatter(i,mean(dS{i,3}),'b');hold on;
    gt_dist = sqrt((impulseLoc{i}(1)-sensorLoc{3}(1))^2+(impulseLoc{i}(2)-sensorLoc{3}(2))^2);
    scatter(i,gt_dist/3,'r');hold on;
end
hold off;

figure;
for i = 1 : 18
    scatter(i,mean(dS{i,4}),'b');hold on;
    gt_dist = sqrt((impulseLoc{i}(1)-sensorLoc{4}(1))^2+(impulseLoc{i}(2)-sensorLoc{4}(2))^2);
    scatter(i,gt_dist/3,'r');hold on;
end
hold off;
% 
% save('alignedSig.mat','sI','sS','locsI','locsS','sensorLoc','impulseLoc','noiseIdx','noiseSig','dI','dS');

