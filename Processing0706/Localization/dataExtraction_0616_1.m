%  import to check
clear all
close all
clc

listing = dir('./Anchor1/copyFromPi1/Footstep/');
count = 0;
for i = 1 : length(listing)
    if sum(ismember(listing(i).name,'txtout')) == 9
        count = count + 1;
        signal1{count} = importdata(['./Anchor1/copyFromPi1/Footstep/' listing(i).name]);
    end
end

listing = dir('./Anchor2/copyFromPi1/Footstep/');
count = 0;
for i = 1 : length(listing)
    if sum(ismember(listing(i).name,'txtout')) == 9
        count = count + 1;
        signal2{count} = importdata(['./Anchor2/copyFromPi1/Footstep/' listing(i).name]);
    end
end

listing = dir('./Anchor3/copyFromPi1/Footstep/');
count = 0;
for i = 1 : length(listing)
    if sum(ismember(listing(i).name,'txtout')) == 9
        count = count + 1;
        signal3{count} = importdata(['./Anchor3/copyFromPi1/Footstep/' listing(i).name]);
    end
end

listing = dir('./Anchor4/copyFromPi1/Footstep/');
count = 0;
for i = 1 : length(listing)
    if sum(ismember(listing(i).name,'txtout')) == 9
        count = count + 1;
        signal4{count} = importdata(['./Anchor4/copyFromPi1/Footstep/' listing(i).name]);
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

noise1 = s1(1000000:2000000);
noise2 = s2(1000000:2000000);
noise3 = s3(1000000:2000000);
noise4 = s4(1000000:2000000);
noise1 = noise1(noise1 < 10000 & noise1 > 300);
noise2 = noise2(noise2 < 10000 & noise2 > 300);
noise3 = noise3(noise3 < 10000 & noise3 > 300);
noise4 = noise4(noise4 < 10000 & noise4 > 300);

s1 = s1(2800000:4300000);
s2 = s2(2800000:4300000);
s3 = s3(2800000:4300000);
s4 = s4(2800000:4300000);

figure;
subplot(4,1,1);
plot(s1);
subplot(4,1,2);
plot(s2);
subplot(4,1,3);
plot(s3);
subplot(4,1,4);
plot(s4);

% s4 = [signal2{4}; signal2{6}];
timestamp1 = s1(s1>10000);
timestampIdx1 = find(s1>10000);
timestamp2 = s2(s2>10000);
timestampIdx2 = find(s2>10000);
timestamp3 = s3(s3>10000);
timestampIdx3 = find(s3>10000);
timestamp4 = s4(s4>10000);
timestampIdx4 = find(s4>10000);

% % sensor 1
deltaTime = timestamp1(2:end)-timestamp1(1:end-1);
deltaTimeIdx = timestampIdx1(2:end)-timestampIdx1(1:end-1);
timeEva = deltaTime./deltaTimeIdx/10^11;
% deltaVal = mean(timeEva(timeEva < 7));
% errorIdx = find(timeEva >= 7);
% placeHolderLen = calculatePlaceHoldLen( deltaTime, deltaTimeIdx, deltaVal, errorIdx );
% if ~isempty(placeHolderLen)
%     s1 = insertTimePlaceHolder( s1, timestampIdx2(errorIdx+1)-1, placeHolderLen );
% end
% 
% % sensor 2
% deltaTime = timestamp2(2:end)-timestamp2(1:end-1);
% deltaTimeIdx = timestampIdx2(2:end)-timestampIdx2(1:end-1);
% timeEva = deltaTime./deltaTimeIdx/10^11;
% deltaVal = mean(timeEva(timeEva < 7));
% errorIdx = find(timeEva >= 7);
% placeHolderLen = calculatePlaceHoldLen( deltaTime, deltaTimeIdx, deltaVal, errorIdx );
% if ~isempty(placeHolderLen)
%     s2 = insertTimePlaceHolder( s2, timestampIdx2(errorIdx+1)-1, placeHolderLen );
% end
% 
% % sensor 3
% deltaTime = timestamp3(2:end)-timestamp3(1:end-1);
% deltaTimeIdx = timestampIdx3(2:end)-timestampIdx3(1:end-1);
% timeEva = deltaTime./deltaTimeIdx/10^11;
% deltaVal = mean(timeEva(timeEva < 7));
% errorIdx = find(timeEva >= 7);
% placeHolderLen = calculatePlaceHoldLen( deltaTime, deltaTimeIdx, deltaVal, errorIdx );
% if ~isempty(placeHolderLen)
%     s3 = insertTimePlaceHolder( s3, timestampIdx3(errorIdx+1)-1, placeHolderLen );
% end
% 
% % sensor 4
% deltaTime = timestamp4(2:end)-timestamp4(1:end-1);
% deltaTimeIdx = timestampIdx4(2:end)-timestampIdx4(1:end-1);
% timeEva = deltaTime./deltaTimeIdx/10^11;
% deltaVal = mean(timeEva(timeEva < 7));
% errorIdx = find(timeEva >= 7);
% placeHolderLen = calculatePlaceHoldLen( deltaTime, deltaTimeIdx, deltaVal, errorIdx );
% if ~isempty(placeHolderLen)
%     s4 = insertTimePlaceHolder( s4, timestampIdx4(errorIdx+1)-1, placeHolderLen );
% end

timestamp1 = s1(s1>10000);
timestampIdx1 = find(s1>10000);
timestamp2 = s2(s2>10000);
timestampIdx2 = find(s2>10000);
timestamp3 = s3(s3>10000);
timestampIdx3 = find(s3>10000);
timestamp4 = s4(s4>10000);
timestampIdx4 = find(s4>10000);

%% alignment
baseline = round(max([timestamp1(1),timestamp2(1),timestamp3(1),timestamp4(1)])/10^16);
cuttingPoint1 = timestampIdx1(find(round(timestamp1/10^16) == baseline));
cuttingPoint2 = timestampIdx2(find(round(timestamp2/10^16) == baseline));
cuttingPoint3 = timestampIdx3(find(round(timestamp3/10^16) == baseline));
cuttingPoint4 = timestampIdx4(find(round(timestamp4/10^16) == baseline));
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

s1 = s1(1:500000);
s2 = s2(1:500000);
s3 = s3(1:500000);
s4 = s4(1:500000);


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
        d1Timestamp(i) = timestamp1(1)+mean(timeEva)*10^11*(d1Idx(i) - timestampIdx1(1));
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
for i = 1 : length(d2Idx)
    for j = 1 : length(timestampIdx2)-1
        if d2Idx(i) > timestampIdx2(j) && d2Idx(i) < timestampIdx2(j+1)
            d2Timestamp(i) = timestamp2(j)+ ...
                (timestamp2(j+1)-timestamp2(j))/(timestampIdx2(j+1)-timestampIdx2(j))*(d2Idx(i) - timestampIdx2(j));
            break;
        end
    end
    if  d2Timestamp(i) == 0
        d2Timestamp(i) = timestamp2(1)+mean(timeEva)*10^11*(d2Idx(i) - timestampIdx2(1));
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
        d3Timestamp(i) = timestamp3(1)+mean(timeEva)*10^11*(d3Idx(i) - timestampIdx3(1));
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
        d4Timestamp(i) = timestamp4(1)+mean(timeEva)*10^11*(d4Idx(i) - timestampIdx4(1));
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
subplot(3,1,1);
ts = s1_vib(s1_vib < 10000);
plot([1:length(ts)]./10^4,ts);
xlabel('Time (s)');ylabel('Amplitude');
subplot(3,1,2);
plot(d{1,1}(:,1)./10^4,d{1,1}(:,2));
xlabel('Time (s)');ylabel('Distance (m)');
subplot(3,1,3);
plot(d{1,2}(:,1)./10^4,d{1,2}(:,2));
xlabel('Time (s)');ylabel('Distance (m)');
figure;
subplot(3,1,1);
ts = s2_vib(s2_vib < 10000);
plot([1:length(ts)]./10^4,ts);
xlabel('Time (s)');ylabel('Amplitude');
subplot(3,1,2);
plot(d{2,1}(:,1)./10^4,d{2,1}(:,2));
xlabel('Time (s)');ylabel('Distance (m)');
subplot(3,1,3);
plot(d{2,2}(:,1)./10^4,d{2,2}(:,2));
xlabel('Time (s)');ylabel('Distance (m)');
figure;
subplot(3,1,1);
ts = s3_vib(s3_vib < 10000);
plot([1:length(ts)]./10^4,ts);
xlabel('Time (s)');ylabel('Amplitude');
subplot(3,1,2);
plot(d{3,1}(:,1)./10^4,d{3,1}(:,2));
xlabel('Time (s)');ylabel('Distance (m)');
subplot(3,1,3);
plot(d{3,2}(:,1)./10^4,d{3,2}(:,2));
xlabel('Time (s)');ylabel('Distance (m)');
figure;
subplot(3,1,1);
ts = s4_vib(s4_vib < 10000);
plot([1:length(ts)]./10^4,ts);
xlabel('Time (s)');ylabel('Amplitude');
subplot(3,1,2);
plot(d{4,1}(:,1)./10^4,d{4,1}(:,2));
xlabel('Time (s)');ylabel('Distance (m)');
subplot(3,1,3);
plot(d{4,2}(:,1)./10^4,d{4,2}(:,2));
xlabel('Time (s)');ylabel('Distance (m)');

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
plot(s1_vib(s1_vib < 10000));
subplot(4,1,2);
plot(s2_vib(s2_vib < 10000));
subplot(4,1,3);
plot(s3_vib(s3_vib < 10000));
subplot(4,1,4);
plot(s4_vib(s4_vib < 10000));

figure;
subplot(4,1,1);
plot(s1_vib(s1_vib > 10000));
subplot(4,1,2);
plot(s2_vib(s2_vib > 10000));
subplot(4,1,3);
plot(s3_vib(s3_vib > 10000));
subplot(4,1,4);
plot(s4_vib(s4_vib > 10000));

figure;
plot(s1_vib(s1_vib < 10000));hold on;
plot(s2_vib(s2_vib < 10000));hold on;
plot(s3_vib(s3_vib < 10000));hold on;
plot(s4_vib(s4_vib < 10000));hold off;

figure;
plot(s1_vib);hold on;
plot(s2_vib);hold on;
plot(s3_vib);hold on;
plot(s4_vib);hold off;
ylim([0,1200]);

s1_vib(1) = [];
s2_vib(1) = [];
s3_vib(1) = [];
s4_vib(1) = [];

timestamp1 = s1_vib(s1_vib>10000);
timestampIdx1 = find(s1_vib>10000);
timestamp2 = s2_vib(s2_vib>10000);
timestampIdx2 = find(s2_vib>10000);
timestamp3 = s3_vib(s3_vib>10000);
timestampIdx3 = find(s3_vib>10000);
timestamp4 = s4_vib(s4_vib>10000);
timestampIdx4 = find(s4_vib>10000);

unifiedTime = [timestamp1; timestamp2; timestamp3; timestamp4];
unifiedTime = unique(unifiedTime);

for timestampID = 1 : length(unifiedTime)
    timestamp1 = s1_vib(s1_vib>10000);
    timestamp2 = s2_vib(s2_vib>10000);
    timestamp3 = s3_vib(s3_vib>10000);
    timestamp4 = s4_vib(s4_vib>10000);
    
    temp1 = find(s1_vib == unifiedTime(timestampID));
    temp2 = find(s2_vib == unifiedTime(timestampID));
    temp3 = find(s3_vib == unifiedTime(timestampID));
    temp4 = find(s4_vib == unifiedTime(timestampID));
    
    maxIdx = max([temp1, temp2, temp3, temp4]);
    tempSig1 = s1_vib(1:temp1-1);
    tempSig2 = s1_vib(temp1:end);
    s1_vib = [tempSig1; ones(maxIdx-temp1,1)*512; tempSig2];
    
    tempSig1 = s2_vib(1:temp2-1);
    tempSig2 = s2_vib(temp2:end);
    s2_vib = [tempSig1; ones(maxIdx-temp2,1)*512; tempSig2];
    
    tempSig1 = s3_vib(1:temp3-1);
    tempSig2 = s3_vib(temp3:end);
    s3_vib = [tempSig1; ones(maxIdx-temp3,1)*512; tempSig2];
    
    tempSig1 = s4_vib(1:temp4-1);
    tempSig2 = s4_vib(temp4:end);
    s4_vib = [tempSig1; ones(maxIdx-temp4,1)*512; tempSig2];
end

timestamp1 = s1_vib(s1_vib>10000);
timestampIdx1 = find(s1_vib>10000);
timestamp2 = s2_vib(s2_vib>10000);
timestampIdx2 = find(s2_vib>10000);
timestamp3 = s3_vib(s3_vib>10000);
timestampIdx3 = find(s3_vib>10000);
timestamp4 = s4_vib(s4_vib>10000);
timestampIdx4 = find(s4_vib>10000);

figure;
plot(s1_vib);hold on;
plot(s2_vib);hold on;
plot(s3_vib);hold on;
plot(s4_vib);hold off;
ylim([0,1200]);

% sig1Padding = zeros(1,length(unifiedTime));
% sig2Padding = zeros(1,length(unifiedTime));
% sig3Padding = zeros(1,length(unifiedTime));
% sig4Padding = zeros(1,length(unifiedTime));

% for timestampID = 1 : length(unifiedTime)
%     timestampIdxGroup = [timestampIdx1(timestampID), timestampIdx2(timestampID), timestampIdx3(timestampID), timestampIdx4(timestampID)];
%     alignedTimestampIdx = max(timestampIdxGroup);
%     sig1Padding(timestampID) = alignedTimestampIdx - timestampIdx1(timestampID);
%     sig2Padding(timestampID) = alignedTimestampIdx - timestampIdx2(timestampID);
%     sig3Padding(timestampID) = alignedTimestampIdx - timestampIdx3(timestampID);
%     sig4Padding(timestampID) = alignedTimestampIdx - timestampIdx4(timestampID);
% end
% 
% s1_vib = insertTimePlaceHolder(s1_vib, timestampIdx1-1, sig1Padding, 512);
% s2_vib = insertTimePlaceHolder(s2_vib, timestampIdx2-1, sig2Padding, 512);
% s3_vib = insertTimePlaceHolder(s3_vib, timestampIdx3-1, sig3Padding, 512);
% s4_vib = insertTimePlaceHolder(s4_vib, timestampIdx4-1, sig4Padding, 512);
% 

% 
% timestamp1 = s1_vib(s1_vib>10000);
% timestampIdx1 = find(s1_vib>10000);
% timestamp2 = s2_vib(s2_vib>10000);
% timestampIdx2 = find(s2_vib>10000);
% timestamp3 = s3_vib(s3_vib>10000);
% timestampIdx3 = find(s3_vib>10000);
% timestamp4 = s4_vib(s4_vib>10000);
% timestampIdx4 = find(s4_vib>10000);

% s1_vib(s1_vib > 10000) = 2000;
% s2_vib(s2_vib > 10000) = 2000;
% s3_vib(s3_vib > 10000) = 2000;
% s4_vib(s4_vib > 10000) = 2000;
% figure;
% plot(s1_vib);hold on;
% plot(s2_vib);hold on;
% plot(s3_vib);hold on;
% plot(s4_vib);hold off;

save('ranging1.mat','s1','s2','s3','s4',...
    'noise1','noise2','noise3','noise4',...
    's1_vib','s2_vib','s3_vib','s4_vib');
