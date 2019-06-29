%  import to check
% clear all
% close all
% clc

listing = dir('./Anchor1/copyFromPi/16/');
count = 0;
for i = 1 : length(listing)
    if sum(ismember(listing(i).name,'txtout')) == 9
        count = count + 1;
        signal1{count} = importdata(['./Anchor1/copyFromPi/16/' listing(i).name]);
    end
end

listing = dir('./Anchor2/copyFromPi/16/');
count = 0;
for i = 1 : length(listing)
    if sum(ismember(listing(i).name,'txtout')) == 9
        count = count + 1;
        signal2{count} = importdata(['./Anchor2/copyFromPi/16/' listing(i).name]);
    end
end

listing = dir('./Anchor3/copyFromPi/16/');
count = 0;
for i = 1 : length(listing)
    if sum(ismember(listing(i).name,'txtout')) == 9
        count = count + 1;
        signal3{count} = importdata(['./Anchor3/copyFromPi/16/' listing(i).name]);
    end
end

listing = dir('./Anchor4/copyFromPi/16/');
count = 0;
for i = 1 : length(listing)
    if sum(ismember(listing(i).name,'txtout')) == 9
        count = count + 1;
        signal4{count} = importdata(['./Anchor4/copyFromPi/16/' listing(i).name]);
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
    if d4(i) == 11 && d4(i+1) < 10
        d{4,1} = [d{4,1}; d4Idx(i+1) d4(i+1) d4Timestamp(i+1)];
    elseif d4(i) == 12 && d4(i+1) < 10
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
figure;
plot(d{1,2}(:,3),d{1,2}(:,2));hold on;
plot(d{2,2}(:,3),d{2,2}(:,2));hold on;
plot(d{3,2}(:,3),d{3,2}(:,2));hold on;
plot(d{4,2}(:,3),d{4,2}(:,2));hold off;

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

[s1S, s1I] = timestampAlignment( s1_vib, timestamp1, timestampIdx1 );
[s2S, s2I] = timestampAlignment( s2_vib, timestamp2, timestampIdx2 );
[s3S, s3I] = timestampAlignment( s3_vib, timestamp3, timestampIdx3 );
[s4S, s4I] = timestampAlignment( s4_vib, timestamp4, timestampIdx4 );
% 
s1S = s1S-mean(s1S);
s2S = s2S-mean(s2S);
s3S = s3S-mean(s3S);
s4S = s4S-mean(s4S);
s3S = -s3S;
s4S = -s4S;

figure;
plot(s1I, s1S);hold on;
plot(s2I, s2S);hold on;
plot(s3I, s3S);hold on;
plot(s4I, s4S);hold off;
ylim([-700,700]);

save('alignedSig.mat','s1I','s1S','s2I','s2S','s3I','s3S','s4I','s4S');


