%  import to check
clear all
close all
clc

listing = dir('../Anchor1/copyFromPi/Footstep/');
count = 0;
for i = 1 : length(listing)
    if sum(ismember(listing(i).name,'txtout')) == 9
        count = count + 1;
        signal1{count} = importdata(['../Anchor1/copyFromPi/Footstep/' listing(i).name]);
    end
end

listing = dir('../Anchor2/copyFromPi/Footstep/');
count = 0;
for i = 1 : length(listing)
    if sum(ismember(listing(i).name,'txtout')) == 9
        count = count + 1;
        signal2{count} = importdata(['../Anchor2/copyFromPi/Footstep/' listing(i).name]);
    end
end

listing = dir('../Anchor3/copyFromPi/Footstep/');
count = 0;
for i = 1 : length(listing)
    if sum(ismember(listing(i).name,'txtout')) == 9
        count = count + 1;
        signal3{count} = importdata(['../Anchor3/copyFromPi/Footstep/' listing(i).name]);
    end
end

listing = dir('../Anchor4/copyFromPi/Footstep/');
count = 0;
for i = 1 : length(listing)
    if sum(ismember(listing(i).name,'txtout')) == 9
        count = count + 1;
        signal4{count} = importdata(['../Anchor4/copyFromPi/Footstep/' listing(i).name]);
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


figure;
for i = 1 : 4
    subplot(4,1,i);
    plot(s{i}(s{i}<10000));
end

for i = 1 : 4
    timestamps{i} = s{i}(s{i}>10000);
    timestampIdxs{i} = find(s{i}>10000);
    deltaTime = timestamps{i}(2:end)-timestamps{i}(1:end-1);
    deltaTimeIdx = timestampIdxs{i}(2:end)-timestampIdxs{i}(1:end-1);
    timeEva{i} = deltaTime./deltaTimeIdx/10^11;
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
% for i = 1 : 4
%     if ~isempty(cuttingPoint{i})
%         s{i}(1:cuttingPoint{i}(1)-1) = [];
%     end
% %     s{i}(endPoint{1}:end) = [];
% end

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
plot(sI{1}(1:100:end), sS{1}(1:100:end));hold on;
plot(sI{2}(1:100:end), sS{2}(1:100:end));hold on;
plot(sI{3}(1:100:end), sS{3}(1:100:end));hold on;
plot(sI{4}(1:100:end), sS{4}(1:100:end));hold off;

% for i = 1 : 100000 : length(sI{1})-100000
%     figure;
%     plot(sI{1}(i:i+100000), sS{1}(i:i+100000));hold on;
%     plot(sI{2}(i:i+100000), sS{2}(i:i+100000));hold on;
%     plot(sI{3}(i:i+100000), sS{3}(i:i+100000));hold on;
%     plot(sI{4}(i:i+100000), sS{4}(i:i+100000));hold off;
%     ylim([-700,700]);
% end
% 
% 
% noiseTime = [2.43*10^18, 2.48*10^18];
% for i = 1 : 4
%     [~, startIdx] = min(abs(sI{i}-noiseTime(1)));
%     [~, stopIdx] = min(abs(sI{i}-noiseTime(2)));
%     noiseIdx{i} = sI{i}(startIdx:stopIdx);
%     noiseSig{i} = sS{i}(startIdx:stopIdx);
% end
% 
% Fs = 15000;            % Sampling frequency
% T = 1/Fs;             % Sampling period
% L = length(noiseSig{i});             % Length of signal
% t = (0:L-1)*T;        % Time vector
% n = 2^nextpow2(L);
% Y = fft(noiseSig{i},n);
% f = Fs*(0:(n/2))/n;
% P = abs(Y/n);
% figure;
% plot(f,P(1:n/2+1))
% title('Frequency Domain')
% xlabel('Frequency (f)')
% ylabel('|P(f)|')
% 
% 
%% TO ADD
% impulseNum = 7;
% locTime{1} = [1.742*10^18, 1.752*10^18];
% locTime{2} = [1.865*10^18, 1.875*10^18];
% locTime{3} = [1.935*10^18, 1.945*10^18];
% locTime{4} = [2.074*10^18, 2.08*10^18];
% locTime{5} = [2.333*10^18, 2.339*10^18];
% locTime{6} = [2.723*10^18, 2.729*10^18];
% locTime{7} = [3.005*10^18, 3.065*10^18];
% % 
% 
% energyRatio = zeros(impulseNum,4);
% 
% for j = 1 :impulseNum
% %     figure;
%     for i = 1 : 4
%         [~, startIdx] = min(abs(sI{i}-locTime{j}(1)));
%         [~, stopIdx] = min(abs(sI{i}-locTime{j}(2)));
%         sRange{j,i} = [startIdx, stopIdx];
%         locsI{j,i} = sI{i}(startIdx:stopIdx);
%         locsS{j,i} = sS{i}(startIdx:stopIdx);
% %         plot(locsI{j,i},locsS{j,i}); hold on;
%         sigEnergy1 = sum(locsS{j,1}.*locsS{j,1});
%         sigEnergy2 = sum(locsS{j,i}.*locsS{j,i});
%         energyRatio(j,i) = sigEnergy2/sigEnergy1;
%     end
% %     hold off;
% end
% distRaio = 1./energyRatio;
% imagesc(energyRatio);
% 
% for j = 1 impulseNum
% %     figure;
%     for i = 1 : 4
%         [~, startIdx] = min(abs(d{i,1}(:,3)-locTime{j}(1)));
%         [~, stopIdx] = min(abs(d{i,1}(:,3)-locTime{j}(2)));
% %         if d{i,1}(startIdx-2,3) < locTime{j}(1) &&  d{i,1}(startIdx+2,3) > locTime{j}(1)
% %             fprintf('startIdx ok\n');
% %         end
% %         if d{i,1}(stopIdx-2,3) < locTime{j}(2) &&  d{i,1}(stopIdx+2,3) > locTime{j}(2)
% %             fprintf('stopIdx ok\n');
% %         end
%         dRange{j,i} = [startIdx, stopIdx];
%         dI{j,i} = d{i,1}(startIdx:stopIdx,1);
%         dS{j,i} = d{i,1}(startIdx:stopIdx,2);
%         dT{j,i} = d{i,1}(startIdx:stopIdx,3);
% %         plot(dI{j,i},dS{j,i}); hold on;
%     end
% %     hold off;
% end
% for i = 1 : 4
%     figure;
%     subplot(2,1,1);
%     plot(d{i,1}(:,3),d{i,1}(:,2));hold on;
%     for j = 1 : impulseNum
%         plot([dT{j,i}(1),dT{j,i}(1)],[0,10],'b');
% %         gt_dist = sqrt((impulseLoc{j}(1)-sensorLoc{i}(1))^2+(impulseLoc{j}(2)-sensorLoc{i}(2))^2);
% %         scatter(sI{i}(sRange{j,i}(1)),gt_dist/3,'r');hold on;
%     end
%     hold off;
%     axis tight
%     xlim([4*10^18,7*10^18]);
%     subplot(2,1,2);
%     plot(sI{i},sS{i});hold on;
%     for j = 1 : impulseNum
%         plot(sI{i}(sRange{j,i}),[-512,512]); 
%     end
%     hold off;
%     axis tight
% %     xlim([4*10^18,7*10^18]);
% end
% % 
% sensorLoc{1} = [0,0];
% sensorLoc{2} = [12,0];
% sensorLoc{3} = [12,10];
% sensorLoc{4} = [0,10];
% 
% impulseLoc{1} = [-4,0];
% impulseLoc{2} = [-2,0];
% impulseLoc{3} = [2,0];
% impulseLoc{4} = [4,0];
% impulseLoc{5} = [6,0];
% impulseLoc{6} = [8,0];
% impulseLoc{7} = [10,0];
% impulseLoc{8} = [12,2];
% impulseLoc{9} = [12,4];
% impulseLoc{10} = [12,6];
% impulseLoc{11} = [12,8];
% impulseLoc{12} = [10,10];
% impulseLoc{13} = [8,10];
% impulseLoc{14} = [6,10];
% impulseLoc{15} = [4,10];
% impulseLoc{16} = [2,10];
% impulseLoc{17} = [-2,10];
% impulseLoc{18} = [-4,10];
% impulseLoc{19} = [-4,8];
% impulseLoc{20} = [-4,6];
% impulseLoc{21} = [-4,4];
% impulseLoc{22} = [-4,2];
% 
% figure;
% for i = 1 : 22
%     scatter(i,mean(dS{i,1}),'b');hold on;
%     gt_dist = sqrt((impulseLoc{i}(1)-sensorLoc{1}(1))^2+(impulseLoc{i}(2)-sensorLoc{1}(2))^2);
%     scatter(i,gt_dist/3,'r');hold on;
% end
% hold off;
% 
% figure;
% for i = 1 : 22
%     scatter(i,mean(dS{i,2}),'b');hold on;
%     gt_dist = sqrt((impulseLoc{i}(1)-sensorLoc{2}(1))^2+(impulseLoc{i}(2)-sensorLoc{2}(2))^2);
%     scatter(i,gt_dist/3,'r');hold on;
% end
% hold off;
% 
% figure;
% for i = 1 : 22
%     scatter(i,mean(dS{i,3}),'b');hold on;
%     gt_dist = sqrt((impulseLoc{i}(1)-sensorLoc{3}(1))^2+(impulseLoc{i}(2)-sensorLoc{3}(2))^2);
%     scatter(i,gt_dist/3,'r');hold on;
% end
% hold off;
% 
% figure;
% for i = 1 : 22
%     scatter(i,mean(dS{i,4}),'b');hold on;
%     gt_dist = sqrt((impulseLoc{i}(1)-sensorLoc{4}(1))^2+(impulseLoc{i}(2)-sensorLoc{4}(2))^2);
%     scatter(i,gt_dist/3,'r');hold on;
% end
% hold off;

% 
% save('alignedSig.mat','sI','sS','locsI','locsS','sensorLoc','impulseLoc','noiseIdx','noiseSig','dI','dS','energyRatio');
% 
% % save('alignedSig.mat','sI','sS','locsI','locsS','sensorLoc','impulseLoc','noiseIdx','noiseSig','energyRatio');
% 
% [mV, mI] = min(energyRatio');
% for j = 1 : 18
%     figure;
%     sensorSet = [1:4];
%     sensorSet(sensorSet == mI(j)) = [];
%     peakTime{j} = [-1,-1,-1,-1];
%     sensorCombo{j} = [];
%     for i = sensorSet%1 : 4
%         sensorCombo{j} = [sensorCombo{j}; sensorLoc{i}];
%         plot(locsI{j,i},locsS{j,i});hold on;
%         if max(abs(locsS{j,i})) > 450
%             firstPeak = findFirstPeak(locsS{j,i}, max(abs(locsS{j,i}))/4);
%     %     elseif max(abs(locsS{j,i})) < 100
%     %         firstPeak = findFirstPeak(locsS{j,i}, max(abs(locsS{j,i}))/10);
%         else 
%             firstPeak = findFirstPeak(locsS{j,i}, max(abs(locsS{j,i}))/10);
%         end
%         peakTime{j}(i) = locsI{j,i}(firstPeak);
%         scatter(locsI{j,i}(firstPeak),locsS{j,i}(firstPeak));
%     end
%     hold off;
% end
% 
% %% pick top strongest three to calculate
% errorRecord = [];
% for speed = 1060
%     figure;
%     tempError = 0;
%     locCount = 0;
%     for i = 1 : 18
%         [locX, locY] = simpleMultilateration( peakTime{i}, sensorCombo{i}, speed );
%         if locX == -1 && locY == -1
%             continue;
%         end
%         scatter(locX, locY);hold on;
%         plot([impulseLoc{i}(1),locX],[impulseLoc{i}(2),locY]);hold on;
%         tempError = tempError + sqrt((impulseLoc{i}(1)-locX)^2+(impulseLoc{i}(2)-locY)^2);
%         locCount = locCount + 1;
%     end
%     hold off;
%     errorRecord = [errorRecord; speed, tempError, locCount];
% end
% 
% 
% 
% figure;
% for i = 1 : 18
%     scatter(i,mean(dS{i,1}));hold on;
% end
% hold off;
% 
% figure;
% for i = 1 : 18
%     scatter(i,mean(dS{i,2}));hold on;
% end
% hold off;
% 
% figure;
% for i = 1 : 18
%     scatter(i,mean(dS{i,3}));hold on;
% end
% hold off;
% 
% figure;
% for i = 1 : 18
%     scatter(i,mean(dS{i,4}));hold on;
% end
% hold off;
% 
% save('alignedSig.mat','sI','sS','locsI','locsS','sensorLoc','impulseLoc','noiseIdx','noiseSig','dI','dS');

