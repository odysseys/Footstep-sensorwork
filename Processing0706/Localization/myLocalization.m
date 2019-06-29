clear all
close all
clc

load('alignedSig.mat');

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
        elseif j == 3
            firstPeak = findFirstPeak(locsS{j,i}, max(abs(locsS{j,i}))/7);
        else 
            firstPeak = findFirstPeak(locsS{j,i}, max(abs(locsS{j,i}))/10);
        end
        if j == 17 && i == 4
            firstPeak = firstPeak + 10;
        elseif j == 16 && i == 4
            firstPeak = firstPeak + 7;
        elseif j == 16 && i == 1
            firstPeak = firstPeak + 4;
        elseif j == 13 && i == 2
            firstPeak = firstPeak + 4;
        elseif j == 13 && i == 3
            firstPeak = firstPeak + 4;
        elseif j == 3 && i == 2
            firstPeak = firstPeak - 75;
        end
        peakTime{j}(i) = locsI{j,i}(firstPeak);
        scatter(locsI{j,i}(firstPeak),locsS{j,i}(firstPeak));
    end
    hold off;
end

%% pick top strongest three to calculate
errorRecord = [];
for speed = 1430
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

aveErrorRate = errorRecord(:,2)./errorRecord(:,3)./3

figure;
for i = 1 : 18
    scatter(i,mean(dS{i,1}-1.5),'b');hold on;
    gt_dist = sqrt((impulseLoc{i}(1)-sensorLoc{1}(1))^2+(impulseLoc{i}(2)-sensorLoc{1}(2))^2);
    scatter(i,gt_dist/3,'r');hold on;
end
hold off;

figure;
for i = 1 : 18
    scatter(i,mean(dS{i,2}-1.5),'b');hold on;
    gt_dist = sqrt((impulseLoc{i}(1)-sensorLoc{2}(1))^2+(impulseLoc{i}(2)-sensorLoc{2}(2))^2);
    scatter(i,gt_dist/3,'r');hold on;
end
hold off;

figure;
for i = 1 : 18
    scatter(i,mean(dS{i,3}-1.5),'b');hold on;
    gt_dist = sqrt((impulseLoc{i}(1)-sensorLoc{3}(1))^2+(impulseLoc{i}(2)-sensorLoc{3}(2))^2);
    scatter(i,gt_dist/3,'r');hold on;
end
hold off;

figure;
for i = 1 : 18
    scatter(i,mean(dS{i,4}-1.5),'b');hold on;
    gt_dist = sqrt((impulseLoc{i}(1)-sensorLoc{4}(1))^2+(impulseLoc{i}(2)-sensorLoc{4}(2))^2);
    scatter(i,gt_dist/3,'r');hold on;
end
hold off;

%% using radio time of flight to localize



% figure;
% for i = 1 : 4
% plot(locsI{1,i},locsS{1,i});hold on;
% firstPeak = findFirstPeak(abs(locsS{1,i}), max(abs(locsS{1,i}))/7);
% scatter(locsI{1,i}(firstPeak),locsS{1,i}(firstPeak));
% end
% hold off;
% 
% figure;
% for i = 1 : 18
%     scatter(i,convert2D(mean(dS{i,1}),1));hold on;
% end
% hold off;
% 
% figure;
% for i = 1 : 18
%     scatter(i,convert2D(mean(dS{i,2}),1));hold on;
% end
% hold off;
% 
% figure;
% for i = 1 : 18
%     scatter(i,convert2D(mean(dS{i,3}),1));hold on;
% end
% hold off;
% 
% figure;
% for i = 1 : 18
%     scatter(i,convert2D(mean(dS{i,4}),1));hold on;
% end
% hold off;