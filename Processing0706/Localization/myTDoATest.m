clear all
close all
clc

load('alignedSig.mat');

for j = 8:11
    figure;
    peakTime{j} = [-1,-1,-1,-1];
    for i = 2:3
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
        elseif j == 9 && i == 2
            firstPeak = firstPeak + 6;
        elseif j == 9 && i == 3
            firstPeak = firstPeak + 4;
        elseif j == 10 && i == 2
            firstPeak = firstPeak + 3;
        elseif j == 10 && i == 3
            firstPeak = firstPeak + 5;
        elseif j == 11 && i == 3
            firstPeak = firstPeak + 3;
        end
        peakTime{j}(i) = locsI{j,i}(firstPeak);
        scatter(locsI{j,i}(firstPeak),locsS{j,i}(firstPeak));
    end
    hold off;
end
tdoa = [];
tdoa = [tdoa, peakTime{8}(3)- peakTime{8}(2)];
tdoa = [tdoa, peakTime{9}(3)- peakTime{9}(2)];
tdoa = [tdoa, peakTime{10}(3)- peakTime{10}(2)];
tdoa = [tdoa, peakTime{11}(3)- peakTime{11}(2)];