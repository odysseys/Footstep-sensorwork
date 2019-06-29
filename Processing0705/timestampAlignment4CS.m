function [ outputSig, outputIdx ] = timestampAlignment4C( signal, timestamps, timestampIdx )
%TIMESTAMPALIGNMENT Summary of this function goes here
%   Detailed explanation goes here
        count = 0;
        aveDelta = [];
        outputSig{1} = [];
        outputSig{2} = [];
        outputSig{3} = [];
        outputSig{4} = [];
        outputIdx{1} = [];
        outputIdx{2} = [];
        outputIdx{3} = [];
        outputIdx{4} = [];
        
        unifiedTime = round([timestamps{1}; timestamps{2}; timestamps{3}; timestamps{4}]./10^16);
        unifiedTime = unique(unifiedTime);
        blackList = [];
        for i = 1 : length(unifiedTime)
            if ~ismember(unifiedTime(i),round(timestamps{1}./10^16)) || ...
                    ~ismember(unifiedTime(i),round(timestamps{2}./10^16)) || ...
                    ~ismember(unifiedTime(i),round(timestamps{3}./10^16)) || ...
                    ~ismember(unifiedTime(i),round(timestamps{4}./10^16)) 
                    blackList = [blackList, i];
            end
        end
        
        for i = 1 : 4
            BL{i} = [];
            BLT{i} = [];
        end
        
        for j = 1 : length(blackList)
            for i = 1 : 4
                blackIdx = find(round(timestamps{i}./10^16) == unifiedTime(blackList(j)));
                BL{i} = [BL{i}, timestampIdx{i}(blackIdx)];
                BLT{i} = [BLT{i}, blackIdx];
            end
        end
        
        for i = 1 : 4
            signal{i}(BL{i}) = [];
            timestamps{i} = signal{i}(signal{i}>10000);
            timestampIdx{i} = find(signal{i}>10000);
        end
        
                
        
        for i = 1 : length(timestampIdx{1})-1
            % find the smallest length
            for j = 1 : 4
                tempSig{j} = signal{j}(timestampIdx{j}(i)+1:timestampIdx{j}(i+1)-1);
                tempDeltaTime{j} = timestamps{j}(i+1)-timestamps{j}(i);
                tempSpeed{j} = tempDeltaTime{j}/length(tempSig{j});
            end
%             targetLen = min([length(tempSig{1}),length(tempSig{2}),length(tempSig{3}),length(tempSig{4})]);
            targetSpeed = max([tempSpeed{1},tempSpeed{2},tempSpeed{3},tempSpeed{4}])
            for j = 1 : 4
                targetLen = round(tempDeltaTime{j}/targetSpeed);
                if length(tempSig{j}) > targetLen
                    tempSig{j} = myDownSample(tempSig{j}, length(tempSig{j})-targetLen);
                end
                deltaTime = (timestamps{j}(i+1)-timestamps{j}(i))/length(tempSig{j}+1);
                tempIdx{j} = zeros(length(tempSig{j}),1);
                tempIdx{j}(1) = timestamps{j}(i);
                for k = 2 : length(tempSig{j})
                    tempIdx{j}(k) = tempIdx{j}(k-1)+deltaTime;
                end
                outputSig{j} = [outputSig{j}; tempSig{j}];
                outputIdx{j} = [outputIdx{j}; tempIdx{j}];
            end
            
        end
        
        figure;
        plot(outputSig{1}(outputSig{1}<10000));hold on;
        plot(outputSig{2}(outputSig{2}<10000));hold on;
        plot(outputSig{3}(outputSig{3}<10000));hold on;
        plot(outputSig{4}(outputSig{4}<10000));hold off;
        
        figure;
        plot(outputIdx{1},outputSig{1});hold on;
        plot(outputIdx{2},outputSig{2});hold on;
        plot(outputIdx{3},outputSig{3});hold on;
        plot(outputIdx{4},outputSig{4});hold off;
        ylim([0,1000]);
        
end

