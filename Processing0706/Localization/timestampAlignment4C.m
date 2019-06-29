function [ outputSig ] = timestampAlignment4C( signal, timestamps, timestampIdx )
%TIMESTAMPALIGNMENT Summary of this function goes here
%   Detailed explanation goes here
        count = 0;
        aveDelta = [];
        outputSig{1} = [];
        outputSig{2} = [];
        outputSig{3} = [];
        outputSig{4} = [];
        outputIdx = [];
        
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
        for j = 1 : length(blackList)
            for i = 1 : 4
                blackIdx = find(round(timestamps{i}./10^16) == unifiedTime(blackList(j)));
                signal{i}(timestampIdx{i}(blackIdx)) = [];
                timestamps{i}(blackIdx) = [];
                timestampIdx{i}(blackIdx) = [];
            end
        end
        
        for i = 1 : length(timestampIdx{1})-1
            % find the smallest length
            for j = 1 : 4
                tempSig{j} = signal{j}(timestampIdx{j}(i)+1:timestampIdx{j}(i+1)-1);
            end
            targetLen = min([length(tempSig{1}),length(tempSig{2}),length(tempSig{3}),length(tempSig{4})]);
            for j = 1 : 4
                if length(tempSig{j}) > targetLen
                    tempSig{j} = myDownSample(tempSig{j}, length(tempSig{j})-targetLen);
                end
                outputSig{j} = [outputSig{j}; tempSig{j}];
            end
            
        end
        
        figure;
        plot(outputSig{1}(outputSig{1}<10000));hold on;
        plot(outputSig{2}(outputSig{2}<10000));hold on;
        plot(outputSig{3}(outputSig{3}<10000));hold on;
        plot(outputSig{4}(outputSig{4}<10000));hold off;
        
        
end

