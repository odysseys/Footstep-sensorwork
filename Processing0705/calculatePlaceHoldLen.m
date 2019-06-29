function [ placeHolderLen ] = calculatePlaceHoldLen( deltaTime, deltaTimeIdx, deltaVal, errorIdx )
%CALCULATEPLACEHOLDLEN Summary of this function goes here
%   Detailed explanation goes here
    placeHolderLen = [];
    for i = 1 : length(errorIdx)
        rightLen = deltaTime(errorIdx(i))/deltaVal/10^11;
        missingLen = rightLen - deltaTimeIdx(errorIdx(i));
        placeHolderLen = [placeHolderLen, missingLen];
    end
end

