function [ d ] = convert2D( dist, height )
%CONVERT2D Summary of this function goes here
%   Detailed explanation goes here

    d = sqrt(dist*dist-height*height);
end

