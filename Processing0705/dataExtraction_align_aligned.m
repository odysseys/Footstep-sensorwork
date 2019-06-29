clear all
close all
clc

load('alignedSig.mat');

beginTimestamp = min([s1I(1),s2I(1),s3I(1),s4I(1)]);
s1ICompare = abs(s1I - beginTimestamp);
s2ICompare = abs(s2I - beginTimestamp);
s3ICompare = abs(s3I - beginTimestamp);
s4ICompare = abs(s4I - beginTimestamp);

