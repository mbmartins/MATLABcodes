function [outputArg1] = phase(z)
%PHASE Summary of this function goes here
%   Detailed explanation goes here
outputArg1 = unwrap(angle(z));
end

