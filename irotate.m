function [ output_args ] = irotate( input_args )
%IROTATE Summary of this function goes here
%   Detailed explanation goes here
T=[ 0.9483    0.3169    0.0178;   -0.3135    0.9262    0.2099;    0.0499   -0.2046    0.9775];
FaiE = 180/pi*atan((T(1,3)-T(3,1))/(T(3,2)-T(2,3)))
LanmdaE = 180/pi*asin((T(2,1)-T(1,2))/sqrt((T(3,2)-T(2,3))^2+ (T(1,3)-T(3,1))^2+(T(2,1)-T(1,2))^2))
Omiga = 180/pi*atan(sqrt((T(3,2)-T(2,3))^2+ (T(1,3)-T(3,1))^2+(T(2,1)-T(1,2))^2)/(T(1,1)+T(2,2)+T(3,3)-1))
