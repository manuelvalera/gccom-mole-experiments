function [T] = expand2D(T)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


    T(2,:) = T(3,:);
    T(:,2) = T(:,3);
    
    T(end-1,:) = T(end-2,:);
    T(:,end-1) = T(:,end-2);   

    T(1,:) = T(2,:);
    T(:,1) = T(:,2);

    T(end,:) = T(end-1,:);
    T(:,end) = T(:,end-1);
end

