function [p] = SSPRK102D(p,RHS,dt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Qt1 = p;
Qt2 = p;
%Qtold = Qt;

for i=1:5
    Qt1(:,:) = Qt1(:,:) + dt*RHS(:,:)/6.0;
end

Qt2(:,:) = 1.0/25.0*Qt2(:,:) + 9.0/25.0*Qt1(:,:);
Qt1(:,:) = 15.0*Qt2(:,:) - 5.0*Qt1(:,:);

for i=6:9
    Qt1(:,:) = Qt1(:,:) + dt*RHS(:,:)/6.0;
end

Qt1(:,:) = Qt2(:,:) + 3.0/5.0*Qt1(:,:) + 1.0/10.0*dt*RHS(:,:);

p = Qt1;

end

