function [adTdx,bdTdy] = MimeticUpwindD2D(Tu,Tv,coeffa,coeffb,Guf,Gub,Gvf,Gvb,m,n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

u_ct = (m+1)*n;


%Du:

Solplus  = Guf*Tu; %forward
Solminus = Gub*Tu; %backwards

aplus = [];
aminus = [];

for i=1:size(coeffa)

        aplus(i)  = max(coeffa(i),0);
        aminus(i) = min(coeffa(i),0);
              
end

adTdx = aplus.*Solminus' + aminus.*Solplus';  

%Upwind wrt x

adTdx = adTdx';

%Dv:
Solvec = [];
Solplus = [];
Solminus = [];

Solplus  = Gvf*Tv; %forward
Solminus = Gvb*Tv; %backwards

aplus = [];
aminus = [];

for i=1:size(coeffb)
        aplus(i)  = max(coeffb(i),0);
        aminus(i) = min(coeffb(i),0);
              
end

bdTdy = aplus.*Solminus' + aminus.*Solplus'; 

bdTdy = bdTdy';






