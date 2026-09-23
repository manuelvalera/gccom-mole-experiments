function [adTdx,bdTdy,cdTdz] = MimeticUpwindD3D(Tu,Tv,Tw,coeffa,coeffb,coeffc,Guf,Gub,Gvf,Gvb,Gwf,Gwb,m,n,o)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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

%Dz:
Solvec = [];
Solplus = [];
Solminus = [];

Solplus  = Gwf*Tw; %forward
Solminus = Gwb*Tw; %backwards

aplus = [];
aminus = [];

for i=1:size(coeffc)
        aplus(i)  = max(coeffc(i),0);
        aminus(i) = min(coeffc(i),0);
              
end

cdTdz = aplus.*Solminus' + aminus.*Solplus'; 

%[m,ind] = min(cdTdz)
%[m,ind] = max(cdTdz)


cdTdz = cdTdz';






