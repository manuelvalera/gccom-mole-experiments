function [adTdx,bdTdy] = MimeticUpwindv3grad(T,coeffa,coeffb,I0l,I1l,I0r,I1r,G,m,n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

u_ct = (m+1)*n;

div = G*[T]; %T is ls-size

%Upwind wrt x

divx = div(1:u_ct);

Solplus = I0l*divx;
Solminus = I1l*divx;

aplus = [];
aminus = [];

for i=1:size(coeffa)  
        aplus(i)  = max(coeffa(i),0);
        aminus(i) = min(coeffa(i),0);
end

Solvec = aplus.*Solminus' + aminus.*Solplus';  
adTdx = Solvec';

Solvec = [];
Solplus = [];
Solminus = [];

%Upwind wrt y

divy = div(u_ct+1:end);

Solplus = I0r*divy;
Solminus = I1r*divy;

aplus = [];
aminus = [];

for i=1:size(coeffb)
    
        aplus(i)  = max(coeffb(i),0);
        aminus(i) = min(coeffb(i),0);
   
end

Solvec = aplus.*Solminus' + aminus.*Solplus';  
bdTdy = Solvec';

