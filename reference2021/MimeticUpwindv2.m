function [divv] = MimeticUpwindv2(T,coeff,I0,I1,G)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Solplus = I0*T;
Solminus = I1*T;



for i=1:size(coeff)
    
        aplus  = max(coeff(i),0);
        aminus = min(coeff(i),0);
        
        Solvec(i) = aplus*Solminus(i) + aminus*Solplus(i);  
   
end

divv = G*[Solvec'];



