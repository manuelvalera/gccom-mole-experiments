function [Solvec] = MimeticUpwind(T,coeff,I0,I1)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Solplus = I0*T;
Solminus = I1*T;

csign = sign(coeff);

for i=1:size(coeff)
    if csign(i) <= 0
        Solvec(i) = Solplus(i);
    elseif csign(i) > 0
        Solvec(i) = Solminus(i);
    end
end

