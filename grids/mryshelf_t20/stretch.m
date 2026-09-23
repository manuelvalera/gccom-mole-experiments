function xp = stretch(s, L, beta)
    D = L/2;
    A = (1/(2*beta))*log((1+(exp(beta)-1)*D/L)/(1+(exp(-beta)-1)*D/L));
    xp = D*(1 + sinh(beta*(s-A))/sinh(beta*A));
end
