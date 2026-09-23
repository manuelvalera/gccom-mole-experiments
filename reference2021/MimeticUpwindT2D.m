function [RHST] = MimeticUpwindT2D(uatu,vatv,Tvec,DcfA,DcbA,m,n)

    uatl = zeros(n+2,m+2);
    vatl = zeros(n+2,m+2);

    

     for i=1:m
     for j=1:n
         uatl(j,i) = (uatu(j,i) + uatu(j,i+1))/2;
         vatl(j,i) = (vatv(j,i) + vatv(j+1,i))/2;
     end
     end
%   
    uatl(2,:) = uatl(3,:);
    uatl(1,:) = uatl(2,:);
    uatl(end-1,:) = uatl(end-2,:);
    uatl(end,:) = uatl(end-1,:);
     
    uatl(:,2) = uatl(:,3);
    uatl(:,1) = uatl(:,2);
    uatl(:,end-1) = uatl(:,end-2);
    uatl(:,end) = uatl(:,end-1);

    vatl(2,:) = vatl(3,:);
    vatl(1,:) = vatl(2,:);
    vatl(end-1,:) = vatl(end-2,:);
    vatl(end,:) = vatl(end-1,:);
     
    vatl(:,2) = vatl(:,3);
    vatl(:,1) = vatl(:,2);
    vatl(:,end-1) = vatl(:,end-2);
    vatl(:,end) = vatl(:,end-1);

    Ul = reshape(uatl.',[], 1);
    Vl = reshape(vatl.',[], 1);
    
    Solplus  = DcfA*Tvec; %forward
    Solminus = DcbA*Tvec; %backwards

    aplus = [];
    aminus = [];

    for i=1:size(Ul)

            aplus(i)  = max(Ul(i),0);
            aminus(i) = min(Ul(i),0);

    end

    adTdx = aplus.*Solminus' + aminus.*Solplus';  

    %Upwind wrt x

    adTdx = adTdx';

    %Dv:
    Solvec = [];
    Solplus = [];
    Solminus = [];

    Solplus  = DcfA*Tvec; %forward
    Solminus = DcbA*Tvec; %backwards

    aplus = [];
    aminus = [];

    for i=1:size(Vl)
            aplus(i)  = max(Vl(i),0);
            aminus(i) = min(Vl(i),0);

    end

    bdTdy = aplus.*Solminus' + aminus.*Solplus'; 

    RHST = -(adTdx + bdTdy');


end

