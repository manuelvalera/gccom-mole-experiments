function [uatv] = InterpVtoV(Q,m,n,o,type)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here  
   
uatv = [];

switch type  
    
case 'utov'
        uatv = zeros(m,n+1,o);
   
        for i=1:m-1
        for j=2:n
        for k=1:o-1
            uatv(i,j,k) = (Q(i,j,k)+Q(i,j-1,k)+Q(i+1,j,k)+Q(i+1,j-1,k))*0.250; %
        end
        end
        end

        
case 'utow'
        
        uatv = zeros(m,n,o+1);
        
        for i=1:m-1
        for j=1:n-1
        for k=2:o
            uatv(i,j,k) = (Q(i,j,k)+Q(i,j,k-1)+Q(i+1,j,k)+Q(i+1,j,k-1))*0.250; %
        end
        end
        end
        
case 'vtou'
    
        uatv = zeros(m+1,n,o);
        
        for i=2:m
        for j=1:n-1
        for k=1:o-1
            uatv(i,j,k) = (Q(i,j,k)+Q(i,j+1,k)+Q(i-1,j+1,k)+Q(i-1,j,k))*0.250; %
        end
        end
        end
        
case 'vtow'
    
        uatv = zeros(m,n,o+1);
        
        for i=1:m-1
        for j=1:n-1
        for k=2:o
            uatv(i,j,k) = (Q(i,j,k)+Q(i,j,k-1)+Q(i,j+1,k)+Q(i,j+1,k-1))*0.250; %
        end
        end
        end
        
case 'wtou'
        
        uatv = zeros(m+1,n,o);
    
        for i=2:m
        for j=1:n-1
        for k=1:o-1
            uatv(i,j,k) = (Q(i,j,k)+Q(i,j,k+1)+Q(i-1,j,k+1)+Q(i-1,j,k))*0.250; %
        end
        end
        end
        
case 'wtov'
    
        uatv = zeros(m,n+1,o);
        
        for i=1:m-1
        for j=2:n
        for k=1:o-1
            uatv(i,j,k) = (Q(i,j,k)+Q(i,j,k+1)+Q(i,j-1,k+1)+Q(i,j-1,k))*0.250; %
        end
        end
        end
        
end

   

end

