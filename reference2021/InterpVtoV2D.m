function [uatv] = InterpVtoV2D(Q,m,n,type)

uatv = [];

switch type  
    
    case 'utov'
        uatv = zeros(m,n+1);
   
        for i=1:m
        for j=2:n
        
            uatv(i,j) = (Q(i,j)+Q(i,j-1)+Q(i+1,j)+Q(i+1,j-1))*0.250; %
        end
        end
        
 
         %uatv(:,1) = uatv(:,2);
         
    case 'vtou'
    
        uatv = zeros(m+1,n);
        
        for i=2:m
        for j=1:n
        
            uatv(i,j) = (Q(i,j)+Q(i,j+1)+Q(i-1,j+1)+Q(i-1,j))*0.250; %
        end
        end
        
        
        %uatv(1,:) = uatv(2,:);
        
end




end
