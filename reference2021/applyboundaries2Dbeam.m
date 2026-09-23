function [u,v,p] = applyboundaries2Dbeam(u,v,p)

%   

%% left-right
      
     u(:,1) = u(:,2);
     u(:,end) = u(:,end-1);
     v(:,1) = v(:,2);
     v(:,end) = v(:,end-1);  

%     u(:,1) = 0.0 ; %u(:,2);
%     u(:,end) = 0.0 ; %u(:,JMax);    
%     v(:,1) = 0.0;
%     v(:,end) = 0.0;  

     

     
     
%% top-bottom    
    
   
    u(1,:) = u(2,:);
    u(end,:) = u(end-1,:);          
%     v(1,:) = v(2,:);    
%     v(end,:) = v(end-1,:);    

%    u(1,:) = 0.0;
%    u(end,:) = 0.0 ; %u(IMax,:);      
   v(1,:) = 0.0;    
    v(end,:) = 0.0;    %must be 0
    
    
%% Pressure    
    
     p(1,:) = 0.0 ; % p(2,:);
     p(end,:) = 0.0; % p(IMax,:);    
     p(:,1) = 0.0; % p(:,2);
     p(:,end) = 0.0 ; %p(:,JMax);   
%     
%   p(2,:) = p(3,:);
%   p(end-1,:) = p(end-2,:);  
%   p(:,2) = p(:,3); 
%   p(:,end-1) = p(:,end-2);   

%   p(1,:) = p(2,:);
%   p(end,:) = p(end-1,:);    
%   p(:,1) = p(:,2);
%   p(:,end) = p(:,end-1);   
    
    
end
