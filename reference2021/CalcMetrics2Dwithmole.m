function [jac,i_x,i_y,j_x,j_y] = CalcMetrics2Dwithmole(x_c,y_c,IMax,JMax)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Grid Derivatives:

[jac, x_i, x_j, y_i, y_j] = jacobian2D(2, x_c, y_c);

jac = reshape(jac,[IMax JMax]).' ;  
jac = jac';
x_i = reshape(x_i,[IMax JMax]).' ;    
x_j = reshape(x_j,[IMax JMax]).' ;    
y_i = reshape(y_i,[IMax JMax]).' ;    
y_j = reshape(y_j,[IMax JMax]).' ;    

x_i = expand2D(x_i'); x_j = expand2D(x_j');
y_i = expand2D(y_i'); y_j = expand2D(y_j');

%% Curvilining:

Idi = 1/(IMax-2); Idj = 1/(JMax-2); 

%Idi = 1.0; Idj= 1.0;

for j=1:JMax
    for i=1:IMax

        % x :
        x_i(i,j) = x_i(i,j)*Idi; %0
        x_j(i,j) = x_j(i,j)*Idi; %1

        % y :
        y_i(i,j) = y_i(i,j)*Idj;%3
        y_j(i,j) = y_j(i,j)*Idj;%4

    end
end

%% Inverse metrics:


i_x = zeros(IMax,JMax); i_y = zeros(IMax,JMax); 
j_x = zeros(IMax,JMax); j_y = zeros(IMax,JMax); 

for j=1:JMax
    for i=1:IMax

        i_x(i,j) =  jac(i,j)*( y_j(i,j)) ;
        i_y(i,j) = -jac(i,j)*( x_j(i,j)) ;

        j_x(i,j) = -jac(i,j)*( y_i(i,j)) ;
        j_y(i,j) =  jac(i,j)*( x_i(i,j)) ;
     
    end
end

i_x = expand2D(i_x); j_x = expand2D(j_x);
i_y = expand2D(i_y); j_y = expand2D(j_y);


