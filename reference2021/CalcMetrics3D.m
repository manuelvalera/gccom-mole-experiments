function [i_x,i_y,i_z,j_x,j_y,j_z,k_x,k_y,k_z] = CalcMetrics3D(x_c,y_c,z_c,IMax,JMax,KMax)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Grid Derivatives:

[x_i,x_j,x_k] =   CalcDeriv1(x_c,IMax,JMax,KMax);

[y_i,y_j,y_k] =   CalcDeriv1(y_c,IMax,JMax,KMax);

[z_i,z_j,z_k] =   CalcDeriv1(z_c,IMax,JMax,KMax);

%% Curvilining:

Idi = 1/(IMax-2); Idj = 1/(JMax-2); Idk = 1/(KMax-2);

for k=1:KMax
    for j=1:JMax
        for i=1:IMax
            
            % x :
            x_i(i,j,k) = x_i(i,j,k)*Idi; %0
            x_j(i,j,k) = x_j(i,j,k)*Idi; %1
            x_k(i,j,k) = x_k(i,j,k)*Idi; %2
            
            % y :
            y_i(i,j,k) = y_i(i,j,k)*Idj;%3
            y_j(i,j,k) = y_j(i,j,k)*Idj;%4
            y_k(i,j,k) = y_k(i,j,k)*Idj;%5
            
            % z :
            z_i(i,j,k) = z_i(i,j,k)*Idk; %6
            z_j(i,j,k) = z_j(i,j,k)*Idk; %7
            z_k(i,j,k) = z_k(i,j,k)*Idk; %8
            
        end
    end
end
%% Jacobian:
jac = zeros(IMax,JMax,KMax);

for k=1:KMax
    for j=1:JMax
        for i=1:IMax
            
            jac(i,j,k) = 1.0/( x_i(i,j,k)*(y_j(i,j,k)*z_k(i,j,k)- y_k(i,j,k)*z_j(i,j,k))...
                - x_j(i,j,k)*(y_i(i,j,k)*z_k(i,j,k)-y_k(i,j,k)*z_i(i,j,k))...
                + x_k(i,j,k)*(y_i(i,j,k)*z_j(i,j,k)-y_j(i,j,k)*z_i(i,j,k)));
            
        end
    end
end

%% Inverse metrics:


i_x = zeros(IMax,JMax,KMax); i_y = zeros(IMax,JMax,KMax); i_z = zeros(IMax,JMax,KMax);
j_x = zeros(IMax,JMax,KMax); j_y = zeros(IMax,JMax,KMax); j_z = zeros(IMax,JMax,KMax);
k_x = zeros(IMax,JMax,KMax); k_y = zeros(IMax,JMax,KMax); k_z = zeros(IMax,JMax,KMax);



for k=1:KMax
    for j=1:JMax
        for i=1:IMax
            
            i_x(i,j,k) = jac(i,j,k)*( y_j(i,j,k)*z_k(i,j,k) - y_k(i,j,k)*z_j(i,j,k)) ;
            i_y(i,j,k) = jac(i,j,k)*( x_k(i,j,k)*z_j(i,j,k) - x_j(i,j,k)*z_k(i,j,k)) ;
            i_z(i,j,k) = jac(i,j,k)*( x_j(i,j,k)*y_k(i,j,k) - x_k(i,j,k)*y_j(i,j,k)) ;
            
            j_x(i,j,k) = jac(i,j,k)*( y_k(i,j,k)*z_i(i,j,k) - y_i(i,j,k)*z_k(i,j,k)) ;
            j_y(i,j,k) = jac(i,j,k)*( x_i(i,j,k)*z_k(i,j,k) - x_k(i,j,k)*z_i(i,j,k)) ;
            j_z(i,j,k) = jac(i,j,k)*( x_k(i,j,k)*y_i(i,j,k) - x_i(i,j,k)*y_k(i,j,k)) ;
            
            k_x(i,j,k) = jac(i,j,k)*( y_i(i,j,k)*z_j(i,j,k) - y_j(i,j,k)*z_i(i,j,k)) ;
            k_y(i,j,k) = jac(i,j,k)*( x_j(i,j,k)*z_i(i,j,k) - x_i(i,j,k)*z_j(i,j,k)) ;
            k_z(i,j,k) = jac(i,j,k)*( x_i(i,j,k)*y_j(i,j,k) - x_j(i,j,k)*y_i(i,j,k)) ;
            
        end
    end
end


end

