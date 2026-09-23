close all; clear all; clc

%IMax=130; JMax = 100;
IMax=200; JMax = 100;


m = IMax-1;
n = JMax-1; 

N=0.007;   %Brunt-Vaisala frequency
%w_N=0.8    %w/N parameter
w_N=0.6    %w/N parameter
%w_N=0.4    %w/N parameter
%w_N=0.2    %w/N parameter
omega=w_N*N %forcing frequency

T20 = 2*pi()/omega

w_N2=(w_N)^2;
Ho = 1000; 
ab = 20; 

phy=atan(sqrt(w_N2/(1-w_N2)));

L = 3000; %4*Ho/tan(phy);
%Lb= 2*1.9;
Lb=L/50;
beta = 12;
ab= -Lb*tan(beta); %in order to make the angle 15
D=Ho;

LStar = 1.0;
UStar = 1.0;

a=-L/2; b=L/2;
c=-Ho ; d=0; 

Dx = (b-a)/m ;
Dy = (d-c)/n ;

%[Y, X] = meshgrid(c:Dy:d, a:Dx:b);
%

xlength = linspace(a,b,m+1);
ylength = linspace(c,d,n+1);
%
%[Y,X] = meshgrid(ylength,xlength);
[X,Y] = meshgrid(xlength,ylength);

for i=1:IMax
    for j=1:JMax
  
            Y(j,i) = abs((JMax-j)/(JMax-1))*(-( Ho - ab*exp( -( X(1,i)^2 )/ (2*(Lb^2))) ));
            
    end
end

% plot(squeeze(Y(:,1)))
% 
% figure(2)
% for j=1:JMax
% plot(squeeze(Y(j,:)),'-*')
% hold on
% end
% hold off
%
%

Xu = (X(1:end-1, :) + X(2:end, :))/2;
Yu = (Y(1:end-1, :) + Y(2:end, :))/2;

Xv = (X(:, 1:end-1) + X(:, 2:end))/2;
Yv = (Y(:, 1:end-1) + Y(:, 2:end))/2;

%
xgrid = [a a+Dx/2 : Dx : b-Dx/2 b];
ygrid = [c c+Dy/2 : Dy : d-Dy/2 d];

[Xls, Yls] = meshgrid(xgrid,ygrid);
%
 
for i=1:IMax+1
    for j=1:JMax+1
  
            Yls(j,i) = abs((JMax-j)/(JMax-1))*(-( Ho - ab*exp( -( Xls(1,i)^2 )/ (2*(Lb^2))) ));
            
    end
end
%
% figure(2)
% for j=1:JMax
% plot(squeeze(Yls(j,:)),'-*')
% hold on
% end
% hold off
% 
% figure(1)
% pcolor(Yls)
% shading interp

%
% clc
% [xu_i,xu_j] = CalcDeriv12D(Xu',m+1,n);
% [yu_i,yu_j] = CalcDeriv12D(Yu',m+1,n);
% [xv_i,xv_j] = CalcDeriv12D(Xv',m,n+1);
% [yv_i,yv_j] = CalcDeriv12D(Yv',m,n+1);

[Ju,i_xu,i_yu,j_xu,j_yu] = CalcMetrics2Dwithmole(Xu,Yu,m+1,n);
[Jv,i_xv,i_yv,j_xv,j_yv] = CalcMetrics2Dwithmole(Xv,Yv,m,n+1);
[Jl,i_xl,i_yl,j_xl,j_yl] = CalcMetrics2Dwithmole(Xls,Yls,m+2,n+2);


% 
% xu_i = xu_i'; xu_j = xu_j';
% xv_i = xv_i'; xv_j = xv_j';
% 
% yu_i = yu_i'; yu_j = yu_j';
% yv_i = yv_i'; yv_j = yv_j';


%
%clc
% [i_xu,i_yu,j_xu,j_yu] = CalcMetrics2D(Xu',Yu',m+1,n);
% [Jv,i_xv,i_yv,j_xv,j_yv] = CalcMetrics2Dv2(Xv',Yv',m,n+1);
% 
% [i_xl,i_yl,j_xl,j_yl] = CalcMetrics2D(Xls',Yls',m+2,n+2);

Idiu = 1/(m-1); Idju = 1/(n-2); 
Idiv = 1/(m-2); Idjv = 1/(n-1); 
Idil = 1/(m+1); Idjl = 1/(n+1); 

%%

figure(3)
pcolor(Xu(:,:),Yu(:,:),j_xu');
%pcolor(Xv(:,:),Yv(:,:),j_xv');
%pcolor(Xls(:,:),Yls(:,:),j_xl');

colorbar



%%
TMAX = 180; TMIN = 179.9705;
Tic = TMAX + Yls*(TMAX-TMIN);


rho_ref = 1027.0;
S_ref   = 35.0;
T_ref   = 10.0;
drho_dS = 0.7810;
drho_dT = -0.17080;
gforce = 9.8;

dens = rho_ref + drho_dT*(Tic - T_ref);
rho0=mean(dens(:));

Ex = 0.066;

uo = Ex*omega
%uo = 0.1;
%uo=0.01;

%T0 = 10.0;
T = Tic;

layu = [m+1 n];
layv = [m n+1];

% fig = figure(jj);
% u_v = squeeze(dens');
% subplot(2,1,1)
% pcolor(u_v);
% axis tight
% title('batv')
% colorbar;
% shading interp;
% subplot(2,1,2)
% u_3 = squeeze(Tic');
% pcolor(u_3);
% axis tight
% title('batv_sq')
% colorbar;
% shading interp;

buoy = gforce*(dens-rho0)/rho0;
%buoy = gforce*(dens-rho_ref)/rho_ref;

%Fbuoy = griddedInterpolant(Xls', Yls', buoy'); 
%batv = Fbuoy(Xv',Yv')';

for i=2:m
  for j=2:n+1
         batv(j,i) = (buoy(j,i) + buoy(j-1,i))/2;
  end
end

%batv(:,2) = batv(:,3);
% batv(:,1) = batv(:,2);
% batv(end,:) = batv(end-1,:);
% batv(1,:) = batv(2,:);

Bb = batv(:);
%batv_sq = reshape(Bb,[m n+1]).';  %reshape(udvdx,[n+1 m]); 
jj=2


%
fig = figure(jj);
u_v = squeeze(batv');
%subplot(2,1,1)
pcolor(u_v');
axis tight
title('batv')
colorbar;
shading flat;

RHSu  = sparse(n*(m+1),1);
RHSv  = sparse((n+1)*m,1);
T_new = sparse(n+2,m+2);
D2T = sparse(n+2,m+2);
p = T_new;

u = zeros(size(Xls));
v = zeros(size(Xls));

dudxc = zeros(size(Xls));
dvdyc = zeros(size(Xls));


Uu = RHSu; Vu = Uu; 
Vv = RHSv; Uv = Vv; 
%
%u,v derivatives:

Uls = reshape(u',[], 1);
Vls = reshape(v',[], 1);   
%
uatu = zeros(n,m+1);
vatu = zeros(n,m+1);

uatv = zeros(n+1,m);
vatv = zeros(n+1,m);

uatl = zeros(n+2,m+2);
vatl = zeros(n+2,m+2);

% fig = figure(jj);
% %set(gcf,'units','inches','position',[10 10 15 10])
% u_v = squeeze(uatu(:,:)*UStar);
% subplot(2,2,1)
% pcolor(Xu(:,:),Yu(:,:),u_v);
% %pcolor(u_v);
% %axis equal
% axis tight
% title('u [m/s]')
% colorbar;
% shading interp;
% subplot(2,2,2)
% u_3 = squeeze(vatv(:,:)*UStar);
% pcolor(Xv(:,:),Yv(:,:),u_3);
% %pcolor(u_3);
% %axis equal
% axis tight
% title('v [m/s]')
% colorbar;
% shading interp;
% subplot(2,2,3)
% u_2 = squeeze(batv(:,:));
% pcolor(Xv(:,:),Yv(:,:),u_2);
% %pcolor(u_2);
% %axis equal
% axis tight 
% title('buoy')
% colorbar;
% shading interp;
% subplot(2,2,4)
% u_2 = squeeze(T(:,:));
% pcolor(Xls(:,:),Yls(:,:),u_2);
% %pcolor(u_2);
% %axis equal
% axis tight 
% title(['T [C] at ' num2str(current_time) 's'])
% colorbar;
% shading interp;
%

Re = 3500;
%Re = 750/80;
IRe = 1/Re;

earth_omega = 7.2921e-5;

latitude = 24;

f_cor = 2*earth_omega*sin(latitude);
f_cor_h = 2*earth_omega*cos(latitude);
alpha = 1.664e-4 ;

nuatv = 1e-6 ; %according to fringer

% Fuatu = griddedInterpolant(Xu', Yu', uatu'); 
% Fvatv = griddedInterpolant(Xv', Yv', vatv'); 
% 
% u = Fuatu(Xls',Yls')';
% v = Fvatv(Xls',Yls')';
%
k = 2; % Mimetic order of accuracy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dxc = Dx;
Dyc = Dy; 

% D =  div2D(k,m,Dxc,n,Dyc)*LStar; 
% G =  grad2D(k,m,Dxc,n,Dyc)*LStar;
% 
% R = LStar^2*robinBC2D(k,m,Dxc,n,Dyc,0,1);
% N = Neumann2DCurv(G,m,n,1);
% 
% R = R-N;

D =  div2DCurv(2,X,Y)*LStar; 
G =  grad2DCurv(2,X,Y)*LStar;

%N = Neumann2DCurv(G,m,n,1);

%L = lap2D(k, m, Dx, n, Dy);

%L = D*G + Neumann2DCurv(G,m,n,1);  

%L = D*G + LStar^2*robinBC2D(2,m,Dxc,n,Dyc,1,0); %Dirichlet

%Robin = LStar^2*robinBC2D(2,m,Dxc,n,Dyc,1,0); 


Robin = Neumann2DCurv(G,m,n,1); 
bdry = find(diag(Robin));
L = D*G + Robin; %Dirichlet
Gp = -1/rho0*G;

Tbc=2*pi/omega; % Tidal Period

tf = 650;
%tf = 35; 
%tf = Tbc*2; 
%tf = 260.0; 
%tf = 22320.0; %GCCOM tf 
%tf = 19.25*Tbc;
%tf = dt*5;
tf = tf*(UStar/LStar);

dt = 5.0;
%dt = Tbc/500;
time = 0;
current_time = (time+1)*dt*LStar/UStar


%CFL dt limit:
%maxv = UStar;
%dt = (min(Dx,Dy)/maxv)/2;

dt = dt*(UStar/LStar);
nf = tf/dt ;
time_steps = linspace(0,tf,nf);

%print_freq = 0.1;
print_freq = dt;
%print_freq = dt*50; 
forcing_freq = 5.0;

velocities_save = [];
temperature_save = [];
momentum_save = [];
momentum_savex = [];
pgrad_save = [];
pgrad_savex = [];

buoy_save = [];

Pr = 1.0;
jj = 2;
time = 0;

%Main2D;
udummy= ones(size(Xu));

RHSu_sq = ones(size(Xu));
RHSv_sq = ones(size(Xv));
RHST_sq = ones(size(Xls));


uones = ones(size(Uu));
vones = ones(size(Uv));

%     Cs = 0.22;
%     Nu = sparse(n,m); 
%   

for time = 1:size(time_steps,2)
    
      %SolveU:
      
      Ui = uatu; %uatu.*i_xu' + vatu.*i_yu';
      Uj = vatu; %uatu.*j_xu' + vatu.*j_yu';
    
      Ui = expand2D(Ui);
      Uj = expand2D(Uj);
      
      for j = 2:m
        for i = 2:n-1
         
              a_pls = max(Ui(i,j),0);
              a_neg = min(Uj(i,j),0);
                     
              u_pls = (uatu(i,j+1) - uatu(i,j))*Idiu;
              u_neg = (uatu(i,j) - uatu(i,j-1))*Idiu;
                           
              udu_di = a_pls*u_neg + a_neg*u_pls;
              
              u_pls = (uatu(i+1,j) - uatu(i,j))*Idju;
              u_neg = (uatu(i,j) - uatu(i-1,j))*Idju;
              
              vdu_dj = a_pls*u_neg + a_neg*u_pls;
                               
%               dispersion = (uatu(i,j+1) -2*uatu(i,j) +uatu(i,j-1))  /(Idiu*Idiu) +...
%                            (uatu(i+1,j) -2*uatu(i,j) +uatu(i-1,j))  /(Idju*Idju);
            
             RHSu_sq(i,j) = -udu_di -vdu_dj; %+ dispersion  ; 
             %RHSu_sq(i,j) = -udu_di*i_xu(j,i)' -vdu_dj; %*works!?
             %RHSu_sq(i,j) = -udu_di*i_xu(j,i)' -vdu_dj*j_xu(j,i)';      
              
        end
      end
       
%       i = 1;
%          for j = 2:m
%               u_pls = (uatu(i,j+1) - uatu(i,j))*Idiu;
%               u_neg = (uatu(i,j) - uatu(i,j-1))*Idiu;
% 
%               udu_di = a_pls*u_neg + a_neg*u_pls;
%               
%               u_pls = (uatu(i+1,j) - uatu(i,j))*Idju;
%               vdu_dj = a_neg*u_pls;
%               RHSu_sq(i,j) = -udu_di*i_xu(j,i)' -vdu_dj; %*j_xu(j,i)';      
%          end
%                                    
%       i = n;
%           for j = 2:m                   
%               u_neg = (uatu(i,j) - uatu(i,j-1))*Idiu;
%               u_pls = (uatu(i,j+1) - uatu(i,j))*Idiu;
% 
%               udu_di = a_pls*u_neg;
%               
%               u_neg = (uatu(i,j) - uatu(i-1,j))*Idju;
%               vdu_dj = a_neg*u_pls;
%               RHSu_sq(i,j) = -udu_di*i_xu(j,i)' -vdu_dj; %*j_xu(j,i)';      
%           end  
%           
%        j = 1;               
%            for i = 2:n-1
%               u_pls = (uatu(i,j+1) - uatu(i,j))*Idiu;
%               udu_di =  a_neg*u_pls;        
%               
%               u_pls = (uatu(i+1,j) - uatu(i,j))*Idju;
%               u_neg = (uatu(i,j) - uatu(i-1,j))*Idju;
%               vdu_dj = a_pls*u_neg + a_neg*u_pls;
%               RHSu_sq(i,j) = -udu_di*i_xu(j,i)' -vdu_dj; %*j_xu(j,i)';      
%            end
%            
%        j = m+1;
%             for i = 2:n-1          
%               u_neg = (uatu(i,j) - uatu(i,j-1))*Idiu;
%               udu_di = a_neg*u_pls;
%                       
%               u_pls = (uatu(i+1,j) - uatu(i,j))*Idju;
%               u_neg = (uatu(i,j) - uatu(i-1,j))*Idju;
%               vdu_dj = a_pls*u_neg + a_neg*u_pls;
%               RHSu_sq(i,j) = -udu_di*i_xu(j,i)' -vdu_dj; %*j_xu(j,i)';      
%                                  
%             end
      
      
     
     clear Ui; clear Uj; 
      
     Ui = uatv; %uatv.*i_xv' + vatv.*i_yv';
     Uj = vatv; %uatv.*j_xv' + vatv.*j_yv';
    
     Ui = expand2D(Ui);
     Uj = expand2D(Uj);
      
     for j = 2:m-1
        for i = 2:n
            
              u_pls = (vatv(i,j+1) - vatv(i,j))*Idiv;
              u_neg = (vatv(i,j) - vatv(i,j-1))*Idiv;
              
              a_pls = max(Ui(i,j),0);
              a_neg = min(Uj(i,j),0);
            
              udv_di = a_pls*u_neg + a_neg*u_pls;
              
              u_pls = (vatv(i+1,j) - vatv(i,j))*Idjv;
              u_neg = (vatv(i,j) - vatv(i-1,j))*Idjv;
              
              vdv_dj = a_pls*u_neg + a_neg*u_pls;
          
%               dispersion = (vatv(i,j+1) -2*vatv(i,j) +vatv(i,j-1))  /(Idiv*Idiv) +...
%                            (vatv(i+1,j) -2*vatv(i,j) +vatv(i-1,j))  /(Idjv*Idjv);          
                        
              RHSv_sq(i,j) = -udv_di -vdv_dj -batv(i,j); %+ dispersion; 
              %RHSv_sq(i,j) = -udv_di*i_yv(j,i)' -vdv_dj*j_yv(j,i)' -batv(i,j)*Jv(j,i)'; 
              %RHSv_sq(i,j) = -udv_di*i_yv(j,i)' -vdv_dj*j_yv(j,i)'; % -batv(i,j); 
              
        end
     end

    RHSu_sq = expand2D(RHSu_sq);
    RHSv_sq = expand2D(RHSv_sq);
      
    RHSu = reshape(RHSu_sq.',[], 1);
    RHSv = reshape(RHSv_sq.',[], 1);
           
    Uu = SSPRK101D(Uu,RHSu,dt);
    Vv = SSPRK101D(Vv,RHSv,dt);   
    
    u_star = reshape(Uu,layu).' ;    
    v_star = reshape(Vv,layv).'; 
    
    %%%%u,v prediction boundaries:
    [u_star,v_star,~] = applyboundaries2Dbeam(u_star,v_star,p);
       
    force = uo*sin(omega*(current_time*LStar/UStar));
    
    %u_star(1,:) =   force*udummy(1,:);
    %u_star(end,:) = force*udummy(end,:);
    
    %u_star(:,1) = u_star(:,1) + force*udummy(:,1);
    %u_star(:,end) = u_star(:,end) + force*udummy(:,end);
    
%     fig = figure(jj);
%     u_v = squeeze(u_star(:,:)*UStar);
%     %u_v = squeeze(udvdx_sq(:,:)*UStar);
%     subplot(2,1,1)
%     %pcolor(Xu(:,:),Yu(:,:),u_v);
%     pcolor(u_v);
%     axis tight
%     title('u [m/s]')
%     colorbar;
%     shading interp;
%     subplot(2,1,2)
%     u_3 = squeeze(v_star(:,:)*UStar);
%     %u_3 = squeeze(vdvdy_sq(:,:)*UStar);
%     %pcolor(Xv(:,:),Yv(:,:),u_3);
%     pcolor(u_3);
%     axis tight
%     title('v [m/s]')
%     colorbar;
%     shading interp;
%    
  %  
    Uu = reshape(u_star.',[], 1);
    Vv = reshape(v_star.',[], 1);
    
    %laplacian solver:
    R = [Uu; Vv];
    
    DR = D*R;
% 
    if current_time ~= 1   
        DR(bdry) = Pp(bdry);    
    end
    
    %Rorig = D*R*rho0/dt;   
    %Pp = L\Rorig;

    Pp = L\(DR); 
    
    p = reshape(Pp,[m+2 n+2]).';    
    [~,~,p] = applyboundaries2Dbeam(uatu,vatv,p);
    Pp = reshape(p.',[], 1);                                                                                  
    
%    figure(888)
%     pcolor(p)
%     shading flat
%     colorbar
    
%    
    %Pp = Pp*LStar;    
    
    %Pp = Pp*rhoStar;    
    
    %gradient of pressure:
    Gpp = Gp*[ Pp ];
    dpdx = Gpp(1:(m+1)*n);
    dpdy = Gpp(((m+1)*n+1):end);   
    
    dpdx_sq = reshape(dpdx,[m+1 n]).' ;
    dpdy_sq = reshape(dpdy,[m n+1]).' ;  
%     
%     fig = figure(jj+1);
%     u_v = squeeze(dpdx_sq(:,:)*UStar);
%     %u_v = squeeze(udvdx_sq(:,:)*UStar);
%     subplot(2,1,1)
%     %pcolor(Xu(:,:),Yu(:,:),u_v);
%     pcolor(u_v);
%     axis tight
%     title('dpdx')
%     colorbar;
%     shading flat; %interp;
%     subplot(2,1,2)
%     u_3 = squeeze(dpdy_sq(:,:)*UStar);
%     %u_3 = squeeze(vdvdy_sq(:,:)*UStar);
%     %pcolor(Xv(:,:),Yv(:,:),u_3);
%     pcolor(u_3);
%     axis tight
%     title('dpdy')
%     colorbar;
%     shading flat %;interp;    
% 
    %u,v correction:
    %Uu = Uu + dt*dpdx;
    %Vv = Vv + dt*dpdy;
   
    
    Uu = SSPRK101D(Uu,dpdx,dt);
    Vv = SSPRK101D(Vv,dpdy,dt);  
    
    uatu = reshape(Uu,layu).' ;
    vatv = reshape(Vv,layv).' ;
      
    %velocities boundaries:
    [uatu,vatv,~] = applyboundaries2Dbeam(uatu,vatv,p);
    
    
%     fig = figure(jj+2);
%     u_v = squeeze(uatu(:,:)*UStar);
%     subplot(2,1,1)
%     pcolor(Xu(:,:),Yu(:,:),u_v);
%     axis tight
%     title('u [m/s]')
%     colorbar;
%     shading interp;
%     subplot(2,1,2)
%     u_3 = squeeze(vatv(:,:)*UStar);
%     pcolor(Xv(:,:),Yv(:,:),u_3);
%     axis tight
%     title('v [m/s]')
%     colorbar;
%     shading interp;
       

%Shock absorbing boundaries:
    %Adding sponge Layers to the code for Beam Experiment
    %Ref: A nonhydrostatic, isopycnal-coordinate ocean model for internal waves. S. Vitousek, O.B. Fringer / Ocean Modelling 83 (2014)
    % equation 53
%     Ls=300.0d0/LStar;
%     %Ls=1; %10.0d0/LStar;
%     taus=100.0d0*UStar/LStar;
     half_x1 = floor(m/6);
%     
%     %West Side
% 
%     for j=1:half_x1-1
%     for i=1:n
%             r(i,j)= Xu(i,j) - Xu(i,1);
%             uatu(i,j)= uatu(i,j) - dt*(uatu(i,j) - uatu(i,end))*exp(-4.0d0*r(i,j)/Ls)/(taus);
%     end
%     end
% 
%     %East Side
%     for j=m-half_x1:m
%     for i=1:n
%             r(i,j)= Xu(i,m) - Xu(i,j);
%             uatu(i,j)= uatu(i,j) - dt*(uatu(i,j) - uatu(i,1))*exp(-4.0d0*r(i,j)/Ls)/(taus);
%     end
%     end
    
    vatu = InterpVtoV2D(vatv',m,n,'vtou');
    uatv = InterpVtoV2D(uatu',m,n,'utov'); 
    
    Vu = reshape(vatu,[], 1);
    Uv = reshape(uatv,[], 1);

    vatu = vatu';
    uatv = uatv';
%     
    %Convection-Difussion equation:
    
    uvec = [Uu;Vv];

    Tvec = reshape(T.',[], 1);

     for i=1:m
     for j=1:n
         uatl(j,i) = (uatu(j,i) + uatu(j,i+1))/2;
         vatl(j,i) = (vatv(j,i) + vatv(j+1,i))/2;
     end
     end
%   
%     uatl(2,:) = uatl(3,:);
%     uatl(1,:) = uatl(2,:);
%     uatl(end-1,:) = uatl(end-2,:);
%     uatl(end,:) = uatl(end-1,:);
%      
%     uatl(:,2) = uatl(:,3);
%     uatl(:,1) = uatl(:,2);
%     uatl(:,end-1) = uatl(:,end-2);
%     uatl(:,end) = uatl(:,end-1);
% 
%     vatl(2,:) = vatl(3,:);
%     vatl(1,:) = vatl(2,:);
%     vatl(end-1,:) = vatl(end-2,:);
%     vatl(end,:) = vatl(end-1,:);
%      
%     vatl(:,2) = vatl(:,3);
%     vatl(:,1) = vatl(:,2);
%     vatl(:,end-1) = vatl(:,end-2);
%     vatl(:,end) = vatl(:,end-1);

% % %     uatl = expand2D(uatl);
% % %     vatl = expand2D(vatl);
% % %     
% % %     Ti = uatl.*i_xl' + vatl.*i_yl';
% % %     Tj = uatl.*j_xl' + vatl.*j_yl';
% % %     
% % %     Ti = expand2D(Ti);
% % %     Tj = expand2D(Tj);
% % %       
% % %     for j = 2:m
% % %      for i = 2:n
% % % 
% % %           T_pls = (T(i,j+1) - T(i,j))*Idil;
% % %           T_neg = (T(i,j) - T(i,j-1))*Idil;
% % % 
% % %           a_pls = max(Ti(i,j),0);
% % %           a_neg = min(Tj(i,j),0);
% % % 
% % %           udT_di = a_pls*T_neg + a_neg*T_pls;
% % % 
% % %           T_pls = (T(i+1,j) - T(i,j))*Idjl;
% % %           T_neg = (T(i,j) - T(i-1,j))*Idjl;
% % % 
% % %           vdT_dj = a_pls*T_neg + a_neg*T_pls;
% % % 
% % %           %RHSu_sq(i,j) = -udT_di*i_xl(j,i)' -vdT_dj*j_xl(j,i)'; 
% % %           %RHST_sq(i,j) = -udT_di*i_xl(j,i)' -vdT_dj; %*works!?
% % % 
% % %           RHST_sq(i,j) = -udT_di -vdT_dj; %*works!?
% % % 
% % %           
% % %       end
% % %     end
% % %     
% % %     RHST = reshape(RHST_sq.',[], 1);
% % %        
% % %     T_new = SSPRK101D(Tvec,RHST,dt);
% % %     T     = reshape(T_new,[m+2 n+2]).';   
% % %        
% % %     %Apply boundaries to T:
% % %     
% % %     T(2,:) = T(3,:);
% % %     T(1,:) = T(2,:);
% % %     T(end-1,:) = T(end-2,:);
% % %     T(end,:) = T(end-1,:);
% % %      
% % %     T(:,2) = T(:,3);
% % %     T(:,1) = T(:,2);
% % %     T(:,end-1) = T(:,end-2);
% % %     T(:,end) = T(:,end-1);
% % %     
% % %     
% % %     %dens = rho_ref + drho_dT*(T*T0-T_ref);  %EOS   
% % %     %rhoprime = dens-rho_P0;
% % %     %buoy = Ri*rhoprime;
% % %     
% % %     dens = rho_ref + drho_dT*(T-T_ref);  %EOS   
% % %     %rhoprime = dens-rho_P0;
% % %     %dens = rho_ref - alpha*rho_ref*(T-T_ref);
% % %     %dens = rho_ref - alpha*rho_ref*(T*T0-T_ref);
% % %     
% % %     %dens = rho_ref*( 1 - alpha*(T-T_ref));
% % %     buoy = gforce*(dens-rho0)/rho0;
% % %         
% % %     for i=2:m
% % %         for j=2:n+1
% % %             batv(j,i) = (buoy(j,i) + buoy(j-1,i))/2;
% % %         end
% % %     end
%     
%     batv(1,:) = batv(2,:);
%     batv(:,1) = batv(:,2);
%     batv(end,:) = batv(end-1,:);
          
    current_time = (time+1)*dt*LStar/UStar
    mod(current_time,print_freq);
if(mod(current_time,print_freq) <= 1e-10 )
    
%     fig = figure(jj+5);
%     %set(gcf,'units','inches','position',[10 10 15 10])
%     u_v = squeeze(uatu(:,:)*UStar);
%     subplot(2,2,1)
%     pcolor(Xu(:,:),Yu(:,:),u_v);
%     %pcolor(u_v);
%     %axis equal
%     axis tight
%     title('u [m/s]')
%     colorbar;
%     shading flat %interp;
%     subplot(2,2,2)
%     u_3 = squeeze(vatv(:,:)*UStar);
%     pcolor(Xv(:,:),Yv(:,:),u_3);
%     %pcolor(u_3);
%     %axis equal
%     axis tight
%     title('v [m/s]')
%     colorbar;
%     shading flat %interp;
%     subplot(2,2,3)
%     u_2 = squeeze(batv(:,:));
%     pcolor(Xv(:,:),Yv(:,:),u_2);
%     %pcolor(u_2);
%     %axis equal
%     axis tight 
%     title('buoy')
%     colorbar;
%     shading flat %interp;
%     subplot(2,2,4)
%     u_2 = squeeze(T(:,:));
%     pcolor(Xls(:,:),Yls(:,:),u_2);
%     %pcolor(u_2);
%     %axis equal
%     axis tight 
%     title(['T [C] at ' num2str(current_time) 's'])
%     colorbar;
%     shading flat %interp;

       
    temperature_save = [temperature_save; {T}];
    velocities_save = [velocities_save; {uatu}];
    
half_x = floor(IMax/2);

tenth_y = floor(JMax/10);
    
% 
    fig = figure(jj+4);
    set(gcf,'units','inches','position',[10 10 15 10])
    subplot(2,1,1)
    u_v = squeeze(uatu(:,:)*UStar);
    pcolor(Xu(:,half_x1:end-half_x1),Yu(:,half_x1:end-half_x1),u_v(:,half_x1:end-half_x1)-force);
    %pcolor(u_v);
    %contourf(u_v,25);
    %axis equal
    axis tight
    title(['u [m/s] at ' num2str(current_time) 's'])
    colorbar;
    shading interp;
    subplot(2,1,2)
    u_v = squeeze(uatu(:,:)*UStar);
    pcolor(Xu(1:tenth_y,half_x-10:half_x+10),Yu(1:tenth_y,half_x-10:half_x+10),u_v(1:tenth_y,half_x-10:half_x+10)-force);
    %pcolor(u_v);
    %contourf(u_v,25);
    %axis equal
    axis tight
    title(['u [m/s] at ' num2str(current_time) 's'])
    colorbar;
    shading interp;

%     fig = figure(jj+4);
%     set(gcf,'units','inches','position',[10 10 15 10])
%     u_v = squeeze(uatu(:,:)*UStar);
%     pcolor(Xu,Yu,u_v-force);
%     %pcolor(u_v);
%     %contourf(u_v,25);
%     %axis equal
%     axis tight
%     title(['u [m/s] at ' num2str(current_time) 's'])
%     colorbar;
%     shading interp;


    
    
    %f_name = ['LE_fullmimetic_' num2str(JMax) '_dt.01_' num2str(current_time) 's.png'];
    %saveas(fig,f_name,'png')


    %current_time
    
 end
    
end

%%

%load('beam2D_101x131_wN_0.2_finished_incorrect.mat')


%%

[N,M] = size(velocities_save)

%v = VideoWriter ('beam_101x131_wN_0.2.avi');
%open(v);  
II = 0;
for i = 1:1:N

vel_n = velocities_save{i,1};

fig = figure(876);
%set(gcf,'units','inches','position',[10 10 15 10])
u_v = squeeze(vel_n(:,:)*UStar);
%pcolor(Xu(:,:),Yu(:,:),u_v-force);
%pcolor(u_v);
contourf(Xu(:,:),Yu(:,:),u_v,25);
%axis equal
axis tight
title(['u [m/s] at ' num2str(i*print_freq*LStar/UStar) 's'])
colorbar;
shading interp;


drawnow;
%frame = getframe (gcf);
%writeVideo (v, frame);  
hold off;


end  

%close(v);

%%

    Uu_sq = reshape(Vu,layu).' ;   
    Uv_sq = reshape(Uv,layv).' ;

    fig = figure(jj+5);
    %u_v = squeeze(u_star(:,:)*UStar);
    u_v = squeeze(Uu_sq(:,:));
    subplot(2,1,1)
    %pcolor(Xu(:,:),Yu(:,:),u_v');
    pcolor(u_v);
    axis tight
    title('u [m/s]')
    colorbar;
    shading interp;
    subplot(2,1,2)
    %u_3 = squeeze(v_star(:,:)*UStar);
    u_3 = squeeze(Uv_sq(:,:));
    %pcolor(Xv(:,:),Yv(:,:),u_3');
    pcolor(u_3);
    axis tight
    title('v [m/s]')
    colorbar;
    shading interp;


%%

    Uu = reshape(uatu.',[], 1);
    Vv = reshape(vatv.',[], 1);
    
    vatu = InterpVtoV2D(vatv',m,n,'vtou')'; 
    %Vu = vatu(:);
    
    uatv = InterpVtoV2D(uatu',m,n,'utov')'; 
    
    Vu = reshape(vatu.',[], 1);
    Uv = reshape(uatv.',[], 1);


    fig = figure(jj+4);
    %set(gcf,'units','inches','position',[10 10 15 10])
    u_v = squeeze(uatu(:,:)*UStar);
    subplot(2,2,1)
    %pcolor(Xu(:,:),Yu(:,:),u_v);
    pcolor(u_v);
    %axis equal
    axis tight
    title('uatu')
    colorbar;
    shading interp;
    subplot(2,2,2)
    u_3 = squeeze(vatv(:,:)*UStar);
    %pcolor(Xv(:,:),Yv(:,:),u_3);
    pcolor(u_3);
    %axis equal
    axis tight
    title('vatv')
    colorbar;
    shading interp;
    subplot(2,2,3)
    u_2 = squeeze(vatu(:,:));
    %pcolor(Xv(:,:),Yv(:,:),u_2);
    pcolor(u_2);
    %axis equal
    axis tight 
    title('vatu')
    colorbar;
    shading interp;
    subplot(2,2,4)
    u_2 = squeeze(uatv(:,:));
    %pcolor(Xls(:,:),Yls(:,:),u_2);
    pcolor(u_2);
    %axis equal
    axis tight 
    title(['uatv'])
    colorbar;
    shading interp;


