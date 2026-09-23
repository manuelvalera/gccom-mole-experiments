%close all; clear all; clc

%IMax=101; JMax = 10; KMax = 38;

IMax=50; JMax = 10; KMax = 38;

m = IMax-1;
o = JMax-1; 
n = KMax-1;

Ho=1; 
ab= 0.5 ;
bb = 0.7;
%bb = 1.;
cc = 1.4;
Lb=8;
D=Ho;
%
LStar = 1000.0;
UStar = 1.0;

%a=-1.8; b=1.8;
a=0; b=3.6;
c=0  ; d=2.8;
%c=-1.4  ; d=1.4;
e=-Ho ; f=0;
%e=Ho ; f=0;

% Dx2 = (b-a)/m ;  
% Dy2 = (f-e)/n ;
% Dz2 = (d-c)/o ;

%Define Grid bathymetry / geometry:

%bathy = @(x,z) (-Ho*( 1 + ab*exp( -Lb*( (x - bb).^2 + z.^2 ) )));
bathy = @(x,z) (-Ho*( 1 + ab*exp( -Lb*( (x - bb).^2 + (z - cc).^2 ) )));

%bathy = @(x,z) (LStar*( -1 + ab*exp( -Lb*( (x - bb).^2 + (z - cc).^2 ) )));


%Create grid:
xspace = linspace(a,b,m);
zspace = linspace(c,d,o);
yspace = linspace(e,f,n);

%[X2, Y2, Z2] = meshgrid(a:Dx2:b, e:Dy2:f, c:Dz2:d);

[X, Y, Z] = meshgrid(xspace,yspace,zspace);

Dx = X(1,2,1) - X(1,1,1);
Dz = Y(2,1,1) - Y(1,1,1);
Dy = Z(1,1,2) - Z(1,1,1);

%Apply bathymetry / geometry:
Y = -Ho*Y ./ ( bathy(X,Z) ); %sigma transformation

%Y = Y.*( bathy(X,Z) );

% %%
%
% figure(2)
% subplot(2,1,1)
% for k=1:KMax-1
% plot(squeeze(X(k,:,floor(JMax/2))),squeeze(Y(k,:,floor(JMax/2))),'-*')
% hold on
% end
% xlabel('i')
% ylabel('k')
% %hold off
% subplot(2,1,2)
% %figure(3)
% for k=1:KMax-1
% plot(squeeze(Z(k,floor(IMax/4),:)),squeeze(Y(k,floor(IMax/4),:)),'-*')
% %plot(squeeze(X(k,floor(IMax/2),:)),squeeze(Y(k,floor(IMax/2),:)),'-*')
% hold on
% end
% xlabel('i')
% ylabel('j')
% hold off
%
% %%
% %half_x = floor(KMax/2);
% half_x = 1;
% half_y = floor(IMax/2);
% half_z = floor(JMax/2);
% 
% figure(101)
%       colormap parula
%       brighten(.7);
%       surf(X(:,:,half_z),Y(:,:,half_z),Z(:,:,half_z)); %plot the bottom layer
%       hold on
%       surf(squeeze(X(half_x,:,:)),squeeze(Y(half_x,:,:)),squeeze(Z(half_x,:,:))); %plot the layer that chops the region in half wrt the x-axis
%       hold on
%       surf(squeeze(X(:,half_y,:)), squeeze(Y(:,half_y,:)), squeeze(Z(:,half_y,:)));
%       xlabel('x')
%       ylabel('y')
%       zlabel('z')
%       %hcb=colorbar('east');
%       %      title(hcb,'T[C]')
%       hold off     
%       fig = gcf;
%       view(170,-72)

%

%X2 = X(1:end-1,:,:);

Ux = (X(1:end-1, :, :) + X(2:end, :, :))/2;
Xu = (Ux(:, :, 1:end-1) + Ux(:, :, 2:end))/2;

u_max = 0.01

forcing = squeeze(Xu(:,1,:));

    for k=1:o-1
        for j=1:n-1
            forcing(j,k) = u_max*j/KMax; %Real forcing
        end
    end

%
X2 = X;
beta = 5;

Lx= 3.6; Ly = 2.8; Lz = 1;
D = 0.5*Lx;
A = (1/(2*beta))*log( (1+(exp(beta)-1)*D/Lx)./(1+(exp(-beta)-1)*D/Lx));

beta2 = linspace(5,0,n+1);
%
for k=1:KMax-1
for j=1:IMax-1
    xi = (j-1)/(IMax-1);
    X2(k,j,:) = D*(1 + sinh(beta2(k)*(xi-A))./sinh(beta2(k)*A)) - Lx/2;
    %X2(k,j,:) = D*(1 + sinh(beta*(xi-A))./sinh(beta*A)) - Lx/2;
    
end
end
%

Z2 = Z;

D = 0.5*Ly;
A = (1/(2*beta))*log( (1+(exp(beta)-1)*D/Ly)./(1+(exp(-beta)-1)*D/Ly));
y1d = zeros(JMax,1);
for k=1:KMax-1
for j=1:JMax-1
    xi = (j-1)/(JMax-1);
    Z2(k,:,j) = D*(1 + sinh(beta2(k)*(xi-A))./sinh(beta2(k)*A)) - Ly/2;
end
end

%%
figure(2)
subplot(2,1,1)
for k=1:KMax-1
plot(squeeze(X2(k,:,floor(JMax/2))),squeeze(Y(k,:,floor(JMax/2))),'-*')
hold on
end
xlabel('i')
ylabel('k')
title(['Seamount 3D grid (front) Y=5 - [X,Y,Z]=[101x7x50]'])
hold off
subplot(2,1,2)
%figure(3)
for k=1:KMax-1
plot(squeeze(Z2(k,floor(IMax/4),:)),squeeze(Y(k,floor(IMax/4),:)),'-*')
%plot(squeeze(X(k,floor(IMax/2),:)),squeeze(Y(k,floor(IMax/2),:)),'-*')
hold on
end
xlabel('j')
ylabel('k')
title(['Seamount 3D grid (side) Y=5 - [X,Y,Z]=[101x7x50]'])
hold off

%%
%Grid Plots:
% figure(2)
% for k=1:KMax-1
% %plot(squeeze(X2(k,floor(IMax/2),:)),squeeze(Y(k,floor(IMax/2),:)),'-*')
% plot(squeeze(X2(k,:,floor(JMax/2))),squeeze(Y(k,:,floor(JMax/2))),'-*')
% hold on
% end
% hold off
% 
% figure(3)
% for k=1:KMax-1
% plot(squeeze(Z2(k,floor(IMax/4),:)),squeeze(Y(k,floor(IMax/4),:)),'-*')    
% %plot(squeeze(Z2(k,floor(IMax/2),:)),squeeze(Y(k,floor(IMax/2),:)),'-*')
% %plot(squeeze(Z2(k,:,floor(JMax/2))),squeeze(Y(k,:,floor(JMax/2))),'-*')
% hold on
% end
% hold off
% % 
% %
% 
% %half_x = floor(KMax/2);
half_x = 1;
half_y = floor(IMax/2);
half_z = floor(JMax/2);

figure(101)
      colormap parula
      brighten(.7);
      surf(X2(:,:,half_z),Y(:,:,half_z),Z2(:,:,half_z)); %plot the bottom layer
      hold on
      surf(squeeze(X2(half_x,:,:)),squeeze(Y(half_x,:,:)),squeeze(Z2(half_x,:,:))); %plot the layer that chops the region in half wrt the x-axis
      hold on
      surf(squeeze(X2(:,half_y,:)), squeeze(Y(:,half_y,:)), squeeze(Z2(:,half_y,:)));
      xlabel('x')
      ylabel('y')
      zlabel('z')
      %hcb=colorbar('east');
      %      title(hcb,'T[C]')
      hold off     
      fig = gcf;
      view(170,-72)
%%
%
%
%Copying curvilinear to be default:  
Z = Z2;
X = X2; 
      

%
%Create logical grid:
% [Xs, Ys, Zs] = meshgrid([a a+Dx/2 : Dx : b-Dx/2 b],...
%                         [e e+Dy/2 : Dy : f-Dy/2 f],...
%                         [c c+Dz/2 : Dz : d-Dz/2 d]);
% %                          
% Xs = permute(Xs, [2, 3, 1]);
% Ys = permute(Ys, [2, 3, 1]);
% Zs = permute(Zs, [2, 3, 1]);
% 
% Ys = -Ho*Ys ./ ( bathy(Xs,Zs) ); %sigma transformation

%
%Create operators:

full_ct = (m+2)*(n+2)*(o+2);

jj = 2;
slice = 3;
time = 0;

m = m-1;
n = n-1;
o = o-1;

LStar = 1;
UStar = 1;

% Get 3D curvilinear mimetic divergence
D = div3DCurv(2, X, Y, Z)*LStar;
G = grad3DCurv(2, X, Y, Z)*LStar;
%N = neumann3DCurv(G, m, n, o, 1);
N = LStar^2*robinBC3D(2,m,Dx,n,Dy,o,Dz,1,0); %Dirichlet
%N = robinBC3D(2,m,Dx,n,Dy,o,Dz,0,1); %Neumann

L = D*G + N;
%cond(L)

%
%Create C-type grids:

Ux = (X(1:end-1, :, :) + X(2:end, :, :))/2;
Xu = (Ux(:, :, 1:end-1) + Ux(:, :, 2:end))/2;
Uy = (Y(1:end-1, :, :) + Y(2:end, :, :))/2;
Yu = (Uy(:, :, 1:end-1) + Uy(:, :, 2:end))/2;
Uz = (Z(1:end-1, :, :) + Z(2:end, :, :))/2;
Zu = (Uz(:, :, 1:end-1) + Uz(:, :, 2:end))/2;
Vx = (X(:, 1:end-1, :) + X(:, 2:end, :))/2;
Xv = (Vx(:, :, 1:end-1) + Vx(:, :, 2:end))/2;
Vy = (Y(:, 1:end-1, :) + Y(:, 2:end, :))/2;
Yv = (Vy(:, :, 1:end-1) + Vy(:, :, 2:end))/2;
Vz = (Z(:, 1:end-1, :) + Z(:, 2:end, :))/2;
Zv = (Vz(:, :, 1:end-1) + Vz(:, :, 2:end))/2;
Wx = (X(1:end-1, :, :) + X(2:end, :, :))/2;
Xw = (Wx(:, 1:end-1, :) + Wx(:, 2:end, :))/2;
Wy = (Y(1:end-1, :, :) + Y(2:end, :, :))/2;
Yw = (Wy(:, 1:end-1, :) + Wy(:, 2:end, :))/2;
Wz = (Z(1:end-1, :, :) + Z(2:end, :, :))/2;
Zw = (Wz(:, 1:end-1, :) + Wz(:, 2:end, :))/2;

% Compute centroids for visualization
Xviz = (X(1:end-1, :, :)+X(2:end, :, :))/2;
Xviz = (Xviz(:, 1:end-1, :)+Xviz(:, 2:end, :))/2;
Xviz = (Xviz(:, :, 1:end-1)+Xviz(:, :, 2:end))/2;
Yviz = (Y(1:end-1, :, :)+Y(2:end, :, :))/2;
Yviz = (Yviz(:, 1:end-1, :)+Yviz(:, 2:end, :))/2;
Yviz = (Yviz(:, :, 1:end-1)+Yviz(:, :, 2:end))/2;
Zviz = (Z(1:end-1, :, :)+Z(2:end, :, :))/2;
Zviz = (Zviz(:, 1:end-1, :)+Zviz(:, 2:end, :))/2;
Zviz = (Zviz(:, :, 1:end-1)+Zviz(:, :, 2:end))/2;

%
%
X = LStar*X; Xu = LStar*Xu; Xv = LStar*Xv; Xw = LStar*Xw;
Y = LStar*Y; Yu = LStar*Yu; Yv = LStar*Yv; Yw = LStar*Yw;
Z = LStar*Z; Zu = LStar*Zu; Zv = LStar*Zv; Zw = LStar*Zw;

%%


% xspace = linspace(a,b,m+1);
% yspace = linspace(c,d,o+1);
% zspace = linspace(e,f,n+1);
% 
% [Xg, Yg, Zg] = meshgrid(xspace,zspace,yspace);
% 
% Yg = -Ho*Yg ./ ( bathy(Xg,Zg) ); %sigma transformation
% 
% 
% Gg = grad3DCurv(2, Xg, Yg, Zg)*LStar;





%%

LStar = 1;
UStar = 1;

%Time discretization (needed for D and G):
tf = 12000.0; 
%tf = 1000.0; 
%tf = Tbc*2; 
%tf = 3500.0; 
%tf = 22320.0; %GCCOM tf 
%tf = 19.25*Tbc;
%tf = dt*5;
tf = tf*(UStar/LStar);
dt = 10.0;
%dt = Tbc/500;

%CFL dt limit:
%maxv = UStar;
%dt = (min(Dx,Dy)/maxv)/2;

dt = dt*(UStar/LStar);
nf = tf/dt ;
time_steps = linspace(0,tf,nf);
rho_ref = 1027.0;



%print_freq = 0.1;
%print_freq = dt;
print_freq = dt*50; 
forcing_freq = 5.0;

%Gp = dt/rho_ref*G;
Gp = -dt/rho_ref*G;
%Gp = dt./G;

%Gp = G;
%
velocities_save = [];
vvelocities_save = [];
wvelocities_save = [];
uatu = zeros(m+1,n,o);

%u_max=1.0; u_min = 0.0;

u_max=0.01; u_min = 0.0;

%Initial conditions for u-v-w:
Ugiven = zeros(size(X)); 
%Ugiven = X.^2;
%Ugiven = u_max*sin(X); %X.^2;

%
% for k=1:o
%     for i=1:n
%         Ugiven(i,1:3,k) = u_max*i/KMax; 
%         %Ugiven(i,1,k) = u_max*X(i,1,k); 
%     end
% end

%

%uatu2 = permute(reshape(Ugiven, m+1, n+1, o+1), [2, 1, 3]);

sz = 3;
sy = 25;
sl = 3;

nlevs = 20;

% figure(200) 
% subplot(2,1,1)
% pcolor(squeeze(X(:, 1:end, sl)), squeeze(Y(:, 1:end, sl)), squeeze(Ugiven(:, 1:end, sz)*UStar))  %visualize X-Y plane
% %contourf(squeeze(X(:, 1:end, sl)), squeeze(Y(:, 1:end, sl)), squeeze(uatu2(:, 1:end, sz)*UStar), nlevs)  %visualize X-Y plane
% xlabel('x'); ylabel('y')
% colorbar;   
% shading interp
% title(['U-velocity X-Y plane (front)'])
% subplot(2,1,2)
% pcolor(squeeze(X(2,1:end,:)), squeeze(Z(2,1:end,:)), squeeze(Ugiven(2,1:end,:)*UStar))  %visualize X-Y plane
% xlabel('x'); ylabel('z')
% colorbar
% title(['U-velocity X-Z plane (top)'])
% shading interp
%    

%
%
% figure(888) 
% pcolor(squeeze(Ugiven(1, :, :)))  %visualize X-Y plane
% xlabel('x'); ylabel('y')
% colorbar
% shading interp
% title(['p X-Y plane (front)'])
%


Vgiven = zeros(size(Y)); 
%Vgiven = u_max*sin(Y);
%Vgiven = 4.*Y;
%Vgiven = Y.^2;


Wgiven = zeros(size(Z)); 
%Wgiven = sin(Z);
%Wgiven = Z.^2;
%Wgiven = 2.*Z;


%Interpolate U,V,W to its own space:
interpolant = scatteredInterpolant([X(:) Y(:) Z(:)], Ugiven(:));
Uu0 = interpolant(Xu, Yu, Zu);
interpolant = scatteredInterpolant([X(:) Y(:) Z(:)], Vgiven(:));
Vv0 = interpolant(Xv, Yv, Zv);
interpolant = scatteredInterpolant([X(:) Y(:) Z(:)], Wgiven(:));
Ww0 = interpolant(Xw, Yw, Zw);


%Reshape to vector:
Uu = reshape(permute(Uu0, [2, 1, 3]), [], 1);
Vv = reshape(permute(Vv0, [2, 1, 3]), [], 1);
Ww = reshape(permute(Ww0, [2, 1, 3]), [], 1);

RHSu = Uu; Vu = Uu; Wu = Uu;
RHSv = Vv; Uv = Vv; Wv = Vv;
RHSw = Ww; Uw = Ww; Vw = Ww;

u_count = (m+1)*n*o;
v_count = u_count + (n+1)*m*o;

%Stratification:

%Tic = ones(IMax+1,JMax+1,KMax+1);


Tic = ones(size(X,1),size(X,2),size(X,3));
TMIN = 180; TMAX = 179.9705;
%Tic = TMIN + Y*(TMIN-TMAX);

Tic = TMIN + Y*(TMIN-TMAX);

S_ref   = 35.0;
T_ref   = 10.0;
drho_dS = 0.7810;
drho_dT = -0.17080;
gforce = 9.8;
%
%
dens = rho_ref + drho_dT*(Tic - T_ref);
rho0=mean(dens(:));     

buoy = gforce*(dens-rho0)/rho0;

%batv = gforce*(dens-rho0)/rho0;
Fbuoy = scatteredInterpolant([X(:) Y(:) Z(:)], buoy(:));
batv = Fbuoy(Xv,Yv,Zv);

%Bb = reshape(permute(batv, [2, 1, 3]), [], 1);

%
Y_in = permute(Y, [2, 1, 3]);

Xu2 = permute(Xu, [2, 1, 3]);
Yu2 = permute(Yu, [2, 1, 3]);

Xv2 = permute(Xv, [2, 1, 3]);
Yv2 = permute(Yv, [2, 1, 3]);

Xw2 = permute(Xw, [2, 1, 3]);
Yw2 = permute(Yw, [2, 1, 3]);
%

% slice = 2;
% figure(888)
% subplot(2,1,1)
% u_3 = squeeze(buoy(:,:,slice));
% pcolor(squeeze(X(:,:,slice)),squeeze(Y(:,:,slice)),u_3);
% shading interp
% colorbar
% subplot(2,1,2)
% u_3 = squeeze(batv(:,:,slice));
% pcolor(squeeze(Xv(:,:,slice)),squeeze(Yv(:,:,slice)),u_3);
% shading interp
% colorbar

% Not changed with new grid construction   
% I0 = interpol3D(m, n, o, 0, 0, 0);
% I1 = interpol3D(m, n, o, 1, 1, 1);
% 
% DI0 = interpolD3D(m, n, o, 0, 0, 0);
% DI1 = interpolD3D(m, n, o, 1, 1, 1);
% 
% DI0_left  = DI0(:,1:u_count);
% DI0_middle = DI0(:,u_count+1:v_count);
% DI0_right = DI0(:,v_count+1:end);
% 
% DI1_left  = DI1(:,1:u_count);
% DI1_middle = DI1(:,u_count+1:v_count);
% DI1_right = DI1(:,v_count+1:end);

temperature_save = [];

[Gub, Gvb, Gwb] = d3D(m, n, o, Dx, Dy, Dz, 0, 0, 0); 
[Guf, Gvf, Gwf] = d3D(m, n, o, Dx, Dy, Dz, 1, 1, 1); 
T_new = zeros(m+2,n+2,o+2); 
p = T_new;
udummy= ones(m+1,n,o);

current_time = (time+1)*dt*LStar/UStar

uatu = reshape(Uu,[m+1 n o]);
vatv = reshape(Vv,[m n+1 o]);
watw = reshape(Ww,[m n o+1]);  

%[sgs_u,sgs_v,sgs_w] = CalcSGS(G,D,uatu,vatv,watw,m,n,o,Dx,Dy,Dz);

% Boussinesq Pred-Corr formulation:
for time = 1:size(time_steps,2) 
    RHST = zeros(m+2,n+2,o+2);
    %Predict velocities:
    %Bb = reshape(permute(batv, [2, 1, 3]), [], 1);
    Bb   = batv(:); 
    
    %'u'
    [udu_dx,vdudy_atv,wdudz_atw] = MimeticUpwindD3D(Uu,Uv,Uw,Uu,Vv,Ww,Guf,Gub,Gvf,Gvb,Gwf,Gwb,m,n,o);   

    vdudy_atv = reshape(full(vdudy_atv),[m n+1 o]);
    vdu_dy    = InterpVtoV(vdudy_atv,m,n,o,'vtou'); 
    

    
    wdudz_atw = reshape(full(wdudz_atw),[m n o+1]);   
    wdu_dz    = InterpVtoV(wdudz_atw,m,n,o,'wtou'); 
    
    %'v'
    [udvdx_atu,vdv_dy,wdvdz_atw] = MimeticUpwindD3D(Vu,Vv,Vw,Uu,Vv,Ww,Guf,Gub,Gvf,Gvb,Gwf,Gwb,m,n,o);
  
    udvdx_atu = reshape(full(udvdx_atu),[m+1 n o]);   
    udv_dx    = InterpVtoV(udvdx_atu,m,n,o,'utov');
 

    wdvdz_atw = reshape(full(wdvdz_atw),[m n o+1]);   
    wdv_dz    = InterpVtoV(wdvdz_atw,m,n,o,'wtov'); 
    
 
    
    %'w'
    [udwdx_atu,vdwdy_atv,wdw_dz] = MimeticUpwindD3D(Wu,Wv,Ww,Uu,Vv,Ww,Guf,Gub,Gvf,Gvb,Gwf,Gwb,m,n,o);
  
    udwdx_atu = reshape(full(udwdx_atu),[m+1 n o]);   
    udw_dx    = InterpVtoV(udwdx_atu,m,n,o,'utow'); 
    


    vdwdy_atv = reshape(full(vdwdy_atv),[m n+1 o]);   
    vdw_dy    = InterpVtoV(vdwdy_atv,m,n,o,'vtow'); 
    
  
    vdu_dy(:,:,1) = vdu_dy(:,:,2);
    vdu_dy(1,:,:) = vdu_dy(2,:,:);
    wdu_dz(1,:,:) = wdu_dz(2,:,:);  %maybe?
    udv_dx(:,:,1) = udv_dx(:,:,2);
    wdv_dz(1,:,:) = wdv_dz(2,:,:);
    wdv_dz(:,1,:) = wdv_dz(:,2,:);
    udw_dx(:,:,1) = udw_dx(:,:,2);     
     vdw_dy(:,:,1) = vdw_dy(:,:,2);
    
    %sl = 25
       
    %RHSu = -(wdu_dz);
    RHSu = -(udu_dx(:)+vdu_dy(:)+wdu_dz(:)); %-sgs_u(:)); %+ DivSgsx; 
    RHSv = -(udv_dx(:)+vdv_dy(:)+wdv_dz(:)); %+ DivSgsy;.
    RHSw = -(udw_dx(:)+vdw_dy(:)+wdw_dz(:)); %-sgs_w(:)); %+ DivSgsx;
%     

%     R = [Uu; Vv; Ww];
%     
%     Uvec_star = D'*I1'*R;
%     
%     RHSu = Uvec_star(1:u_count);
%     RHSv = Uvec_star(u_count+1:v_count);
%     RHSw = Uvec_star(v_count+1:end);
    
    u_star = SSPRK101D(Uu,RHSu,dt);
    v_star = SSPRK101D(Vv,RHSv,dt);        
    %v_star = SSPRK101D(Vv,RHSv-Bb,dt);        
    w_star = SSPRK101D(Ww,RHSw,dt);
    %w_star = SSPRK101D(Ww,RHSw-Bb,dt);
    
%     u_star = Uu + dt*RHSu;
%     v_star = Vv + dt*RHSv;
%     w_star = Ww + dt*(RHSw) ; %-Bb); 
  
     
     u_star = reshape(full(u_star),[m+1 n o]); 
     v_star = reshape(full(v_star),[m n+1 o]); 
     w_star = reshape(full(w_star),[m n o+1]); 

%      figure(325)
%      subplot(3,1,1)
%      pcolor(u_star(:,:,sl));
%      title('u')
%      shading interp
%      subplot(3,1,2)
%      pcolor(v_star(:,:,sl));
%      title('v')
%      shading interp
%      subplot(3,1,3)
%      pcolor(w_star(:,:,sl));
%      title('w')
%      shading interp
     
%      figure(300)
%      pcolor(w_star(:,:,10)');
%      shading interp
     
     
%
  %  u_star = permute(reshape(full(u_star),[m+1 n o]),[2,1,3]); 
  %  v_star = permute(reshape(full(v_star),[m n+1 o]),[2,1,3]); 
  %  w_star = permute(reshape(full(w_star),[m n o+1]),[2,1,3]); 

    
   %[u_star,v_star,w_star,~] = applyboundaries3Dnewgrid(u_star,v_star,w_star,p);
    %[u_star,v_star,w_star,~] = applyboundaries3D(u_star,v_star,w_star,p);
    [u_star,v_star,w_star] = applyboundaries3Dbeam(u_star,v_star,w_star);
% 
    
%     for k=1:o
%         for j=1:n
%             u_star(1,j,k) = u_max ; %u_max*(1-Y_in(1,j,k)); %Real forcing
%             %vatv(:,j,k) = u_max*j/KMax;  %different tryout
%         end
%     end

    u_star(1,:,:) = forcing;
   
%     fig = figure(5);
%     set(gcf,'units','inches','position',[10 10 10 10])
%     u_v = squeeze(u_star(:,:,3));
%     subplot(2,1,1)
%     %contourf(u_v',15);
%     pcolor(u_v');
%     title('FRONT')
%     axis equal
%     axis tight
%     colorbar;
%     shading flat;
%     subplot(2,1,2)
%     u_v = squeeze(u_star(:,25,:));
%     pcolor(u_v');
%     title('TOP')
%     axis equal
%     axis tight
%     colorbar;
%     shading flat;
     
%     Calculate pressure:
    Uu = u_star(:);
    Vv = v_star(:);
    Ww = w_star(:);
    
%    Uu = reshape(permute(u_star, [2, 1, 3]), [], 1);
%    Vv = reshape(permute(v_star, [2, 1, 3]), [], 1);
%    Ww = reshape(permute(w_star, [2, 1, 3]), [], 1);
    
%       
    R = [Uu; Vv; Ww];
    Pp = L\(D*R*rho_ref/dt);
    %Pp = L\(D*R/dt);  
    %Pp = L\(D*R);  
    
    p = reshape(Pp,[m+2 n+2 o+2]);  
    %p = permute(reshape(full(Pp),[m+2 n+2 o+2]),[1,2,3]);
    %p = permute(reshape(full(Pp),[m+2 n+2 o+2]),[2,1,3]);
    
    [~,~,~,p] = applyboundaries3Dstar(u_star,v_star,w_star,p);
    Pp = p(:);

%    
%     figure(888)
%     subplot(2,1,1)
%     u_3 = squeeze(p(:,:,slice));
%     %pcolor(squeeze(Xs(:, :, slice)), squeeze(Zs(:, :, slice)),u_3);
%     pcolor(u_3);
%     shading interp
%     subplot(2,1,2)
%     u_3 = squeeze(p(:,slice,:));
%     pcolor(u_3);
%     shading interp
% 

    %Pp = reshape(permute(p, [1, 2, 3]), [], 1);
    %Pp = reshape(permute(p, [2, 1, 3]), [], 1);

    
    
    %Derivate pressure:
    Gpp = Gp*[ Pp ];
    %Gpp = Gg*[ Pp ];
    dpdx = Gpp(1:u_count);
    dpdy = Gpp(u_count+1:v_count);
    dpdz = Gpp(v_count+1:end);
%
    
%    
    %Correct velocities:
    Uu = Uu + dpdx;                                     
    Vv = Vv + dpdy;
    Ww = Ww + dpdz;

    
%     Uu = Uu - dt*dpdx;                                     
%     Vv = Vv - dt*dpdy;
%     Ww = Ww - dt*dpdz;
      
    uatu = reshape(Uu,[m+1 n o]);
    vatv = reshape(Vv,[m n+1 o]);
    watw = reshape(Ww,[m n o+1]);     
     
%   uatu = permute(reshape(full(dpdx),[m+1 n o]),[2,1,3]); 
%   vatv = permute(reshape(full(dpdy),[m n+1 o]),[2,1,3]); 
%   watw = permute(reshape(full(dpdz),[m n o+1]),[2,1,3]); 
%     
%    sl = 1;  
%    
%     uatu2 = permute(reshape(uatu, m+1, n, o), [2, 1, 3]);
%     vatv2 = permute(reshape(vatv, m, n+1, o), [2, 1, 3]);
%     watw2 = permute(reshape(watw, m, n, o+1), [2, 1, 3]);
% 
% 
%      figure(300)
%      subplot(3,1,1)
%      pcolor(squeeze(Xu(:,:,sl)),squeeze(Yu(:,:,sl)),squeeze(uatu2(:,:,sl)));
%      %pcolor(squeeze(uatu(:,:,sl)));
%      title('u')
%      shading interp
%      subplot(3,1,2)
%      pcolor(squeeze(Xv(:,:,sl)),squeeze(Yv(:,:,sl)),squeeze(vatv2(:,:,sl)));
% %     pcolor(vatv(:,:,sl));
%      title('v')
%      shading interp
%      subplot(3,1,3)
%      pcolor(squeeze(Xw(:,:,sl)),squeeze(Yw(:,:,sl)),squeeze(watw2(:,:,sl)));
% %     pcolor(watw(:,:,sl));
%      title('w')
%      shading interp
%      
%

    
    
      
    %[uatu,vatv,watw,~] = applyboundaries3Dnewgrid(uatu,vatv,watw,p);
    [uatu,vatv,watw] = applyboundaries3Dbeam(uatu,vatv,watw);
    
%    plot_velocities(uatu,vatv,watw,X,Z,IMax)

%   

%     for k=1:o
%         for j=1:n
%             uatu(1,j,k) = u_max*j/KMax; %Real forcing
%             %vatv(:,j,k) = u_max*j/KMax;  %different tryout
%         end
%     end
    
    %uatu(2,:,:) = forcing;

 
%      figure(200) 
%     set(gcf,'units','inches','position',[10 10 8 8])
%     subplot(2,1,1)
%     pcolor(squeeze(vatv(:, 3:end, sl)*UStar))  %visualize X-Y plane
%     %contourf(squeeze(Xu(:, 3:end, sl)), squeeze(Yu(:, 3:end, sl)), squeeze(uatu2(:, 3:end, sz)*UStar), nlevs)  %visualize X-Y plane
%     xlabel('x'); ylabel('y')
%     colorbar;   
%     shading interp
%     title(['U-velocity X-Y plane (front)'])
%     subplot(2,1,2)
%     %pcolor(squeeze(Xw(sl,3:end,:)), squeeze(Zw(sl,3:end,:)), squeeze(watw2(sl,3:end,:)))  %visualize X-Y plane
%     pcolor(squeeze(vatv(sz,3:end,:)*UStar))  %visualize X-Y plane
%     %contourf(squeeze(Xu(2,3:end,:)), squeeze(Zu(2,3:end,:)), squeeze(uatu2(2,3:end,:)*UStar),nlevs)  %visualize X-Y plane
%     xlabel('x'); ylabel('z')
%     colorbar
%     title(['U-velocity X-Z plane (top)'])
%     shading interp
%
     
   % [sgs_u,sgs_v,sgs_w] = CalcSGS(G,D,uatu,vatv,watw,m,n,o,Dx,Dy,Dz);



%         figure(326)
%         subplot(3,1,1)
%         pcolor(squeeze(sgs_u(:,:,3)))
%         shading interp
%         colorbar
%         subplot(3,1,2)
%         pcolor(squeeze(sgs_v(:,:,3)))
%         shading interp
%         colorbar
%         subplot(3,1,3)
%         pcolor(squeeze(sgs_w(:,:,3)))
%         shading interp
%         colorbar
     
     Uu = uatu(:);
     Vv = vatv(:);
     Ww = watw(:);

%    Uu = reshape(permute(uatu, [2, 1, 3]), [], 1);
%    Vv = reshape(permute(vatv, [2, 1, 3]), [], 1);
%    Ww = reshape(permute(watw, [2, 1, 3]), [], 1);
        
    vatu = InterpVtoV(vatv,m,n,o,'vtou'); 
    watu = InterpVtoV(watw,m,n,o,'wtou');

    Vu = vatu(:);
    Wu = watu(:);
    
    uatv = InterpVtoV(uatu,m,n,o,'utov'); 
    watv = InterpVtoV(watw,m,n,o,'wtov');
    
    Uv = uatv(:);
    Wv = watv(:);
    
    uatw = InterpVtoV(uatu,m,n,o,'utow'); 
    vatw = InterpVtoV(vatv,m,n,o,'vtow');
        
    Uw = uatw(:);   
    Vw = vatw(:);   
 
     nondim_time = (time+1)*dt;
     current_time = (time+1)*dt*LStar/UStar
    
%     ohalf = floor(o/2);
%      
%     fig = figure(jj+2);
%     set(gcf,'units','inches','position',[10 10 width height])
%     u_v = squeeze(batv(:,slice,:));
%     subplot(1,2,1)
%     %contourf(u_v',15);
%     pcolor(u_v');
%     axis equal
%     axis tight
%     colorbar;
%     shading flat;
%     subplot(1,2,2)
%     u_v = squeeze(batv(:,:,ohalf));
%     pcolor(u_v');
%     axis equal
%     axis tight
%     colorbar;
%     shading flat;
     

%Visualization: 
 if(mod(nondim_time,print_freq) <= 1e-10 )
    
    width = 8;
    height = 7;
    slice = 3;
    %st = 20;
    %nd = 110;
    st = 1;
    nd = IMax-1;
    
    sz = 2;
    sl = floor(JMax/2);
    
    uatu2 = permute(reshape(uatu, m+1, n, o), [2, 1, 3]);
    vatu2 = permute(reshape(vatu, m+1, n, o), [2, 1, 3]);
    watu2 = permute(reshape(watu, m+1, n, o), [2, 1, 3]);

    velocities_save = [velocities_save; {uatu2(:,:,:)}];
    vvelocities_save = [vvelocities_save; {vatu2(:,:,:)}];
    wvelocities_save = [wvelocities_save; {watu2(:,:,:)}];

%     nlevs = 20;
% 
%     figure(200) 
%     %set(gcf,'units','inches','position',[10 10 8 8])
%     subplot(2,1,1)
%     pcolor(squeeze(Xu(:, 3:end, sl)), squeeze(Yu(:, 3:end, sl)), squeeze(uatu2(:, 3:end, sl)*UStar))  %visualize X-Y plane
%     %contourf(squeeze(Xu(:, 3:end, sl)), squeeze(Yu(:, 3:end, sl)), squeeze(uatu2(:, 3:end, sz)*UStar), nlevs)  %visualize X-Y plane
%     xlabel('x'); ylabel('y')
%     colorbar;   
%     shading interp
%     title(['U-velocity X-Y plane (front)'])
%     subplot(2,1,2)
%     %pcolor(squeeze(Xw(sl,3:end,:)), squeeze(Zw(sl,3:end,:)), squeeze(watw2(sl,3:end,:)))  %visualize X-Y plane
%     pcolor(squeeze(Xu(sz,3:end,:)), squeeze(Zu(sz,3:end,:)), squeeze(uatu2(sz,3:end,:)*UStar))  %visualize X-Y plane
%     %contourf(squeeze(Xu(2,3:end,:)), squeeze(Zu(2,3:end,:)), squeeze(uatu2(2,3:end,:)*UStar),nlevs)  %visualize X-Y plane
%     xlabel('x'); ylabel('z')
%     colorbar
%     title(['U-velocity X-Z plane (top)'])
%     shading interp
%     %drawnow;
%     hold off

    %save('curv_seamount_40000s_101-10-38.mat')

 end
     
end

return
%
%%
%save('curv_seamount_3500s_30-6-50')
%load('sgs_curv_seamount_40000s_50-6-50.mat')
load('curv_seamount_40000s_101-10-38.mat')
%%
% nn = 1200;
% 
% vel_n = velocities_save{n,1};
% %
% figure(876)
% u_3 = squeeze(vel_n(:,:,slice));
% pcolor(u_3');
% shading interp


%

%fList = matlab.codetools.requiredFilesAndProducts('make_beam_grid_gccommolev2.m')

%save('seamount_v4_20000s_100-6-100.mat')

%

%load('seamount_v4_20000s.mat')
%
slice = floor(JMax/2)-1;
%
nlevs = 15;
[N,M] = size(velocities_save);

v = VideoWriter ('curv_seamount_u01_15000s-101x6x50_bwr.avi');
v.FrameRate = 3;
open (v);  
II = 0;

i = N
%
for i = 1:60

vel_n = velocities_save{i,1};

figure(876)
set(gcf,'units','inches','position',[10 10 15 8])
subplot(2,1,1)
u_3 = squeeze(vel_n(:,:,slice)*UStar);
pcolor(squeeze(Xu(:, :, slice)), squeeze(Yu(:, :, slice)),u_3);
%contourf(squeeze(Xu(:, :, slice)*LStar), squeeze(Yu(:, :, slice)*LStar),u_3,nlevs);
axis equal
colormap bluewhitered
hcb=colorbar;
%caxis([-0.01 0.01])
caxis([0 u_max])
title(hcb,'u [m/s]')
title(['U-velocity X-Z plane (front) - Y=5 - [X,Y,Z]=[101x7x50] time: ' num2str(i*print_freq*LStar/UStar) 's'])
xlabel('x [m]');ylabel('z [m]');
shading interp
%shading flat
%writeVideo (v, frame);  
hold off;
subplot(2,1,2)
u_3 = squeeze(vel_n(:,:,slice)*UStar);
pcolor(squeeze(Xu(:, :, slice)), squeeze(Yu(:, :, slice)),u_3);
%contourf(squeeze(Xu(:, :, slice)*LStar), squeeze(Yu(:, :, slice)*LStar),u_3,nlevs);
axis equal
colormap bluewhitered
hcb=colorbar;
%caxis([-0.01 0.01])
caxis([0 u_max])
title(hcb,'u [m/s]')
%title(['U-velocity X-Z plane (front) - Y=5 - [X,Y,Z]=[38x7x50] time: ' num2str(i*print_freq*LStar/UStar) 's'])
xlabel('x [m]');ylabel('z [m]');
%shading interp
%shading flat

hold off;



drawnow;
frame = getframe (gcf);
writeVideo (v, frame);  
end  
%
close (v)
%%

close all
[N,M] = size(velocities_save)

tempc = 0.001;
%tempc = 0.005;
tempc2 = 14;

slice = 1;

temp_n = velocities_save{1,1};

p = patch(isosurface(Ux(:,:,1:end-1),Uy(:,:,1:end-1),Uz(:,:,1:end-1),temp_n,tempc)); 
p2 = patch(isosurface(Ux(:,:,1:end-1),Uy(:,:,1:end-1),Uz(:,:,1:end-1),temp_n,tempc2)); 

%v = VideoWriter ('Seamount_101x10x38.avi');
%open(v);  

i=1

for  i = 1:1:24
    
    %p = patch(isosurface(Ux(:,:,1:end-1),Uy(:,:,1:end-1),Uz(:,:,1:end-1),temp_n,tempc)); 

    temp_n = velocities_save{i,1};

    fig = figure(77);
    set(p2,'Visible','Off')
    set(p,'Visible','Off')
    set(gcf,'units','inches','position',[10 10 10 20])
    %subplot(2,1,1)
    
    figure;
[sx sy sz] = meshgrid(xmin,ymin:10:ymax,zmin:2:zmax);
h = streamline(Ux(:,:,1:end-1),Uy(:,:,1:end-1),Uz(:,:,1:end-1),u,v,w,sx,sy,sz);
set(h,'LineWidth',1,'Color',purple);
view([-40 50]);
axis on;
grid off;
box on;
    
    
    p = patch(isosurface(Ux(:,:,1:end-1),Uy(:,:,1:end-1),Uz(:,:,1:end-1),temp_n,tempc)); 
    p.FaceColor = 'blue';
    %p.EdgeColor = 'none';
    %hold on 
    %p2 = patch(isosurface(Ux(:,:,1:end-1),Uy(:,:,1:end-1),Uz(:,:,1:end-1),temp_n,tempc2)); 
    %p2.FaceColor = 'cyan';
    %p2.EdgeColor = 'none';
    %hold off
    set(gca, 'XLim',[a b]);
    view(-31, 25)
    camlight('left') 
    %axis equal ; 
    axis tight ;
    %title(['T [C] =' num2str(tempc*T0) ' at ' num2str(i*print_freq*LStar/UStar) 's'])
    hLg = legend(['T[C]=' num2str(tempc)],['T[C]=' num2str(tempc2)],'Location','southeast');
    hLg.Box='on';
    hLg.EdgeColor='k';
    %subplot(2,1,2)
    %u_3 = squeeze(temp_n(:,:,slice));
    %pcolor(u_3');
    %axis equal
    %axis tight
    %title(['T  max=' num2str(max(max(max(temp_n)))) ' at ' num2str(i*print_freq*LStar/UStar) 's' ])
    %colorbar;
    %shading interp;
    
%    frame = getframe (gcf);
%    writeVideo (v, frame); 

        %pause(.1)
    
end

%close(v);

%%

nlevs = 15;
[N,M] = size(velocities_save);

figure(876)
%set(gcf,'units','inches','position',[10 10 20 15])
subplot(2,2,1)
i = floor(N/4);
vel_n = velocities_save{i,1};
u_3 = squeeze(vel_n(:,:,slice)*UStar);
%pcolor(squeeze(Xu(:, :, slice)), squeeze(Yu(:, :, slice)),u_3);
contourf(squeeze(Xu(:, :, slice)*LStar), squeeze(Yu(:, :, slice)*LStar),u_3,nlevs);
axis equal
colormap bluewhitered
%hcb=colorbar;
%caxis([-0.01 0.01])
caxis([0 u_max])
%title(hcb,'u [m/s]')
title(['U-velocity Y=5 - [X,Y,Z]=[101x7x50] time: ' num2str(i*print_freq*LStar/UStar) 's'])
xlabel('x [m]');ylabel('z [m]');
subplot(2,2,2)
i = floor(N/2);
vel_n = velocities_save{i,1};
u_3 = squeeze(vel_n(:,:,slice)*UStar);
%pcolor(squeeze(Xu(:, :, slice)), squeeze(Yu(:, :, slice)),u_3);
contourf(squeeze(Xu(:, :, slice)*LStar), squeeze(Yu(:, :, slice)*LStar),u_3,nlevs);
axis equal
colormap bluewhitered
hcb=colorbar;
%caxis([-0.01 0.01])
caxis([0 u_max])
title(hcb,'u [m/s]')
title(['U-velocity Y=5 - [X,Y,Z]=[101x7x50] time: ' num2str(i*print_freq*LStar/UStar) 's'])
xlabel('x [m]');ylabel('z [m]');
subplot(2,2,3)
i = floor(N/2)+floor(N/4);
vel_n = velocities_save{i,1};
u_3 = squeeze(vel_n(:,:,slice)*UStar);
%pcolor(squeeze(Xu(:, :, slice)), squeeze(Yu(:, :, slice)),u_3);
contourf(squeeze(Xu(:, :, slice)*LStar), squeeze(Yu(:, :, slice)*LStar),u_3,nlevs);
axis equal
colormap bluewhitered
%hcb=colorbar;
%caxis([-0.01 0.01])
caxis([0 u_max])
%title(hcb,'u [m/s]')
title(['U-velocity Y=5 - [X,Y,Z]=[101x7x50] time: ' num2str(i*print_freq*LStar/UStar) 's'])
xlabel('x [m]');ylabel('z [m]');
subplot(2,2,4)
i = N;
vel_n = velocities_save{i,1};
u_3 = squeeze(vel_n(:,:,slice)*UStar);
%pcolor(squeeze(Xu(:, :, slice)), squeeze(Yu(:, :, slice)),u_3);
contourf(squeeze(Xu(:, :, slice)*LStar), squeeze(Yu(:, :, slice)*LStar),u_3,nlevs);
axis equal
colormap bluewhitered
%hcb=colorbar;
%caxis([-0.01 0.01])
caxis([0 u_max])
%title(hcb,'u [m/s]')
title(['U-velocity Y=5 - [X,Y,Z]=[101x7x50] time: ' num2str(i*print_freq*LStar/UStar) 's'])
xlabel('x [m]');ylabel('z [m]');





