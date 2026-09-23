clc;clear all; close all; format short
%Nodes:
%Horizontal x, horizontal y, vertical z.
IMax = 401; JMax = 6; KMax = 101; 
%IMax = 201; JMax = 6; KMax = 101; 
%IMax = 101; JMax = 6; KMax = 50; 

width = 16;
height = 7;

%xl = 0.4  ; %0.8m length 
xl = 0.8  ;
yl = 0.1 ;%0.01m depth
%zl = 0.685   ; %0.1m height
zl = 0.0656   ; %0.1m height
%zl = 0.0502;
%zl = 0.0427;
LStar = 1.0; 
%UStar = 1.0; 
UStar = 0.8227 ;
%UStar = 1.5 ;
%UStar = 0.86735 ;
%UStar = 0.91474 ;

%nondimensionalization constants:
rho_ref = 1027.d0; 
S_ref   = 35.d0;
T_ref   = 10.d0;
drho_dS = 0.781d0;
drho_dT = -0.1708d0;
rhoStar = 1000.35;
gforce = 9.80665;
gStar = LStar/UStar^2;

p_atm = 1.01325;
mu = 0.00141;
nu = mu/rho_ref;

%Cells:
m = IMax-1;
n = JMax-1; 
o = KMax-1;

a = -xl;   % West
b = xl;   % East
c = -yl;   % South
d = yl;   % North
e = -zl;
f = zl;

Dx = (b-a)/m ;
Dy = (d-c)/n ;
Dz = (f-e)/o ;

xgrid = [a a+Dx/2 : Dx : b-Dx/2 b];
ygrid = [c c+Dy/2 : Dy : d-Dy/2 d];
zgrid = [e e+Dz/2 : Dz : f-Dz/2 f];
[Xls, Yls, Zls] = meshgrid(xgrid,ygrid,zgrid);
[X, Y, Z] = meshgrid(a:Dx:b, c:Dy:d, e:Dz:f);
Xls = permute(Xls, [2 1 3]);
Yls = permute(Yls, [2 1 3]);
Zls = permute(Zls, [2 1 3]);

X = permute(X, [2 1 3]);
Y = permute(Y, [2 1 3]);
Z = permute(Z, [2 1 3]);
%

xlength = linspace(-xl-Dx,xl+Dx,m+2);
ylength = linspace(-yl,yl,n+1);
zlength = linspace(-zl,zl,o+1);

% % %Special GHOSTED extra dimension just to interpolate to staggered grids
[Yg,Xg,Zg] = ndgrid(ylength,xlength,zlength);
%[Xg, Yg, Zg] = meshgrid(a-Dx:Dx:b, c:Dy:d, e:Dz:f);

% 
Xg = permute(Xg, [2 1 3]);
Yg = permute(Yg, [2 1 3]);
Zg = permute(Zg, [2 1 3]);        
%             

Xu = (Xg(1:end-1,1:n,1:o)+Xg(2:end,1:n,1:o))*0.50;
Yu = (Yg(1:end-1,1:n,1:o)+Yg(2:end,1:n,1:o))*0.50;
Zu = (Zg(1:end-1,1:n,1:o)+Zg(2:end,1:n,1:o))*0.50;

xlength = linspace(-xl,xl,m+1);
ylength = linspace(-yl-Dy,yl+Dy,n+2);
zlength = linspace(-zl,zl,o+1);

% %Special GHOSTED extra dimension just to interpolate to staggered grids
[Yg,Xg,Zg] = ndgrid(ylength,xlength,zlength);
%[Xg, Yg, Zg] = meshgrid(a:Dx:b, c:Dy:d+Dy, e:Dz:f);

% 
 Xg = permute(Xg, [2 1 3]);
 Yg = permute(Yg, [2 1 3]);
 Zg = permute(Zg, [2 1 3]);        

Xv = (Xg(1:m,1:end-1,1:o)+Xg(1:m,2:end,1:o))*0.50;
Yv = (Yg(1:m,1:end-1,1:o)+Yg(1:m,2:end,1:o))*0.50;
Zv = (Zg(1:m,1:end-1,1:o)+Zg(1:m,2:end,1:o))*0.50;
% 
xlength = linspace(-xl,xl,m+1);
ylength = linspace(-yl,yl,n+1);
zlength = linspace(-zl-Dz,zl-Dz,o+2);

% %Special GHOSTED extra dimension just to interpolate to staggered grids
[Yg,Xg,Zg] = ndgrid(ylength,xlength,zlength);
%[Xg, Yg, Zg] = meshgrid(a:Dx:b, c:Dy:d, e:Dz:f+Dz);
% 
Xg = permute(Xg, [2 1 3]);
Yg = permute(Yg, [2 1 3]);
Zg = permute(Zg, [2 1 3]);      

Xw = (Xg(1:m,1:n,1:end-1)+Xg(1:m,1:n,2:end))*0.50;
Yw = (Yg(1:m,1:n,1:end-1)+Yg(1:m,1:n,2:end))*0.50;
Zw = (Zg(1:m,1:n,1:end-1)+Zg(1:m,1:n,2:end))*0.50;

xlength = linspace(-xl,xl,m);
ylength = linspace(-yl,yl,n);
zlength = linspace(-zl,zl,o);

[Xc,Yc,Zc] = meshgrid(xlength,ylength,zlength);

Xc = permute(Xc, [2 1 3]);
Yc = permute(Yc, [2 1 3]);
Zc = permute(Zc, [2 1 3]);

[Xls, Yls, Zls] = ndgrid([a a+Dx/2 : Dx : b-Dx/2 b],...
                      [c c+Dy/2 : Dy : d-Dy/2 d],...
                      [e e+Dz/2 : Dz : f-Dz/2 f]);
%
%
dw = 80.0;
dp = 0.05;

uic =   zeros(m+1,n+1,o+1);
vic =   zeros(m+1,n+1,o+1);
pic =   zeros(m+1,n+1,o+1);
Tic =   zeros(m,n,o);
Sic =   zeros(m+1,n+1,o+1);
T_new = zeros(m+1,n+1,o+1);

p = zeros(m+2,n+2,o+2);
u = zeros(m+2,n+2,o+2);
v = zeros(m+2,n+2,o+2);

%SGS scale:

slice = floor(JMax/2);

%length = (Dx*Dy*Dz)^(1.0/3.0)/4.0d0;
length = (Dx^2+Dy^2+Dy^2)^(1.0/3.0);
%length = (Dx*Dy*Dz)^(1.0/2.0);

%Temperature
%TMAX = 16.1324636; TMIN = 10.0;

TMAX = 16.1285; TMIN = 10.0;

%TMAX = 15.2; TMIN = 10.0;

T0 = TMIN;

Tic = TMIN + (TMAX-TMIN)*(0.5-0.5*erf(Xls/0.01)); 
%Tic = TMIN + (TMAX-TMIN)*(1-erf((Xls)/0.0000001))/2; 

Tic = Tic;

% figure(11);
% u_v = squeeze(Tic(:,slice,:));
% pcolor(u_v');
% title('T')
% colorbar;
% shading interp;


% fig = figure(78);
% p = patch(isosurface(Xls,Yls,Zls,Tic*T0,15)); 
% p.FaceColor = 'white';
% p.EdgeColor = 'none';
% set(gca, 'XLim',[a b]);
% view(-31, 25)
% camlight('left') 
%title(['T [C] =' num2str(tempc*T0) ' at ' num2str(fr_1) 's'])

%
T = Tic;

%Density
%dens = rho_ref + drho_dT*(T*T0-T_ref);  %EOS        
%dens = dens/rhoStar;
alpha = 1958.0e-7 ;
%dens = rho_ref + drho_dT*(T*T0-T_ref);  %EOS        
dens = rho_ref - rho_ref*alpha*(T-T_ref);  %EOS        

DMax = max(max(max(dens)));
DMin = min(min(min(dens)));
rho_P0 = (DMax+DMin)/2.0;

%deltarho = dens(floor(n/2),1) - dens(floor(n/2),end)
r_g = gforce*(DMax-DMin)/rho_ref %reduced gravity 

u_b = sqrt(r_g*zl)

zl2 = (0.0224)^2/r_g

mu = 0.00141;
nu = mu/rho_ref;
Re = 710;
Gr = ((u_b*zl)/nu)^2


Fr = Re/(1.1*sqrt(Gr))

Gr = 1.25e6;
Fr = Re/(1.1*sqrt(Gr))
Fr_an = 1/sqrt(2)
zl3 = nu*sqrt(Gr)/0.0224
Gr = (Re/(1.1*Fr_an))^2
zl4 = nu*sqrt(Gr)/0.0224
%
%%
%
%buoy = zeros(m+2,n+2,o+2);
buoy = gforce*(dens-rho_ref)/rho_ref;
%
% figure(6);
% u_v = squeeze(dens(:,slice,:));
% pcolor(u_v');
% title('rho')
% colorbar;f
% shading interp;

intalg = 'linear';

%Buoyancy
Fbuoy = griddedInterpolant(Xls, Yls, Zls, buoy, intalg, intalg); 
batv = Fbuoy(Xw,Yw,Zw);

size(batv)

%buoy = gforce*gStar*(datw-rho_P0)./rho_P0;
%
% figure(6);
% u_v = squeeze(buoy(:,slice,:));
% pcolor(u_v');
% title('(rho-rho_0)/rho_0')
% colorbar;
% shading interp;

%  Creating operators
k = 2; % Mimetic order of accuracy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dxc = (b-a)/m ;
Dyc = (d-c)/n ;
Dzc = (f-e)/o ;

D =  div3D(k,m,Dxc,n,Dyc,o,Dzc)*LStar; 
G =  grad3D(k,m,Dxc,n,Dyc,o,Dzc)*LStar;

% D =  div3DCurv(k,X,Y,Z)*LStar; 
% G =  grad3DCurv(k,X,Y,Z)*LStar;

%D = div2DCurv(k,X,Y);
%G = grad2DCurv(k,X,Y);

L = D*G;  
%L = lap3D(k,m,Dxc,n,Dyc,o,Dzc);

L = L + LStar^2*robinBC3D(k,m,Dxc,n,Dyc,o,Dzc,0,1); %Neumann boundaries
%L = L + robinBC2D(k,m,Dxc,n,Dyc,1,0); %Dirichlet boundaries

% Preparing t=0 arrays

uic =   zeros(m+1,n+1,o+1);
vic =   zeros(m+1,n+1,o+1);
wic =   zeros(m+1,n+1,o+1);

Fuic = griddedInterpolant(X, Y, Z, uic,intalg,intalg); 
Fvic = griddedInterpolant(X, Y, Z, vic,intalg,intalg); 
Fwic = griddedInterpolant(X, Y, Z, wic,intalg,intalg); 
Fpic = griddedInterpolant(X, Y, Z, pic,intalg,intalg); 
Ftic = griddedInterpolant(Xls, Yls, Zls, T,intalg,intalg); 

uatu = Fuic(Xu,Yu,Zu);
vatu = Fvic(Xu,Yu,Zu);
watu = Fwic(Xu,Yu,Zu);
patu = Fpic(Xu,Yu,Zu);

uatv = Fuic(Xv,Yv,Zv);
vatv = Fvic(Xv,Yv,Zv);
watv = Fwic(Xv,Yv,Zv);
patv = Fpic(Xv,Yv,Zv);

uatw = Fuic(Xw,Yw,Zw);
vatw = Fvic(Xw,Yw,Zw);
watw = Fwic(Xw,Yw,Zw);
patw = Fpic(Xw,Yw,Zw);

u  = Fuic(Xls,Yls,Zls);
v  = Fvic(Xls,Yls,Zls);
w  = Fwic(Xls,Yls,Zls);
p  = Fpic(Xls,Yls,Zls);
Tp = Ftic(Xls,Yls,Zls);
%densp = Fdens(Xls,Yls,Zls);

T = Tic;

% Reshaping

Wv = watv(:);
Vv = vatv(:);
Uv = uatv(:);
Pv = patv(:);

Wu = watu(:);
Vu = vatu(:);
Uu = uatu(:);
Pu = patu(:);

Ww = watw(:);
Vw = vatw(:);
Uw = uatw(:);
Pw = patw(:);

Wls = w(:);
Vls = v(:);
Uls = u(:);
Pp = p(:);

width = 16;
height = 7;

% cs = 1447.0;
% beta = 1/cs;
% Re = 750;
% %Re = 100;
% IRe = 1.0/Re;
% nu = 0.0014;

%[uatu,vatv,watv,p] = applyboundaries3D(uatu,vatv,watv,p); 

Fuatu = griddedInterpolant(Xu, Yu, Zu, uatu,intalg,intalg); 
Fvatv = griddedInterpolant(Xv, Yv, Zv, vatv,intalg,intalg); 
Fwatw = griddedInterpolant(Xw, Yw, Zw, watw,intalg,intalg); 

u = Fuatu(Xls,Yls,Zls);
v = Fvatv(Xls,Yls,Zls);
uatc = Fuatu(Xc,Yc,Zc);
vatc = Fvatv(Xc,Yc,Zc);

vatu = Fvatv(Xu,Yu,Zu);
watu = Fwatw(Xu,Yu,Zu);

uatv = Fuatu(Xv,Yv,Zv);
watv = Fwatw(Xv,Yv,Zv);

uatw = Fuatu(Xw,Yw,Zw);
vatw = Fvatv(Xw,Yw,Zw);

%SGS calculation:

%Calc velocity derivatives:
dudx = zeros(m+1,n,o);
dudy = zeros(m+1,n,o);
dudz = zeros(m+1,n,o);

dvdx = zeros(m,n+1,o);
dvdy = zeros(m,n+1,o);
dvdz = zeros(m,n+1,o);

dwdx = zeros(m,n,o+1);
dwdy = zeros(m,n,o+1);
dwdz = zeros(m,n,o+1);

%dvdxy:

Uls = u(:); 
Vls = v(:); 
Wls = w(:); 

jj = 2;
time = 0;

%Time discretization (needed for D and G):
tf = 64.0; 
%tf = 0.10; 
tf = tf*(UStar/LStar);
%dt = 0.05;
dt = 0.01;

%CFL dt limit:
%maxv = UStar;
%dt = (min(Dx,Dy)/maxv)/2;

dt = dt*(UStar/LStar);
nf = tf/dt ;
time_steps = linspace(0,tf,nf);

%print_freq = 50;
print_freq = 0.5;

Gp = -dt/rho_ref*G;

%Gp = G;

velocities_save = [];
velocities_savev = [];
velocities_savew = [];
density_save = [];

RHSu  = zeros(n*(m+1)*o,1);
RHSv  = zeros((n+1)*m*o,1);
RHSw  = zeros(n*m*(o+1),1);
dP_dx = zeros(m+1,n+1,o+1);
dP_dy = zeros(m+1,n+1,o+1);
T_new = zeros(m+1,n+1,o+1);   
    
u_count = (m+1)*n*o;
v_count = u_count + (n+1)*m*o;
    
I0 = interpol3D(m, n, o, 0, 0, 0);
I1 = interpol3D(m, n, o, 1, 1, 1);

DI0 = interpolD3D(m, n, o, 0, 0, 0);
DI1 = interpolD3D(m, n, o, 1, 1, 1);

DI0_left  = DI0(:,1:u_count);
DI0_middle = DI0(:,u_count+1:v_count);
DI0_right = DI0(:,v_count+1:end);

DI1_left  = DI1(:,1:u_count);
DI1_middle = DI1(:,u_count+1:v_count);
DI1_right = DI1(:,v_count+1:end);

temperature_save = [];

timedim = [];

[Gub, Gvb, Gwb] = d3D(m, n, o, Dx, Dy, Dz, 0, 0, 0); 
[Guf, Gvf, Gwf] = d3D(m, n, o, Dx, Dy, Dz, 1, 1, 1); 
T_new = zeros(m+2,n+2,o+2);   
%
timedim = [];
Main;

%%
%save('LE3D-400x6x100-60s-08142021.mat');
%return
%%
%load('LE3D-400x6x100-complete.mat');

%%

Fr_an = 1/sqrt(2);
u_f = [];
%u_b = sqrt(abs(r_g*2*zl/2))
u_b = sqrt(r_g*zl)
%u_b = 0.0224;

start_frame = 1;

for i = start_frame:100 ;%numel(velocities_save)
    u_frame = velocities_save{i,1};
    u_f(i-start_frame+1) = max(max((abs(u_frame(:,:,2)))))*UStar;
    %u_f(i) = max(max((abs(u_frame(:,2,:)))))*UStar;
    %u_f(i) = max(max((abs(u_frame(2,:,:)))))*UStar;
end

Fr1 = u_f./u_b;

mFr = mean(Fr1);

for i = start_frame:100 %numel(velocities_save)
    u_frame = velocities_save{i,1};
    %u_f(i) = max((abs(u_frame(end-1,:,:))))*UStar;
    u_f(i-start_frame+1) = max(max((abs(u_frame(:,:,end-1)))))*UStar;
end

%Fr2 = Fr1;
Fr2 = u_f./u_b;

mFr2 = mean(Fr2);

mFravg = (mFr + mFr2)/2.0;
errf = 100-mFravg/Fr_an*100;

diff = Fr2 - Fr1;


fig = figure(77);
set(gcf,'units','inches','position',[10 10 width height]);
figure(77);
semilogy(Fr1,'+-','LineWidth',2);
hold on
semilogy(Fr2,'LineWidth',2);
%hold on
%semilogy(diff,'LineWidth',2);
yline(Fr_an,'--');
hold off
legend(['Bottom=' num2str(mFr)],['Top=' num2str(mFr2)],'Theoretical','Location','southeast')
%legend(['Bottom=' num2str(mFr)],['Top=' num2str(mFr2)],['Difference'],'Theoretical','Location','southeast')
title(['3D LE ' num2str(IMax) 'x' num2str(JMax) 'x' num2str(KMax) ' cells, Froude number mean = ' num2str(mFravg) ' err(%)=' num2str(errf)])
shading flat;

f_name = ['LE3D_' num2str(IMax) '_dt1_Froude_401x6x101.png'];
%saveas(fig,f_name,'png')
%%

Fr_an = 1/sqrt(2);
u_f = [];
%u_b = sqrt(abs(r_g*2*zl/2))
u_b = sqrt(r_g*zl)
%u_b = 0.0224;

start_frame = 1;

for i = start_frame:92 ;%numel(velocities_save)
    u_frame = velocities_save{i,1};
    u_f(i-start_frame+1) = max(max((abs(u_frame(:,:,2)))))*UStar;
    %u_f(i) = max(max((abs(u_frame(:,2,:)))))*UStar;
    %u_f(i) = max(max((abs(u_frame(2,:,:)))))*UStar;
end

Fr1 = u_f./u_b;

mFr = mean(Fr1);

errf = 100-mFr/Fr_an*100;

fig = figure(77);
set(gcf,'units','inches','position',[10 10 width height]);
figure(77);
semilogy(Fr1,'+-','LineWidth',2);
%hold on
%semilogy(Fr2,'LineWidth',2);
%hold on
%semilogy(diff,'LineWidth',2);
yline(Fr_an,'--');
hold off
%legend(['Bottom=' num2str(mFr)],['Top=' num2str(mFr2)],'Theoretical','Location','southeast')
legend(['Bottom=' num2str(mFr)],'Theoretical','Location','southeast')
%legend(['Bottom=' num2str(mFr)],['Top=' num2str(mFr2)],['Difference'],'Theoretical','Location','southeast')
title(['3D LE ' num2str(IMax) 'x' num2str(JMax) 'x' num2str(KMax) ' cells, Froude number mean = ' num2str(mFravg) ' err(%)=' num2str(errf)])
shading flat;

f_name = ['LE3D_' num2str(IMax) '_dt1_Froude.png'];

%%
% n_levels = 7;
% fig = figure(3);
% set(gcf,'units','inches','position',[10 10 20 20])
% subplot(4,1,1)
% fr_1 = print_freq*20;
% u_frame = temperature_save{fr_1,1};
% u_2 = squeeze(u_frame(2:end-1,2:end-1)*T0);
% [c h]= contour(Xls(2:end-1,2:end-1),Yls(2:end-1,2:end-1),u_2,n_levels);  h.LineColor= 'black'; h.LineWidth = 1.0;
% xlabel('x [m]')
% ylabel('y [m]')
% title(['T [C] contours at ' num2str(fr_1/2 ) 's'])
% axis equal
% axis tight 
% subplot(4,1,2)
% fr_1 = print_freq*60;
% u_frame = temperature_save{fr_1,1};
% u_2 = squeeze(u_frame(2:end-1,2:end-1)*T0);
% [c h]= contour(Xls(2:end-1,2:end-1),Yls(2:end-1,2:end-1),u_2,n_levels);  h.LineColor= 'black'; h.LineWidth = 1.0;
% title(['T [C] contours at ' num2str(fr_1/2 ) 's'])
% xlabel('x [m]')
% ylabel('y [m]')
% axis equal
% axis tight 
% subplot(4,1,3)
% fr_1 = print_freq*80;
% u_frame = temperature_save{fr_1,1};
% u_2 = squeeze(u_frame(2:end-1,2:end-1)*T0);
% [c h]= contour(Xls(2:end-1,2:end-1),Yls(2:end-1,2:end-1),u_2,n_levels);  h.LineColor= 'black'; h.LineWidth = 1.0;
% title(['T [C] contours at ' num2str(fr_1/2 ) 's'])
% xlabel('x [m]')
% ylabel('y [m]')
% axis equal
% axis tight 
% subplot(4,1,4)
% fr_1 = print_freq*160;
% u_frame = temperature_save{fr_1,1};
% u_2 = squeeze(u_frame(2:end-1,2:end-1)*T0);
% [c h]= contour(Xls(2:end-1,2:end-1),Yls(2:end-1,2:end-1),u_2,n_levels);  h.LineColor= 'black'; h.LineWidth = 1.0;
% title(['T [C] contours at ' num2str(fr_1/2 ) 's'])
% xlabel('x [m]')
% ylabel('y [m]')
% axis equal
% axis tight 
% f_name = ['LE3D_fullmimetic_contour_' num2str(JMax) '_dt.1_' num2str(current_time) 's.png'];
% saveas(fig,f_name,'png')
%
%

fin = size(temperature_save,1);

fin34 = floor(fin/2) + floor(fin/4);

fin12 = floor(fin/2) - floor(fin/4);

frame_times = [1 fin12 fin34 fin-1];
fr_1 = frame_times(1);
fr_2 = frame_times(2);
fr_3 = frame_times(3);
fr_4 = frame_times(4);

% xgrid = [a a+Dx/2 : Dx : b-Dx/2 b];
% ygrid = [c c+Dy/2 : Dy : d-Dy/2 d];
% zgrid = [e e+Dz/2 : Dz : f-Dz/2 f];
% [Xls2, Yls2, Zls2] = meshgrid(xgrid,ygrid,zgrid);
% Xls2 = permute(Xls2, [2 1 3]);
% Yls2 = permute(Yls2, [2 1 3]);
% Zls2 = permute(Zls2, [2 1 3]);
% 
% Xls2 = Xls2(2:end-1,2:end-1,2:end-1);
% Yls2 = Yls2(2:end-1,2:end-1,2:end-1);
% Zls2 = Zls2(2:end-1,2:end-1,2:end-1);

tempc = 12;
tempc2 = 14;
%close all

fig = figure(78);
set(gcf,'units','inches','position',[10 10 10 20])
 subplot(4,1,1)
u_frame = temperature_save{fr_1,1};
p = patch(isosurface(Xls,Yls,Zls,u_frame,tempc)); 
p.FaceColor = 'white';
p.EdgeColor = 'none';
hold on 
p = patch(isosurface(Xls,Yls,Zls,u_frame,tempc2)); 
p.FaceColor = 'cyan';
p.EdgeColor = 'none';
hold off
set(gca, 'XLim',[a b]);
view(-31, 25)
camlight('left') 
title(['T [C] =' num2str(tempc*T0) ' at ' num2str(fr_1*dt) 's'])
subplot(4,1,2)
u_frame = temperature_save{fr_2,1};
p = patch(isosurface(Xls,Yls,Zls,u_frame,tempc)); 
p.FaceColor = 'white';
p.EdgeColor = 'none';
hold on 
p = patch(isosurface(Xls,Yls,Zls,u_frame,tempc2)); 
p.FaceColor = 'cyan';
p.EdgeColor = 'none';
hold off
set(gca, 'XLim',[a b]);
view(-31, 25)
camlight('left') 
 title(['T [C] =' num2str(tempc*T0) ' at ' num2str(fr_2*dt) 's'])
subplot(4,1,3)
u_frame = temperature_save{fr_3,1};
p = patch(isosurface(Xls,Yls,Zls,u_frame,tempc)); 
p.FaceColor = 'white';
p.EdgeColor = 'none';
hold on 
p = patch(isosurface(Xls,Yls,Zls,u_frame,tempc2)); 
p.FaceColor = 'cyan';
p.EdgeColor = 'none';
hold off
set(gca, 'XLim',[a b]);
camlight('left') 
view(-31, 25)
%axis equal ; axis tight
title(['T [C] =' num2str(tempc*T0) ' at ' num2str(fr_3*dt) 's'])
subplot(4,1,4)
u_frame = temperature_save{fr_4,1};
p = patch(isosurface(Xls,Yls,Zls,u_frame,tempc)); 
p.FaceColor = 'white';
p.EdgeColor = 'none';
hold on 
p = patch(isosurface(Xls,Yls,Zls,u_frame,tempc2)); 
p.FaceColor = 'cyan';
p.EdgeColor = 'none';
hold off
set(gca, 'XLim',[a b]);
view(-31, 25)
camlight('left') 
%axis equal ; axis tight
title(['T [C] =' num2str(tempc*T0) ' at ' num2str(fr_4*dt) 's'])
hLg = legend(['T[C]=' num2str(tempc)],['T[C]=' num2str(tempc2)],'Location','southeast');
hLg.Box='on';
hLg.EdgeColor='k';
f_name = ['LE3Disosurface_' num2str(IMax) '_dt' num2str(dt) '_' num2str(current_time) 's.png'];
%saveas(fig,f_name,'png')
%%
close all
[N,M] = size(temperature_save)

tempc = 12;
tempc2 = 14;

slice = 1;

temp_n = temperature_save{1,1};

p = patch(isosurface(Xls,Yls,Zls,temp_n,tempc)); 
p2 = patch(isosurface(Xls,Yls,Zls,temp_n,tempc2)); 

v = VideoWriter ('LE3D_401x101.avi');
open(v);  

for  i = 1:1:125
    

    temp_n = temperature_save{i,1};

    fig = figure(77);
    set(p2,'Visible','Off')
    set(p,'Visible','Off')
    set(gcf,'units','inches','position',[10 10 10 20])
    subplot(2,1,1)
    p = patch(isosurface(Xls,Yls,Zls,temp_n,tempc)); 
    p.FaceColor = 'blue';
    %p.EdgeColor = 'none';
    hold on 
    p2 = patch(isosurface(Xls,Yls,Zls,temp_n,tempc2)); 
    p2.FaceColor = 'cyan';
    %p2.EdgeColor = 'none';
    hold off
    set(gca, 'XLim',[a b]);
    view(-31, 25)
    camlight('left') 
    %axis equal ; 
    axis tight ;
    title(['T [C] =' num2str(tempc*T0) ' at ' num2str(i*print_freq*LStar/UStar) 's'])
    hLg = legend(['T[C]=' num2str(tempc)],['T[C]=' num2str(tempc2)],'Location','southeast');
    hLg.Box='on';
    hLg.EdgeColor='k';
    subplot(2,1,2)
    u_3 = squeeze(temp_n(:,slice,:));
    pcolor(u_3');
    axis equal
    axis tight
    title(['T  max=' num2str(max(max(max(temp_n)))) ' at ' num2str(i*print_freq*LStar/UStar) 's' ])
    colorbar;
    shading interp;
    
    frame = getframe (gcf);
    writeVideo (v, frame); 

        %pause(.1)
    
end

close(v);
%%

figure

for i = 1:5
    temp_n = temperature_save{i,1};
    h = isosurface(Xls,Yls,Zls,temp_n,tempc2);
    p = patch(h); 
    view(-31, 25)
    pause(.5)
    set(p,'Visible','Off')
end

%%
tempc = 1.2;
tempc2 = 1.5;
close all
figure(78)
u_frame = temperature_save{10,1};
p = patch(isosurface(Xls,Yls,Zls,u_frame,tempc)); 
p.FaceColor = 'w';
p.EdgeColor = 'none';
p = patch(isosurface(Xls,Yls,Zls,u_frame,tempc2)); 
p.FaceColor = 'c';
p.EdgeColor = 'none';
set(gca, 'XLim',[a b])
hLg = legend(['T[C]=' num2str(tempc)],['T[C]=' num2str(tempc2)],'Location','southeast');
hLg.Box='on';
hLg.EdgeColor='k';
hold off

%ylim([1 IMax])
view(-31, 25)
%axis equal ; axis tight
%title(['T [C] =' num2str(tempc*T0) ' at ' num2str(fr_4) 's'])
%zlim([1 IMax])
%shading interp

 
%axis tight
camlight('left') 
%lighting gouraud
hold off
%%

%Energy Conservation
% Based on get_volume() from the gccom matlab routine
dVol=zeros(m,n,o,'double');

%%
for k=1:o-2
    for j=1:n-1
        for i=1:m-1

            xi(1)=(Xu(i+1,j,k)-Xu(i,j,k));
            xi(2)=(Yu(i+1,j,k)-Yu(i,j,k));
            xi(3)=(Zu(i+1,j,k)-Zu(i,j,k));
            
            xj(1)=(Xv(i,j+1,k)-Xv(i,j,k));
            xj(2)=(Yv(i,j+1,k)-Yv(i,j,k));
            xj(3)=(Zv(i,j+1,k)-Zv(i,j,k));

            
            xk(1)=(Xw(i,j,k+1)-Xw(i,j,k));
            xk(2)=(Yw(i,j,k+1)-Yw(i,j,k));
            xk(3)=(Zw(i,j,k+1)-Zw(i,j,k));            
            
            dVol(i,j,k)= dot(xk,cross(xi,xj));
        end
    end
end

Vol=sum(sum(sum(dVol)));
dVolBar=sum(dVol,3);
%%
%Based in gccom_energy() 
%
[N,M] = size(velocities_save);
zcc=Zw+max(Zw(:));
zcc = zcc(:,:,1:end-1);

for t=1:N
    %t

    u_n = velocities_save{t,1};
    v_n = velocities_savev{t,1};
    w_n = velocities_savew{t,1};
    rho_n = density_save{t,1};
    
    ucc=u_n(1:end-1,:,:); temp_max_u=max(ucc(:));
    vcc=v_n(:,1:end-1,:); temp_max_v=max(vcc(:));
    wcc=w_n(:,:,1:end-1); temp_max_w=max(wcc(:));
    
    dcc= rho_n(2:end-1,2:end-1,2:end-1);
     uBar=sum((dVol.*ucc),3)./dVolBar;
     vBar=sum((dVol.*vcc),3)./dVolBar;
     wBar=sum((dVol.*wcc),3)./dVolBar;

%     uBar=sum((dVol.*ucc),3)./dVol;
%     vBar=sum((dVol.*vcc),3)./dVol;
%     wBar=sum((dVol.*wcc),3)./dVol;


    eKin(t)=sum(sum(sum(((ucc.^2+vcc.^2 +wcc.^2).*dVol))))/(2.0*Vol);
    
    eKinDen(t)=sum(sum(sum(dcc.*(ucc.^2+vcc.^2 +wcc.^2).*dVol)))/(2.0*Vol);
    
    %eKinDen(t)=sum(sum(sum((ucc.^2+vcc.^2 +wcc.^2).*dVol)))/2.0;

    ePont(t)=sum(sum(sum(gforce*dcc.*zcc.*dVol)))/Vol; 
     
    eBar(t)=sum(sum( (uBar.^2+vBar.^2+wBar.^2).*dVolBar ))/(2.0*Vol);
    
    eKineBar(t)=eBar(t)/eKin(t);
  
   % VMax(t-STime+1)=max(max(max( max(ucc,vcc,wcc) )));
   VMax(t,:) = [temp_max_u,temp_max_v,temp_max_w];
   
   VVmax(t,:)= [ max(max(max(ucc-repmat(uBar,[1 1,KMax-1])))), ...
                 max(max(max(vcc-repmat(vBar,[1 1,KMax-1])))), ...
                 max(max(max(wcc-repmat(wBar,[1 1,KMax-1])))) ];    
end

%%
figure(998)
ax=gca
ax.FontSize = 14
subplot(2,1,1)
yyaxis left
%ax.YColor = 'r'
%plot(timedim,ePont,'--','linewidth',2)
plot(ePont,'--','linewidth',2)
ylabel('PE')
yyaxis right
%ax.YColor = 'k'
%plot(timedim,eKinDen,'','linewidth',2)
plot(eKinDen,'','linewidth',2)
ylabel('KE')
subplot(2,1,2)
TT=(ePont+eKinDen)/ePont(1);
%plot(timedim,TT,'k','linewidth',2)
plot(TT,'k','linewidth',2)
ylabel('TE/TE_0','FontSize',14)
%ylim([0.95 1.01])
ylim([0.9995 1.0005])
%ylim([0.99995 1.00005])
%xlim([timedim(1) timedim(end)])
f_name = ['LE3D_energy_conserv' num2str(IMax) '_dt' num2str(dt) '_' num2str(current_time) 's.png'];
%saveas(gca,f_name,'png')

