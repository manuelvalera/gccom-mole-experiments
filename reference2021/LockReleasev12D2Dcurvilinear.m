clc;clear all; close all; format short

IMax = 100; JMax = 200;

width = 16;
height = 7;
% 
% xl = 16  ; %0.8m length 
% yl = 1   ; %0.1m depth
% LStar = 0.05; 
% UStar = 1.0;

% xl = 40.0;  
% yl = 5.0 ; 
% LStar = 0.01;
% UStar = 1.0 ;

% xl = 8  ;    %GCCOM uses these?
% yl = 1   ; 
% LStar = 0.05; 
% UStar = 1.0;
 
xl = 0.8  ; 
yl = 0.0656; 
LStar = 1.0; 
UStar = 0.8277;

% xl = 50  ; 
% yl = 10; 
% LStar = 1.0; 
% UStar = 1.0;

m = JMax-1;
n = IMax-1; 

a = -xl;   % West
b = xl;   % East
c = -yl;   % South
d = yl;   % North

Dx = (b-a)/m ;
Dy = (d-c)/n ;

% xgrid = [a a+Dx/2 : Dx : b-Dx/2 b];
% ygrid = [c c+Dy/2 : Dy : d-Dy/2 d];
% 
% [Xls, Yls] = ndgrid(xgrid,ygrid);

% Xls = Xls';
% Yls = Yls';
% 
% xlength = linspace(a,b,JMax);
% ylength = linspace(c,d,IMax);
% 
% [Y,X] = ndgrid(ylength,xlength);

xgrid = [a a+Dx/2 : Dx : b-Dx/2 b];
ygrid = [c c+Dy/2 : Dy : d-Dy/2 d];
[Xls, Yls] = meshgrid(xgrid,ygrid);
[X, Y] = meshgrid(a:Dx:b, c:Dy:d);

for j=1:JMax
 for i=1:IMax

%      X2(i,j) = X(i,j) + cos(50*pi()*Y(i,j)/10);
%      X2(i,j) = X2(i,j) - xl;
    
     X2(i,j) = X(i,j) + sin(5*pi()*Y(i,j)/10);

 end
 end
 
for j=1:JMax+1
 for i=1:IMax+1
    
%      Xls2(i,j) = Xls(i,j) + cos(50*pi()*Yls(i,j)/10);
%      Xls2(i,j) = Xls2(i,j) - xl;

      Xls2(i,j) = Xls(i,j) + sin(pi()*Yls(i,j)/10);

 end
end
    
X = X2;
Xls = Xls2;
    
%
figure(01)
plot(Xls2,Yls,Xls2.',Yls.')
title gridded, xlabel x, ylabel y
%axis equal


Xu = (X(1:IMax-1,:)+X(2:IMax,:))*0.50;
Yu = (Y(1:IMax-1,:)+Y(2:IMax,:))*0.50;

Xv = (X(:,1:JMax-1)+X(:,2:JMax))*0.50;
Yv = (Y(:,1:JMax-1)+Y(:,2:JMax))*0.50;

%%
%EOS and buoyancy parameters:
rho_ref = 1027.5d0; 
S_ref   = 35.d0;
T_ref   = 10.d0;
drho_dS = 0.781d0;
drho_dT = -0.1708d0;
%drho_dT = -0.1154d0;
rhoStar = 1000.35;
gforce = 9.8;
gStar = LStar/UStar^2;
%1013.25 mbar (101.325 kPa; 29.921 inHg; 760.00 mmHg)
p_atm = 1.01325;
%mu = 3.38e-3; %at 15C
%mu = 4.03e-3; %at 10C
%mu = 1.4;

mu = 0.00141;

nu = mu/rho_ref

%IRedim = mu/(rho_ref*UStar*LStar);

%Redim = 1/IRedim

Cs = 0.22;
Nu = sparse(n,m); 
length = (Dx^2+Dy^2)^(1.0/2.0);

%gStar = gforce*LStar/(UStar^2);

%gStar = gforce*(UStar^2)/LStar;

uatu = sparse(m+1,n)';
vatv = sparse(m,n+1)';

p = sparse(IMax+1,JMax+1);
u = sparse(IMax+1,JMax+1);
v = sparse(IMax+1,JMax+1);

T = sparse(IMax+1,JMax+1);
buoy = sparse(IMax+1,JMax+1);

%TMAX = 16.1324636; TMIN = 10.0;
TMAX = 16.1285; TMIN = 10.0;

T0 = TMIN;

T = TMIN + (TMAX-TMIN)*(1-erf(Xls/0.01))/2;      %Acc. to paper
%T = fliplr(T);
 
close all
T = T/T0;

%Density
%dens = rho_ref + drho_dT*(T*T0-T_ref);  %EOS        
%dens = dens/rhoStar;
%dens = rho_ref + drho_dT*(T-T_ref);  %EOS        

%alpha = 1.664e-4 ;
%dens = rho_ref*( 1 - alpha*(T*T0-T_ref)); %Validation paper

alpha = 1958.0e-7 ;
dens = rho_ref - alpha*rho_ref*(T*T0-T_ref);

DMax = max(max(dens));
DMin = min(min(dens));
rho_P0 = (DMax+DMin)/2.0;

% figure(77);
% imagesc(dens);
% title('dens')
% colorbar;
% shading flat;


deltarho = dens(floor(n/2),1) - dens(floor(n/2),end);
r_g = -gforce*deltarho/rho_ref %reduced gravity 
u_b = sqrt(r_g*yl)
%u_b = 0.0224
Gr = ((u_b*yl)/nu)^2

Gr = 1.5e6;
yl = (Gr*(nu^2)/0.01)^(1/3)
yl2 = (0.0224)^2/0.01


Ri = gforce*LStar*(DMax-DMin)/(rho_P0*UStar*UStar)

%rhoprime = dens-rho_P0;
%buoy = Ri*rhoprime;
%buoy = gforce*gStar*(dens-rho_P0)/(rho_P0);

buoy = gforce*(dens-rhoStar)/(rhoStar);
%buoy = 2*(dens-rhoStar)/(rhoStar);
buoy2 = buoy;

buoy = reshape(buoy.',[],1);
Xlsv = reshape(Xls.',[],1);
Ylsv = reshape(Yls.',[],1);
Fbuoy = scatteredInterpolant(Xlsv, Ylsv, buoy); 
batv = Fbuoy(Xv',Yv')';

%buoy = gforce*gStar*(datv-rho_P0)/(2*rho_P0);
%buoy3 = gStar*(datv-rho_P0)/rho_P0;

%buoy = 0.00024;

% figure(77);
% imagesc(buoy);
% title('buoy')
% colorbar;
% shading flat;

u = reshape(u.',[], 1);
v = reshape(v.',[], 1);

Uls = reshape(u.',[], 1);
Vls = reshape(v.',[], 1);    

Fuatu = scatteredInterpolant(Xlsv, Ylsv, full(u));
uatu = Fuatu(Xu',Yu')';
uatv = Fuatu(Xv',Yv')';

Fvatv = scatteredInterpolant(Xlsv, Ylsv, full(v));
vatv = Fvatv(Xv',Yv')';
vatu = Fvatv(Xu',Yu')';

Uv = reshape(uatv.',[], 1);
Vu = reshape(vatu.',[], 1);  

Re = 750;
%Re = 750/80;
IRe = 1/Re;

% fig = figure(1);
%     set(gcf,'units','inches','position',[10 10 width height])
%     u_v = squeeze(batv(2:end-1,2:end-1)*UStar);
%     subplot(2,2,1)
%     pcolor(u_v);
%     axis equal
%     axis tight
%     title('batv')
%     colorbar;
%     shading interp;
%     subplot(2,2,2)
%     u_3 = squeeze(buoy2(2:end-1,2:end-1)*UStar);
%     pcolor(u_3);
%     axis equal
%     axis tight
%     title('buoy')
%     colorbar;
%     shading interp;
%     subplot(2,2,3)
%     u_2 = squeeze(dens(2:end-1,2:end-1));
%     pcolor(u_2);
%     axis equal
%     axis tight 
%     title('dens')
%     colorbar;
%     shading interp;
%     subplot(2,2,4)
%     u_2 = squeeze(T(2:end-1,2:end-1)*T0);
%     pcolor(u_2);
%     axis equal
%     axis tight 
%     title('T')
%     colorbar;
%     shading interp;

Xuv = reshape(Xu.',[],1);
Yuv = reshape(Yu.',[],1);

Xvv = reshape(Xv.',[],1);
Yvv = reshape(Yv.',[],1);

Vv = reshape(vatv.',[], 1);
Uu = reshape(uatu.',[], 1);

Fuatu = scatteredInterpolant(Xuv, Yuv, Uu); 
Fvatv = scatteredInterpolant(Xvv, Yvv, Vv); 

%u = Fuatu(Xls',Yls')';
%v = Fvatv(Xls',Yls')';
uatv = Fuatu(Xv',Yv')';
vatu = Fvatv(Xu',Yu')';

k = 2; % Mimetic order of accuracy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dxc = Dx;
Dyc = Dy; 

%D =  div2D(k,m,Dxc,n,Dyc)*LStar; 
%G =  grad2D(k,m,Dxc,n,Dyc)*LStar;

D =  div2DCurv(k,X,Y)*LStar; 
G =  grad2DCurv(k,X,Y)*LStar;

L = D*G;  


%L = L + LStar^2*robinBC2D(k,m,Dxc,n,Dyc,1,0); %Neumann boundaries

%L = L + LStar^2*robinBC2D(k,m,Dxc,n,Dyc,0,1); %Neumann boundaries

Robin = robinBC2D(k,m,Dxc,n,Dyc,1,0); %Neumann boundaries
L = L + LStar^2*Robin; 

I0 = interpol2D(m, n, 0, 0);
I1 = interpol2D(m, n, 1, 1);

DI0 = interpolD2D(m, n, 0, 0);

DI0_left  = DI0(:,1:(m+1)*n);
DI0_right = DI0(:,(m+1)*n+1:end);

DI1 = interpolD2D(m, n, 1, 1);

DI1_left  = DI1(:,1:(m+1)*n);
DI1_right = DI1(:,(m+1)*n+1:end);

velocities_save = [];
temperature_save = [];
momentum_save = [];
momentum_savex = [];
pgrad_save = [];
pgrad_savex = [];

buoy_save = [];

%Time discretization :
tf = 5.00;
tf = tf*(UStar/LStar);
dt = 0.1;

dt = dt*(UStar/LStar);
nf = tf/dt ;
time_steps = linspace(0,tf,nf);

print_freq = 0.1;
Pr = 1.0;
jj = 2;
time = 0;

Gp = -dt/rho_ref*G;
%Gp = G;

%[Guf, Gvf] = d2D(m, n, Dx, Dy, 1, 1); 
%[Gub, Gvb] = d2D(m, n, Dx, Dy, 0, 0); 

%
% close all
% fig = figure(100);
% set(gcf,'units','inches','position',[4 10 5 7])
% u_v = squeeze(batv(2:end-1,2:end-1));
% pcolor(u_v);
% axis tight
% title('Bb')
% colorbar;
% shading interp;
%
for time = 1:size(time_steps,2)
    
    RHSu  = sparse(n*(m+1),1);
    RHSv  = sparse((n+1)*m,1);
    T_new = sparse(n+2,m+2);
    D2T = sparse(n+2,m+2);

    Vv_star = reshape(vatv.',[], 1);
    Uu_star = reshape(uatu.',[], 1);

    Bb = reshape(batv.',[], 1);
    
    %u,v derivatives:
    
    %[ududx,vdudy_atv] = MimeticUpwindD2D(Uu,Uv,Uu,Vv,Guf,Gub,Gvf,Gvb,m,n);
    %[udvdx_atu,vdvdy] = MimeticUpwindD2D(Vu,Vv,Uu,Vv,Guf,Gub,Gvf,Gvb,m,n);
    
%    Fdudyatv = scatteredInterpolant(Xvv, Yvv, full(vdudy_atv));
%    vdudy = Fdudyatv(Xu',Yu')'; 
%      
%    Fdvdxatx = scatteredInterpolant(Xuv, Yuv, full(udvdx_atu));
%    udvdx = Fdvdxatx(Xv',Yv')';
    
%     [ududx_ls,vdudy_ls] = MimeticUpwindv3grad(Uls,Uls,Vls,DI0_left,DI1_left,DI0_right,DI1_right,G,m,n);
%     [udvdx_ls,vdvdy_ls] = MimeticUpwindv3grad(Vls,Uls,Vls,DI0_left,DI1_left,DI0_right,DI1_right,G,m,n);
%     
%     Fdudxatv = scatteredInterpolant(Xlsv, Ylsv, full(ududx_ls));
%     ududx = Fdudxatv(Xu',Yu')'; 
%     
%     Fdudyatv = scatteredInterpolant(Xlsv, Ylsv, full(vdudy_ls));
%     vdudy = Fdudyatv(Xu',Yu')'; 
%     
%     Fdvdxatx = scatteredInterpolant(Xlsv, Ylsv, full(udvdx_ls));
%     udvdx = Fdvdxatx(Xv',Yv')';
%         
%     Fdvdyatx = scatteredInterpolant(Xlsv, Ylsv, full(vdvdy_ls));
%     vdvdy = Fdvdyatx(Xv',Yv')';


    divu = G*[Uls];
    dudx = divu(1:((m+1)*n));
    ududx = Uu.*dudx;

    du_dy_atv = divu(((m+1)*n)+1:end);
    Fdudyatv = scatteredInterpolant(Xvv, Yvv, full(du_dy_atv));
    dudy = Fdudyatv(Xu',Yu')';
    dudy = reshape(dudy.',[],1);
    vdudy = Vu.*dudy;
    
    divv = G*[Vls];
    dv_dx_atu = divv(1:((m+1)*n));
    Fdvdxatx = scatteredInterpolant(Xuv, Yuv, full(dv_dx_atu));
    dvdx = Fdvdxatx(Xv',Yv')';
    dvdx = reshape(dvdx.',[],1);
    udvdx = Uv.*dvdx;
    
    dv_dy_atv = divv(((m+1)*n)+1:end);
    vdvdy = Vv.*dv_dy_atv;
    
     
    %SGS calculation:
for i=1:1   
% % % % %     dudx = reshape(dudx,[m+1 n]).' ;
% % % % %     dvdy = reshape(dvdy,[m n+1]).' ;
% % % % %     
% % % % %     Fdudx = scatteredInterpolant(Xu', Yu', dudx'); 
% % % % % %     Fdudy = scatteredInterpolant(Xu', Yu', dudy'); 
% % % % % %     Fdudy.ExtrapolationMethod = 'linear';
% % % % % %     Fdvdx = scatteredInterpolant(Xv', Yv', dvdx'); 
% % % % % %     Fdvdx.ExtrapolationMethod = 'linear';
% % % % %     Fdvdy = scatteredInterpolant(Xv', Yv', dvdy'); 
% % % % %     
% % % % %     dudxc = Fdudx(Xls',Yls')';
% % % % % %     dudyc = Fdudy(Xls',Yls')';
% % % % % %     dvdxc = Fdvdx(Xls',Yls')';
% % % % %     dvdyc = Fdvdy(Xls',Yls')';
% % % % % %     
% % % % % %     Cs = 0.22;
% % % % % %     Nu = sparse(n,m); 
% % % % % %     
% % % % % %     %u,v dispersion:
% % % % %     Uxls = reshape(dudxc',[], 1);
% % % % % % 
% % % % %     divdudx = G*[Uxls];
% % % % %     d2udx2 = divdudx(1:((m+1)*n));
% % % % %     d2udx2 = reshape(d2udx2,[JMax IMax-1]).' ; 
% % % % %     Fd2udx2atu = scatteredInterpolant(Xu', Yu', d2udx2');
% % % % %     d2u_dx2_atv = Fd2udx2atu(Xv',Yv')';
% % % % %     d2u_dx2_atv = reshape(d2u_dx2_atv',[],1);
% % % % %     d2udx2 = reshape(d2udx2',[],1);
% % % % %     
% % % % %     Vyls = reshape(dvdyc',[], 1);
% % % % %     
% % % % %     divdvdy = G*[Vyls];
% % % % %     d2vdy2 = divdvdy(((m+1)*n+1):end);
% % % % %     d2vdy2 = reshape(d2vdy2,[JMax-1 IMax]).' ;  
% % % % %     Fd2vdy2atv = scatteredInterpolant(Xv', Yv', d2vdy2');
% % % % %     d2v_dy2_atu = Fd2vdy2atv(Xu',Yu')';
% % % % %     d2v_dy2_atu = reshape(d2v_dy2_atu',[],1);
% % % % %     d2vdy2 = reshape(d2vdy2',[],1);
    
    
    
% 
%     for j = 1:m+2
%         for i = 1:n+2
% 
%             tmp = 2.0*( dudxc(i,j)*dudxc(i,j) + dvdyc(i,j)*dvdyc(i,j) +...
%                         dudxc(i,j)*dvdxc(i,j) + dudyc(i,j)*dvdyc(i,j) +...
%                         dudyc(i,j)*dudyc(i,j) + dvdxc(i,j)*dvdxc(i,j));
% 
%             Nu(i,j) = Cs^2*length^2*sqrt(tmp) ;
% 
%         end
%     end
% 
%     tau_xx = -Nu.*(dudxc+dudxc);
%     tau_xy = -Nu.*(dudyc+dvdxc);
%     tau_yy = -Nu.*(dvdyc+dvdyc);
% 
%     tau_xxp = reshape(tau_xx',[], 1);
%     tau_xyp = reshape(tau_xy',[], 1);
%     tau_yyp = reshape(tau_yy',[], 1);
% 
%     divtau_xx = G*[tau_xxp];
%     dtau_xx_dx_atu = divtau_xx(1:((m+1)*n));
%     dtau_xx_dy_atv = divtau_xx(((m+1)*n)+1:end);
% 
%     divtau_xy = G*[tau_xyp];
%     dtau_xy_dx_atu = divtau_xy(1:((m+1)*n));
%     dtau_xy_dy_atv = divtau_xy(((m+1)*n)+1:end);
% 
%     divtau_yy = G*[tau_yyp];
%     dtau_yy_dx_atu = divtau_yy(1:((m+1)*n));
%     dtau_yy_dy_atv = divtau_yy(((m+1)*n)+1:end);
% 
%     dtau_xy_dx_atu = reshape(dtau_xy_dx_atu,[JMax IMax-1]).' ; 
%     dtau_xy_dy_atv = reshape(dtau_xy_dy_atv,[JMax-1 IMax]).' ; 
% 
%     Fdtau_xy_dy_atv = scatteredInterpolant(Xv', Yv', dtau_xy_dy_atv'); 
%     Fdtau_xy_dy_atv.ExtrapolationMethod = 'linear';
%     Fdtau_xy_dx_atu = scatteredInterpolant(Xu', Yu', dtau_xy_dx_atu');
%     Fdtau_xy_dx_atu.ExtrapolationMethod = 'linear';
%     
%     dtau_xy_dx_atv = Fdtau_xy_dx_atu(Xv',Yv')';
%     dtau_xy_dx_atv = reshape(dtau_xy_dx_atv',[], 1);
%     dtau_xy_dy_atu = Fdtau_xy_dy_atv(Xu',Yu')';
%     dtau_xy_dy_atu = reshape(dtau_xy_dy_atu',[], 1);
% 
%     DivSgsx = dtau_xx_dx_atu + dtau_xy_dy_atu;
%     DivSgsy = dtau_xy_dx_atv + dtau_yy_dy_atv;
%
end

    ududx = reshape(ududx.',[], 1);
    vdudy = reshape(vdudy.',[], 1);
    udvdx = reshape(udvdx.',[], 1);
    vdvdy = reshape(vdvdy.',[], 1);
    
        
    %Predict velocities:
    u_star = uatu;
    v_star = vatv;  
    
    %u-momentum prediction:
   
    %RHSu = -(ududx);
    RHSu = -(ududx+vdudy); % + mu/(rho_ref*UStar*LStar)*(d2udx2 + d2v_dy2_atu);    
    %RHSu = -(udu_dx+vdu_dy) + IRe*(d2udx2 + d2v_dy2_atu)    ;
    
    Uu_star = SSPRK101D(Uu_star,RHSu,dt);
    
    u_star = reshape(Uu_star,[JMax IMax-1])' ;

    %v-momentum prediction:
    %RHSv = -(udvdx) - Bb ;
    RHSv =  -(udvdx+vdvdy) - Bb ; % + mu/(rho_ref*UStar*LStar)*(d2u_dx2_atv + d2vdy2) ;   
    %RHSv = -(udv_dx+vdv_dy) - Bb  + IRe*(d2u_dx2_atv + d2vdy2) ;
    
    Vv_star = SSPRK101D(Vv_star,RHSv,dt).';
    
    v_star = reshape(Vv_star,[m n+1])'; 
    
    %u,v prediction boundaries:
    [u_star,v_star,~] = applyboundaries2Dstarcurv(u_star,v_star,p);
    %[u_star,v_star,~] = applyboundaries2Dstar(u_star,v_star,p);
    
    Uu = reshape(u_star.',[], 1);
    Vv = reshape(v_star.',[], 1);
    
    %laplacian solver:
    
    R = [Uu; Vv];
    
    Rorig = D*R*rho_ref/dt;  
    %Rorig = D*R;
    Pp = L\Rorig + p_atm;
    
    %Pp = L\Rorig;
         
    p = reshape(Pp,[m+2 n+2])';    
    % [~,~,p] = applyboundaries2Dcurv(uatu,vatv,p);
     %[~,~,p] = applyboundaries2D(uatu,vatv,p);
     Pp = reshape(p',[], 1);                                                                                  
%     
    %Pp = Pp*LStar;    
    
    %Pp = Pp*rhoStar;    
    
    %gradient of pressure:
    Gpp = Gp*[ Pp ];
    dpdx = Gpp(1:(m+1)*n);
    dpdy = Gpp((m+1)*n+1:end);   
    
    %u,v correction:
    Uu = Uu + dpdx;
    Vv = Vv + dpdy;
    
    uatu = reshape(Uu,[m+1 n]).' ;
    vatv = reshape(Vv,[m n+1]).' ;
    
    dpdx = reshape(dpdx,[m+1 n])' ;
    dpdy = reshape(dpdy,[m n+1])' ;
    
    %velocities boundaries:
    [uatu,vatv,~] = applyboundaries2Dcurv(uatu,vatv,p);
    %[uatu,vatv,~] = applyboundaries2D(uatu,vatv,p);
    
%     fig = figure(jj-1);
%     set(gcf,'units','inches','position',[10 0 width height])
%     u_v = squeeze(dpdx(2:end-1,2:end-1)*UStar);
%     subplot(2,2,1)
%     %pcolor(Xu(2:end-1,2:end-1),Yu(2:end-1,2:end-1),u_v);
%     pcolor(u_v);
%     %axis equal
%     axis tight
%     title('dpdx ')
%     colorbar;
%     shading interp;
%     subplot(2,2,2)
%     u_3 = squeeze(dpdy(2:end-1,2:end-1)*UStar);
%     %pcolor(Xv(2:end-1,2:end-1),Yv(2:end-1,2:end-1),u_3);
%     pcolor(u_3);
%     %axis equal
%     axis tight
%     title('dpdy')
%     colorbar;
%     shading interp;
%     subplot(2,2,3)
%     u_2 = squeeze(p(2:end-1,2:end-1));
%     pcolor(Xls(2:end-1,2:end-1),Yls(2:end-1,2:end-1),u_2);
%     %pcolor(u_2);
%     %axis equal
%     axis tight 
%     title('p [bar]')
%     colorbar;
%     shading interp;
%     subplot(2,2,4)
%     u_2 = squeeze(v_star(2:end-1,2:end-1));
%     pcolor(Xv(2:end-1,2:end-1),Yv(2:end-1,2:end-1),u_2);
%     %pcolor(Xls(2:end-1,2:end-1),Yls(2:end-1,2:end-1),u_2);
%     pcolor(u_2);
%     %axis equal
%     axis tight 
%     title('v star')
%     %title(['T [C] at ' num2str(current_time) 's'])
%     colorbar;
%     shading interp;
    
    Uu = reshape(uatu.',[], 1);
    Vv = reshape(vatv.',[], 1);
    
    %interpolating u,v to center for temperature eq:
    Fuatu = scatteredInterpolant(Xuv, Yuv, full(Uu));
    Fvatv = scatteredInterpolant(Xvv, Yvv, full(Vv));
    u = Fuatu(Xls',Yls')';
    v = Fvatv(Xls',Yls')';
    
    Uls = reshape(u.',[], 1);
    Vls = reshape(v.',[], 1);    
    
    uatc = Fuatu(Xls',Yls')';
    vatc = Fvatv(Xls',Yls')'; 
              
    uatv = Fuatu(Xv',Yv')';
    vatu = Fvatv(Xu',Yu')';

    Uv = reshape(uatv.',[], 1);
    Vu = reshape(vatu.',[], 1);  
    
    %Convection-Difussion equation:
    uvec = [Uu;Vv];

    Tvec = reshape(T.',[], 1);

%     SoluT = MimeticUpwind(Tvec,uvec,I0,I1);
%     uDT   = D*(uvec.*SoluT');  
%     RHST  = -uDT ; %+ D2T ;
     
    [udTdx_ls,vdTdy_ls] = MimeticUpwindv3grad(Tvec,Uls,Vls,DI0_left,DI1_left,DI0_right,DI1_right,G,m,n);
    %RHST =  -(udTdx_ls ); %+ vdTdy_ls);
    RHST =  -(udTdx_ls + vdTdy_ls);
    
    T_new = SSPRK101D(Tvec,RHST,dt);
    T     = reshape(T_new,[JMax+1 IMax+1])';   
       
    %Apply boundaries to T:
    T(2,:) = T(3,:);
    T(end-1,:) = T(end-2,:);
    T(:,2) = T(:,3);
    T(:,end-1) = T(:,end-2);
    %
    T(1,:) = T(2,:);
    T(end,:) = T(end-1,:);
    T(:,1) = T(:,2);
    T(:,end) = T(:,end-1);
   
    dens = rho_ref*( 1 - alpha*(T*T0-T_ref));
    buoy = gforce*(dens-rho_P0)/rho_P0;
    
    buoy = reshape(buoy.',[],1);   
    
    Fbuoy = scatteredInterpolant(Xlsv, Ylsv, buoy);
    batv = Fbuoy(Xv',Yv')';   
    
    current_time = (time+1)*dt*LStar/UStar
        
    p = reshape(Pp,[JMax+1 IMax+1])';
    
 if(mod(current_time,print_freq) <= 1e-10 )
    
    fig = figure(jj);
    set(gcf,'units','inches','position',[10 10 width height])
    u_v = squeeze(uatu(2:end-1,2:end-1)*UStar);
    subplot(2,2,1)
    pcolor(Xu(2:end-1,2:end-1),Yu(2:end-1,2:end-1),u_v);
    %pcolor(u_v);
    %axis equal
    axis tight
    title('u [m/s]')
    colorbar;
    shading interp;
    subplot(2,2,2)
    u_3 = squeeze(vatv(2:end-1,2:end-1)*UStar);
    pcolor(Xv(2:end-1,2:end-1),Yv(2:end-1,2:end-1),u_3);
    %pcolor(u_3);
    %axis equal
    axis tight
    title('v [m/s]')
    colorbar;
    shading interp;
    subplot(2,2,3)
    u_2 = squeeze(p(2:end-1,2:end-1));
    pcolor(Xls(2:end-1,2:end-1),Yls(2:end-1,2:end-1),u_2);
    %pcolor(u_2);
    %axis equal
    axis tight 
    title('p [bar]')
    colorbar;
    shading interp;
    subplot(2,2,4)
    u_2 = squeeze(T(2:end-1,2:end-1)*T0);
    pcolor(Xls(2:end-1,2:end-1),Yls(2:end-1,2:end-1),u_2);
    %pcolor(u_2);
    %axis equal
    axis tight 
    title(['T [C] at ' num2str(current_time) 's'])
    colorbar;
    shading interp;
    
    temperature_save = [temperature_save; {T}];
    velocities_save = [velocities_save; {uatu}];
    momentum_save = [momentum_save, {udvdx+vdvdy}];
    momentum_savex = [momentum_savex, {ududx+vdudy}];
    pgrad_save = [pgrad_save, {dpdy}];
    pgrad_savex = [pgrad_savex, {dpdx}];
    buoy_save = [buoy_save, {Bb}];
    
    f_name = ['LE_fullmimetic_' num2str(JMax) '_dt.01_' num2str(current_time) 's.png'];
    %saveas(fig,f_name,'png')


    current_time
    
 end
    
end


%%
n_levels = 15;

frame_times = [5 10 30 50];

fr_1 = frame_times(1);
fr_2 = frame_times(2);
fr_3 = frame_times(3);
fr_4 = frame_times(4);

fig = figure(3)
set(gcf,'units','inches','position',[10 10 20 20])
subplot(4,1,1)
u_frame = temperature_save{fr_1,1};
u_2 = squeeze(u_frame(2:end-1,2:end-1)*T0);
[cont h]= contour(Xls(2:end-1,2:end-1),Yls(2:end-1,2:end-1),u_2,n_levels);  h.LineColor= 'black'; h.LineWidth = 1.0;
xlabel('x [m]')
ylabel('y [m]')
title(['T [C] contours at ' num2str(fr_1) 's'])
axis equal
axis tight 
subplot(4,1,2)
u_frame = temperature_save{fr_2,1};
u_2 = squeeze(u_frame(2:end-1,2:end-1)*T0);
[cont h]= contour(Xls(2:end-1,2:end-1),Yls(2:end-1,2:end-1),u_2,n_levels);  h.LineColor= 'black'; h.LineWidth = 1.0;
title(['T [C] contours at ' num2str(fr_2) 's'])
xlabel('x [m]')
ylabel('y [m]')
axis equal
axis tight 
subplot(4,1,3)
u_frame = temperature_save{fr_3,1};
u_2 = squeeze(u_frame(2:end-1,2:end-1)*T0);
[cont h]= contour(Xls(2:end-1,2:end-1),Yls(2:end-1,2:end-1),u_2,n_levels);  h.LineColor= 'black'; h.LineWidth = 1.0;
title(['T [C] contours at ' num2str(fr_3) 's'])
xlabel('x [m]')
ylabel('y [m]')
axis equal
axis tight 
subplot(4,1,4)
u_frame = temperature_save{fr_4,1};
u_2 = squeeze(u_frame(2:end-1,2:end-1)*T0);
[cont h]= contour(Xls(2:end-1,2:end-1),Yls(2:end-1,2:end-1),u_2,n_levels);  h.LineColor= 'black'; h.LineWidth = 1.0;
title(['T [C] contours at ' num2str(fr_4) 's'])
xlabel('x [m]')
ylabel('y [m]')
axis equal
axis tight 
f_name = ['LE_fullmimetic_contour_' num2str(JMax) '_dt.1_' num2str(current_time) 's.png'];
%saveas(fig,f_name,'png')

%%
u_f = [];
u_b = sqrt(abs(r_g*(d-c)/2));
%u_b = 0.0224;

for i = 1:numel(velocities_save)-1
    u_frame = velocities_save{i,1};
    u_f(i) = max((abs(u_frame(2,:))))*UStar;
end

Fr1 = u_f./u_b;

mFr = mean(Fr1);

for i = 1:numel(velocities_save)-1
    u_frame = velocities_save{i,1};
    u_f(i) = max((abs(u_frame(end-1,:))))*UStar;
end

Fr2 = u_f./u_b;

mFr2 = mean(Fr2);

mFravg = (mFr + mFr2)/2.0;
diff = 100-mFravg/0.7012*100;
width = 12;
height = 5;

fig = figure(77);
set(gcf,'units','inches','position',[5 2 width height])
plot(Fr1,'+-','LineWidth',2);
hold on
plot(Fr2,'LineWidth',2);
yline(0.7012,'--')
hold off
legend(['Bottom=' num2str(mFr)],['Top=' num2str(mFr2)],'Theoretical','Location','southeast')
title(['LE ' num2str(JMax) 'x' num2str(IMax) ' cells, Froude number mean = ' num2str(mFravg) ' err(%)=' num2str(diff)])
xlabel('time [s]')
shading flat;

f_name = ['LE_' num2str(JMax) '_dt01_Froude'];
%saveas(fig,f_name,'png')

%%

mFr = mean(Fr1(15:end-5));
mFr2 = mean(Fr2(15:end-5));

mFravg = (mFr + mFr2)/2.0;

diff = 100-mFravg/0.7012*100;

fig = figure(77);
set(gcf,'units','inches','position',[5 2 width height])
plot(Fr1(15:end-5),'+-','LineWidth',2);
hold on
plot(Fr2(15:end-5),'LineWidth',2);
yline(0.7012,'--')
hold off
legend(['Bottom=' num2str(mFr)],['Top=' num2str(mFr2)],'Theoretical','Location','southeast')
title(['LE ' num2str(JMax) 'x' num2str(IMax) ' cells, Froude number mean = ' num2str(mFravg) ' err(%)=' num2str(diff)])
xlabel('time [s]')
shading flat;

f_name = ['LE_stable_' num2str(JMax) '_dt01_Froude'];
%saveas(fig,f_name,'png')


%%
mom_f = [];
mom_fx = [];
pg_f = [];
pg_fx = [];

b_f = [];

for i = 1:numel(momentum_save)
    mom_f(i)  = mean(abs(momentum_save{:,i}));
    mom_fx(i) = mean(abs(momentum_savex{:,i}));
    pg_f(i)   = mean(abs(pgrad_save{:,i}));
    pg_fx(i)  = mean(abs(pgrad_savex{:,i}));
    b_f(i)    = mean(abs(buoy_save{:,i}));
end

f_name = figure(78);
set(gcf,'units','inches','position',[10 0 width height])
figure(78);
semilogy(mom_f,'+-');
hold on
semilogy(mom_fx,'--');
hold on
semilogy(pg_f,'^-');
hold on
semilogy(pg_fx,'v-');
hold on
semilogy(b_f,'s-');
hold off
legend(['v-advection=' num2str(mom_f(end))],['u-advection=' num2str(mom_fx(end))],...
['dpdy=' num2str(pg_f(end))],['dpdx=' num2str(pg_fx(end))],['buoyancy=' num2str(b_f(end))],'Location','southeast')
title(['Magnitude balance for LE'])
xlabel('time (s)')
shading flat;
f_name = ['LE_' num2str(JMax) '_magnitudes'];
saveas(fig,f_name,'png')

