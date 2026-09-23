% Boussinesq Pred-Corr formulation:
for time = 1:size(time_steps,2) 
    RHST = zeros(m+2,n+2,o+2);
    %Predict velocities:
    Bb   = batv(:); 
   
%     %'u'
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
     
    %RHSu = -(wdu_dz);
    RHSu = -(udu_dx(:)+vdu_dy(:)+wdu_dz(:)); %+ DivSgsx; 
    RHSv = -(udv_dx(:)+vdv_dy(:)+wdv_dz(:)); %+ DivSgsy;.
    RHSw = -(udw_dx(:)+vdw_dy(:)+wdw_dz(:)); %+ DivSgsx;
%     
    u_star = SSPRK101D(Uu,RHSu,dt);
    v_star = SSPRK101D(Vv,RHSv,dt);        
    %w_star = SSPRK101D(Ww,RHSw,dt);
    w_star = SSPRK101D(Ww,RHSw-Bb,dt);
    
%     u_star = Uu + dt*RHSu;
%     v_star = Vv + dt*RHSv;
%     w_star = Ww + dt*(RHSw-Bb); 
    
    u_star = reshape(full(u_star),[m+1 n o]); 
    v_star = reshape(full(v_star),[m n+1 o]); 
    w_star = reshape(full(w_star),[m n o+1]); 
       
   % [u_star,v_star,w_star,~] = applyboundaries3Dstarfreeslip(u_star,v_star,w_star,p);
   % [u_star,v_star,w_star,~] = applyboundaries3Dstar(u_star,v_star,w_star,p);
   [u_star,v_star,w_star,~] = applyboundaries3D(u_star,v_star,w_star,p);
   % [u_star,v_star,w_star,~] = applyboundaries3Dfreeslip(u_star,v_star,w_star,p);

        
    %Calculate pressure:
    Uu = u_star(:);
    Vv = v_star(:);
    Ww = w_star(:);
       
    R = [Uu; Vv; Ww];
    Pp = L\(D*R*rho_ref/dt);
    %Pp = L\(D*R/dt);  
    
    p = reshape(Pp,[m+2 n+2 o+2]);  
    %[~,~,~,p] = applyboundaries3Dstarfreeslip(u_star,v_star,w_star,p);
    %[~,~,~,p] = applyboundaries3Dstar(u_star,v_star,w_star,p);
    %Pp = p(:);
       
    %Derivate pressure:
    Gpp = Gp*[ Pp ];
    dpdx = Gpp(1:u_count);
    dpdy = Gpp(u_count+1:v_count);
    dpdz = Gpp(v_count+1:end);
    
    %Correct velocities:
    Uu = Uu + dpdx;
    Vv = Vv + dpdy;
    Ww = Ww + dpdz;
      
    uatu = reshape(Uu,[m+1 n o]);
    vatv = reshape(Vv,[m n+1 o]);
    watw = reshape(Ww,[m n o+1]);
       
    current_time = (time+1)*dt*LStar/UStar
    
   %beam_external;
   [uatu,vatv,watw,~] = applyboundaries3D(uatu,vatv,watw,p);

    Uu = uatu(:);
    Vv = vatv(:);
    Ww = watw(:);
        
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
 
    Tvec = T(:);
    uvec = [Uu;Vv;Ww];
        
    Solplus  = I0*Tvec;
    Solminus = I1*Tvec;
     
    aplus  = max(uvec,0)';
    aminus = min(uvec,0)';
     
     SoluT = aplus.*Solminus' + aminus.*Solplus'; 
     uDT   = D*(SoluT');
     
     RHST  = -uDT ; 
    
     T_new = SSPRK101D(Tvec,RHST,dt);
 
    %T_new = Tvec + dt*RHST;
    
     T     = reshape(T_new,[m+2 n+2 o+2]);   
    
     T(2,:,:) = T(3,:,:);
     T(:,2,:) = T(:,3,:);
     T(:,:,2) = T(:,:,3);
     
     T(1,:,:) = T(2,:,:);
     T(:,1,:) = T(:,2,:);
     T(:,:,1) = T(:,:,2);
    
     T(end-1,:,:) = T(end-2,:,:);
     T(:,end-1,:) = T(:,end-2,:);   
     T(:,:,end-1) = T(:,:,end-2);
   
     T(end,:,:) = T(end-1,:,:);
     T(:,end,:) = T(:,end-1,:);
     T(:,:,end) = T(:,:,end-1);

     %Calculate density and buoyancy:
     %dens = rho_ref + drho_dT*(T*T0-T_ref);  %EOS  
     %dens = dens/rhoStar;
    
     dens = rho_ref + drho_dT*(T-T_ref);  %EOS  
     %dens = rho_ref*( 1 - alpha*(T-T_ref));
     buoy = gforce*(dens-rho_ref)/rho_ref; 
    
     %Fbuoy = scatteredInterpolant(Xls(:), Yls(:), Zls(:), buoy(:));
     %batv = Fbuoy(Xw,Yw,Zw);

     for i=2:m
     for j=2:n
     for k=2:o+1
         batv(i,j,k) = (buoy(i,j,k) + buoy(i,j,k-1))/2;
     end
     end
     end
     
    
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
     
    %velocities_save = [velocities_save; {uatu(:,:,:)}];
    %temperature_save = [temperature_save; {T(:,:,:)}];

%Visualization: 
 if(mod(current_time,print_freq) <= 1e-10 )
    
    width = 16;
    height = 7;
    slice = 3;
    
    fig = figure(jj);
    set(gcf,'units','inches','position',[10 10 width height])
    u_v = squeeze(uatu(:,slice,:)*UStar);
    subplot(2,2,1)
    %contourf(u_v',15);
    pcolor(u_v');
    axis equal
    axis tight
    title(['u max=' num2str(max(max(max(uatu))))])
    colorbar;
    shading interp;
    subplot(2,2,2)
    u_3 = squeeze(watw(:,slice,:)*UStar);
    pcolor(u_3');
    %contourf(u_3',15);
    axis equal
    axis tight
    title(['w max=' num2str(max(max(max(watw))))])
    colorbar;
    shading interp;
    subplot(2,2,3)
    u_2 = squeeze(p(:,slice,:));
    pcolor(u_2');
    axis equal
    axis tight
    title(['p max=' num2str(max(max(max(p))))])
    colorbar;
    shading interp;
    subplot(2,2,4)
    u_3 = squeeze(T(:,slice,:));
    pcolor(u_3');
    axis equal
    axis tight
    title(['T  max=' num2str(max(max(max(T)))) ' at ' num2str(current_time) ])
    colorbar;
    shading interp;
    velocities_save  = [velocities_save ; {uatu(:,:,:)}];
    velocities_savev = [velocities_savev; {vatv(:,:,:)}];
    velocities_savew = [velocities_savew; {watw(:,:,:)}];
    density_save     = [density_save    ; {dens(:,:,:)}];
    temperature_save = [temperature_save; {T(:,:,:)}];
    timedim = [timedim; {current_time}];

    
    %jj = jj+1;  
    f_name = ['LE3D_' num2str(IMax) '_dt' num2str(dt) '_' num2str(current_time) 's.png'];
    %saveas(fig,f_name,'png')
 end
     
end