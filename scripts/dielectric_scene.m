%% Flag for target type

target_type_surfnvol=1;     %%%%%%%%% Dielectric rough surface and heterogeneous volume
target_type_point=0;        %%%%%%%%% Single-grid-cell point target
target_type_cylinder=1;     %%%%%%%%% Dielectric Cylinder
target_type_custom=0;       %%%%%%%%% Customized dielectric scene, i.e. user-defined

epsilon1 = epsilon;
conductivity1 = conductivity;

%% Dielectric rough surface and heterogeneous volume

if target_type_surfnvol == 1
    %%%%%%%%%                                                      %%%%%%%%%
    %%%%%%%%%                   INPUT PARAMETERS                   %%%%%%%%%
    %%%%%%%%%                                                      %%%%%%%%%
    surface_position=-40;       %%%%%%%%% mean elevation of rough surface in meters
    rmsh_1D=0.25;               %%%%%%%%% rms height of rough surface in terms of wavelength
    lc_1D=2;                    %%%%%%%%% correlation length of rough surface in terms of wavelength
    type_1D='norm';             %%%%%%%%% type of correlation function for the rough surface: 'norm' for Gaussian, 'exp' for exponential
    std_eps=5/100;              %%%%%%%%% percentile rms permittivity of heterogeneous volume
    lc_2D=5;                    %%%%%%%%% correlation length of heterogeneous volume in terms of meters
    type_2D='exp';              %%%%%%%%% type of correlation function for the heterogeneous volume: 'norm' for Gaussian, 'exp' for exponential
    relative_permittivity=2.3;  %%%%%%%%% mean permittivity of heterogeneous volume (could be complex, e.g. formatted as "a-1jb")
    sigma=0e2;                  %%%%%%%%% mean conductivity of heterogeneous volume
    %%%%%%%%%                                                      %%%%%%%%%
    %%%%%%%%%                                                      %%%%%%%%%
    %%%%%%%%%                                                      %%%%%%%%%
    
    surface_position=-round(surface_position/delta)+center_height;
    
    [epsilon1,conductivity1] = dielectric_scene_surfnvol(surface_position,rmsh_1D,lc_1D,type_1D,std_eps,lc_2D,type_2D,relative_permittivity,sigma,...
        epsilon1,conductivity1,delta,epsilon0,omega,wavelength);
end

%% Single-grid-cell point target

if target_type_point == 1
    %%%%%%%%%                                                      %%%%%%%%%
    %%%%%%%%%                   INPUT PARAMETERS                   %%%%%%%%%
    %%%%%%%%%                                                      %%%%%%%%%
    target_height=-80;          %%%%%%%%% height of point target in meters
    target_location=100;        %%%%%%%%% horizontal location of point target in meters
    relative_permittivity=3;    %%%%%%%%% permittivity of point target (could be complex, e.g. formatted as "a-1jb")
    sigma=2e2;                  %%%%%%%%% conductivity of point target
    %%%%%%%%%                                                      %%%%%%%%%
    %%%%%%%%%                                                      %%%%%%%%%
    %%%%%%%%%                                                      %%%%%%%%%
    
    target_height=-round(target_height/delta)+center_height;target_location=round(target_location/delta)+center_location;
    
    epsilon1(target_height:target_height+1,target_location:target_location+1) = relative_permittivity;
    conductivity1(target_height:target_height+1,target_location:target_location+1) = sigma;
end

%% Dielectric Cylinder

if target_type_cylinder == 1
    %%%%%%%%%                                                      %%%%%%%%%
    %%%%%%%%%                   INPUT PARAMETERS                   %%%%%%%%%
    %%%%%%%%%                                                      %%%%%%%%%
    target_height=-80;          %%%%%%%%% height of cylinder in meters
    target_location=0;          %%%%%%%%% horizontal location of cylinder in meters
    radius=11.5;                %%%%%%%%% radius of cylinder in meters
    relative_permittivity=3;    %%%%%%%%% permittivity of cylinder (could be complex, e.g. formatted as "a-1jb")
    sigma=2e2;                  %%%%%%%%% conductivity of cylinder
    %%%%%%%%%                                                      %%%%%%%%%
    %%%%%%%%%                                                      %%%%%%%%%
    %%%%%%%%%                                                      %%%%%%%%%
    
    target_height=-round(target_height/delta)+center_height;target_location=round(target_location/delta)+center_location;
    
    [epsilon1,conductivity1] = dielectric_scene_cylinder(target_height,target_location,radius,relative_permittivity,sigma,...
        epsilon1,conductivity1,delta,epsilon0,omega);
end

%% Customized dielectric scene, i.e. user-defined

if target_type_custom == 1
    load('dielectric_scene.mat','epsilon1','conductivity1');
end

%% 

figure;
ax(1)=subplot(2,1,1);imagesc(delta*1*(1:1:xdim),(delta*1*(1:1:ydim))',epsilon1./epsilon0,[1 3]);
xlabel('Horizontal direction (m)'); ylabel('Depth direction (m)'); 
cbarlabel = colorbar; ylabel(cbarlabel,'Relative Permitivity')
axis equal
axis([0 x_dist 0 y_dist])

% figure;
ax(2)=subplot(2,1,2);imagesc(delta*1*(1:1:xdim),(delta*1*(1:1:ydim))',conductivity1,[0 1e-6]);
xlabel('Horizontal direction (m)'); ylabel('Depth direction (m)'); 
cbarlabel = colorbar; ylabel(cbarlabel,'Conductivity')
axis equal
axis([0 x_dist 0 y_dist])

linkaxes(ax);