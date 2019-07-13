function [envelope,pha,fy1,fy2,peak,numt0r,numt0,deltat,delta]=PSTD(fcenter,cellsperwavelength,BW,sigType,time_tot,time_shift,x_dist,y_dist,source_height,surface_height,TFSF_type,PML_type,...
    obs_x,obs_y,obs_x_r,obs_y_r,plane_wave,focus,ref_corr)

close all;

%%                    INPUT                   

%%%%%%%%% fcenter (in Hz): center frequency being simulated
%%%%%%%%% cellsperwavelength: number of grid cells per wavelength at center frequency
%%%%%%%%% BW (in Hz): transmitted bandwidth (for Gaussian pulse only)
%%%%%%%%% sigType: flag of the transmitted signal type, i.e. 0 for Gaussian pulse, 1 for Blackman-Harris (BH) pulse, and 2 for sinc pulse
%%%%%%%%% time_tot (in s): the total time period being simulated, roughly given by 2*time_shift + (source2interface)/c + (interface2surface)/c
%%%%%%%%% time_shift (in s): single-sided pulse width under 3-sigma rule, roughly given by round(0.966/BW/2*3/deltat) for Gaussian pulse
%%%%%%%%% x_dist (in m): horizontal dimension of the domain
%%%%%%%%% y_dist (in m): vertical dimension of the domain
%%%%%%%%% source_height (in m): distance between TF/SF boundary and domain boundary
%%%%%%%%% surface_height (in m): distance between Huygens' contour and domain boundary
%%%%%%%%% TFSF_type: flag of the TF/SF domain type, i.e. 0 for 4 walls (point target), 1 for 1 back wall (distributed target under normal plane-wave incidence), and 2 for 1 back wall (distributed target under arbitrary wave incidence)
%%%%%%%%% PML_type: flag of the PML type, i.e. 0 for Berenger's PML (lossless background medium), 1 for Generalized PML (GPML; lossy background medium)
%%%%%%%%% obs_x (in m): horizontal coordinate of the radar transmitter (must be scalar)
%%%%%%%%% obs_y (in m): vertical coordinate of the radar transmitter (must be scalar)
%%%%%%%%% obs_x_r (in m): horizontal coordinate of the radar receiver (can be scalar or vector; always same dimension as obs_y_r)
%%%%%%%%% obs_y_r (in m): vertical coordinate of the radar receiver (can be scalar or vector; always same dimension as obs_x_r)
%%%%%%%%% plane_wave: flag identifying plane-wave source, i.e. 0 not plane wave, 1 plane wave
%%%%%%%%% focus: flag performing focusing, i.e. 0 no focusing, non-zero value focusing with the value representing azimuth resolution in meters, '-' means default azimuth resolution equals to range resolution. When "obs_x_r" and "obs_y_r" are vectors, any non-zero value (or '-') of "focus" will apply focusing on the user-specified receiver locations.
%%%%%%%%% ref_corr: flag of refraction correction in the focused radargram, i.e. 0 no correction, non-zero-value vector "[a,b]" means correction with "a" representing the surface elevation in meters, "b" representing the dielectric constant of volume under surface

%%                   OUTPUT                   

%%%%%%%%% envelope: time-domain amplitude of the far field (can be vector or matrix; number of rows always same as length of obs_y_r)
%%%%%%%%% pha: time-domain phase of the far field (can be vector or matrix; number of rows always same as length of obs_y_r)
%%%%%%%%% fy1: frequency-domain amplitude of the far field (can be vector or matrix; number of rows always same as length of obs_y_r)
%%%%%%%%% fy2: frequency-domain amplitude of the incident field (can be vector or matrix; number of rows always same as length of obs_y_r)
%%%%%%%%% peak: ratio of fy1 and fy2 at center frequency for Radar Cross Width (RCW; defined as "2*pi*sqrt(obs_x_r.^2+obs_y_r.^2).*peak.^2") calcualtions (can be scalar or column vector; size always same as length of obs_y_r)
%%%%%%%%% numt0r: number of temporal grid cells propagating from Near2Far Huygens' contour to each receiver point (can be scalar or column vector; size always same as length of obs_y_r)
%%%%%%%%% numt0: number of temporal grid cells propagating from transmitter point to TF/SF boundary (scalar; only for active focusing, column vector with size same as length of obs_y_r)
%%%%%%%%% deltat: temporal grid sampling resolution (in s; scalar)
%%%%%%%%% delta: spatial grid sampling resolution (in m; scalar)

%%%%%%%%% NOTE: "numt0r" and "numt0" are the timing variables used to reconstruct the time series of far field as well as the focused radargram.
%%%%%%%%%       For example, when accessing values in "envelope" and/or "pha", the time stamp for a given temporal grid index is calculated as:
%%%%%%%%%       "(temporal_grid_index + numt0 + numt0r)*deltat";
%%%%%%%%%       in reverse, the temporal grid index for a given target (where the envelope peak is) is expressed as below:
%%%%%%%%%       "(modeled_round_trip_distance_between_target_and_radar/speed_of_light/deltat + time_shift/deltat - numt0 - numt0r)".
%%%%%%%%%       The calculated temporal grid index is then used to access the corresponding far-field value from "envelope" and/or "pha", where the complex far-field value is represented as "envelope.*exp(1j*pha)".

%%                    sample run                   
%%%%%%%%%   passive focusing (automatically generated receiver locations by specifiying the receiver array elevation and desired azimuth resolution)
%%%%%%%%%   [envelope,pha,fy1,fy2,peak,numt0r,numt0,deltat,delta]=PSTD(20e6,16,10e6,0,1.5e-6,0.3e-6,400,160,40,25,1,0,0,6000,0,6000,1,'-',[-40,2.3])
%%%%%%%%%   
%%%%%%%%%   passive focusing (user-specified receiver locations)
%%%%%%%%%   [envelope,pha,fy1,fy2,peak,numt0r,numt0,deltat,delta]=PSTD(20e6,16,10e6,0,1.5e-6,0.3e-6,400,160,40,25,1,0,0,6000,[-3000,0,3000],[6000,6000,6000],1,'-',[-40,2.3])
%%%%%%%%%   
%%%%%%%%%   active focusing (automatically generated receiver locations by specifiying the receiver array elevation and desired azimuth resolution)
%%%%%%%%%   [envelope,pha,fy1,fy2,peak,numt0r,numt0,deltat,delta]=PSTD(20e6,16,10e6,0,1.5e-6,0.3e-6,400,160,60,45,2,1,0,6000,0,6000,0,'-',[-40,2.3])
%%%%%%%%%   
%%%%%%%%%   active focusing (user-specified receiver locations)
%%%%%%%%%   [envelope,pha,fy1,fy2,peak,numt0r,numt0,deltat,delta]=PSTD(20e6,16,10e6,0,1.5e-6,0.3e-6,400,160,60,45,2,1,0,6000,[-3000,0,3000],[6000,6000,6000],0,'-',[-40,2.3])

%%
tic

%%  Define phsyical attributes of PSTD and waveforms

simspace_waveforms;

%%  Plot waveforms and spectral characteristics

waveform_characteristics;

%%  Define simulation canvas

canvas;

%%  Near2Far Huygens' contour based on type of TF/SF formulation

source_height = round(source_height/delta);
surface_height = round(surface_height/delta);
center_location = round(x_dist/2/delta);
center_height = source_height;

if TFSF_type==0
    surface_location=surface_height:xdim-surface_height;
    surface_location_L=surface_height;
    surface_location_R=xdim-surface_height;
    surface_location_B=surface_height:xdim-surface_height;
    surface_height_L=(ydim-surface_height):-1:surface_height;
    surface_height_R=(ydim-surface_height):-1:surface_height;
    surface_height_B=ydim-surface_height;
    Ez_data = zeros(time_tot,length(surface_location));
    Ez_data_L = zeros(time_tot,length(surface_height_L));
    Ez_data_R = zeros(time_tot,length(surface_height_R));
    Ez_data_B = zeros(time_tot,length(surface_location_B));
    Hx_data = zeros(time_tot,length(surface_location));
    Hx_data_L = zeros(time_tot,length(surface_height_L));
    Hx_data_R = zeros(time_tot,length(surface_height_R));
    Hx_data_B = zeros(time_tot,length(surface_location_B));
elseif TFSF_type==1
    surface_location=1:xdim;
    Ez_data = zeros(time_tot,length(surface_location));
    Hx_data = zeros(time_tot,length(surface_location));
elseif TFSF_type==2
    surface_location=1:xdim;
    Ez_data = zeros(time_tot,length(surface_location));
    Hx_data = zeros(time_tot,length(surface_location));
else
    error('TF/SF formulation NOT defined!!!')
end

F(time_tot) = struct('cdata',[],'colormap',[]);

%%  Define dielectric scene

% Here is the only place where human interference is needed: the dielectric scene parameters should be inputted (prior to the PSTD run) in order to define the scene
dielectric_scene;

%%  View canvas

view_canvas;

%%  Define PML layers

if PML_type==0
    PML_define;
elseif PML_type==1
    GPML_define;
else
    error('PML type NOT defined!!!')
end

%%  Update Maxwell's coefficients

update_Maxwell;

%%  Define spectral propagators for PSTD

k_space;

%%  Solve Maxwell's equation via PSTD update

if (TFSF_type==0)&&(plane_wave==1)&&(PML_type==0)&&(focus==0)&&(obs_x==0)       %%%%%%%%% Normal-incidence plane-wave point target (input "obs_x=0" for normal incidence) 
    if sigType==0
        TFSF_PW_box_CR;
    elseif sigType==1
        TFSF_PW_BH_box_CR;
    else
        error('Waveform NOT defined!!!')
    end
    if length(obs_y_r)==1
        obs_y_r_sub=obs_y_r;obs_x_r_sub=obs_x_r;
        N2F_TFSF;
    else
        vectorized_receivers;
    end
elseif (TFSF_type==1)&&(plane_wave==1)&&(PML_type==0)&&(focus==0)&&(obs_x==0)   %%%%%%%%% Normal-incidence plane-wave distributed target (input "obs_x=0" for normal incidence)
    if sigType==0
        TFSF_PW_surface_CR;
    elseif sigType==1
        TFSF_PW_BH_surface_CR;
    else
        error('Waveform NOT defined!!!')
    end
    if length(obs_y_r)==1
        obs_y_r_sub=obs_y_r;obs_x_r_sub=obs_x_r;
        N2F_TFSF_surface;
    else
        vectorized_receivers_surface;
    end
elseif (TFSF_type==2)&&(PML_type==1)&&(focus==0)&&(length(obs_x)==1)            %%%%%%%%% Oblique-incidence arbitrary wave distributed target
    if sigType==0
        TFSF_surface_CR_GPML;
    else
        error('Waveform NOT defined!!!')
    end
    if length(obs_y_r)==1
        obs_y_r_sub=obs_y_r;obs_x_r_sub=obs_x_r;
        N2F_TFSF_surface;
    else
        vectorized_receivers_surface;
    end
elseif (TFSF_type==1)&&(plane_wave==1)&&(PML_type==0)&&(focus~=0)&&(obs_x==0)   %%%%%%%%% Passive (normal-incidence plane-wave) focusing (input "obs_x=0" for normal incidence)
    if sigType==0
        TFSF_PW_surface_CR;
    else
        error('Waveform NOT defined!!!')
    end
    radargram_result;
    passive_focus_result;
    if ref_corr~=0
        passive_focus_corr_result;
    end
elseif (TFSF_type==2)&&(PML_type==1)&&(focus~=0)&&(length(obs_x)==1)            %%%%%%%%% Active focusing (inputs "obs_x" and "obs_y" are meaningless, because they are same as "obs_x_r" and "obs_y_r" for active focusing)
    active_focus_result;
    if ref_corr~=0
        active_focus_corr_result;
    end
else
    envelope=0;pha=0;fy1=0;fy2=0;peak=0;numt0r=0;numt0=0;deltat=0;delta=0;
    error('TF/SF formulation NOT defined!!!')
end

%%

toc