%%%%%%%%%%%%%%%%%%%%%%%%%%  Connecting Region 


ibh_th=10;
ibh_x=(0:1/(ibh_th-1):1);
ibh_y=1*sin(2*pi*eps*ibh_x)/(2*pi*eps)-1.3611*sin(2*pi*1*ibh_x)/(2*pi*1)+0.3938*sin(2*pi*2*ibh_x)/(2*pi*2)-0.0326*sin(2*pi*3*ibh_x)/(2*pi*3);

% zeta_x=[zeros(1,source_height-ibh_th) ibh_y ones(1,xdim-2*source_height) fliplr(ibh_y) zeros(1,source_height-ibh_th)];
% zeta_x=ones(1,xdim);
% zeta_x=[ibh_y ones(1,xdim-2*ibh_th) fliplr(ibh_y)];
% zeta_x=[ibh_y ones(1,bound_width_x-ibh_th) ones(1,xdim-2*bound_width_x) ones(1,bound_width_x-ibh_th) fliplr(ibh_y)];
zeta_x=[zeros(1,10) ibh_y ones(1,bound_width_x-ibh_th-10) ones(1,xdim-2*bound_width_x) ones(1,bound_width_x-ibh_th-10) fliplr(ibh_y) zeros(1,10)];
% zeta_y=[zeros(1,source_height-ibh_th) ibh_y ones(1,ydim-2*source_height) fliplr(ibh_y) zeros(1,source_height-ibh_th)];
zeta_y=[zeros(1,source_height-ibh_th) ibh_y ones(1,ydim-source_height-ibh_th) fliplr(ibh_y)];

zeta=zeta_y.'*zeta_x;



%%%%%%%%%%%%%%%%%%%%%%%% Integral components


EzyI=zeros(size(zeta));
EzxI=zeros(size(zeta));
HyI=zeros(size(zeta));
HxI=zeros(size(zeta));





%%%%%%%%%%%%%%%%%%%%%%%%%% TFSF box definition


TFSF_height_T = source_height:-1:(source_height-ibh_th+1);
% TFSF_location_T=(source_height-ibh_th+1):(xdim-source_height+ibh_th-1);
TFSF_location_T=(1):(xdim);
[TX,TY]=meshgrid(TFSF_location_T,TFSF_height_T);

% TFSF_height_B = (ydim-source_height+ibh_th-1):-1:(ydim-source_height);
% TFSF_location_B=(source_height-ibh_th+1):(xdim-source_height+ibh_th-1);
% [BX,BY]=meshgrid(TFSF_location_B,TFSF_height_B);
% 
% TFSF_height_L = (ydim-source_height-1):-1:(source_height+1);
% TFSF_location_L=(source_height-ibh_th+1):source_height;
% [LX,LY]=meshgrid(TFSF_location_L,TFSF_height_L);
% 
% TFSF_height_R = (ydim-source_height-1):-1:(source_height+1);
% TFSF_location_R=(xdim-source_height):(xdim-source_height+ibh_th-1);
% [RX,RY]=meshgrid(TFSF_location_R,TFSF_height_R);




%%%%%%%%%%%%%%%%%%%%% TFSF triangle positions

TFSF_x_tri_T = (TX - center_location) * delta;
TFSF_y_tri_T = -(TY - center_height) * delta;

% TFSF_x_tri_B = (BX - center_location) * delta;
% TFSF_y_tri_B = -(BY - center_height) * delta;
% 
% TFSF_y_tri_L = -(LY - center_height) * delta;
% TFSF_x_tri_L = (LX - center_location) * delta;
% 
% TFSF_y_tri_R = -(RY - center_height) * delta;
% TFSF_x_tri_R = (RX - center_location) * delta;



%%%%%%%%%%%%%%%%%%%% TFSF triangles geometry (distance, phase)

TFSF_dist_tri_T = sqrt((TFSF_x_tri_T-obs_x).^2+(TFSF_y_tri_T-obs_y).^2);
TFSF_phi_tri_T = atan2(-TFSF_y_tri_T+obs_y,-TFSF_x_tri_T+obs_x);

% TFSF_dist_tri_L = sqrt((TFSF_x_tri_L-obs_x).^2+(TFSF_y_tri_L-obs_y).^2);
% TFSF_phi_tri_L = atan2(-TFSF_y_tri_L+obs_y,-TFSF_x_tri_L+obs_x);
% 
% TFSF_dist_tri_R = sqrt((TFSF_x_tri_R-obs_x).^2+(TFSF_y_tri_R-obs_y).^2);
% TFSF_phi_tri_R = atan2(-TFSF_y_tri_R+obs_y,-TFSF_x_tri_R+obs_x);
% 
% TFSF_dist_tri_B = sqrt((TFSF_x_tri_B-obs_x).^2+(TFSF_y_tri_B-obs_y).^2);
% TFSF_phi_tri_B = atan2(-TFSF_y_tri_B+obs_y,-TFSF_x_tri_B+obs_x);


%%%%%%%%%%%%%%%%%%%%%%% place-wave injection


% TFSF_dist_tri_T = TFSF_dist_tri_T .* sin(TFSF_phi_tri_T);
vect_y=TFSF_y_tri_T-obs_y;
vect_x=TFSF_x_tri_T-obs_x;
veco_y=(0-obs_y)*ones(size(TFSF_y_tri_T));
veco_x=(0-obs_x)*ones(size(TFSF_x_tri_T));
veco_abs=sqrt(veco_y.^2+veco_x.^2);
veco_y=veco_y./veco_abs;veco_x=veco_x./veco_abs;
if plane_wave==1
    TFSF_dist_tri_T=veco_x.*vect_x+veco_y.*vect_y;
else
    TFSF_dist_tri_T=vect_x.*vect_x+vect_y.*vect_y;
    TFSF_dist_tri_T=sqrt(TFSF_dist_tri_T);
end






% TFSF_dist_tri_L = TFSF_dist_tri_L .* sin(TFSF_phi_tri_L);
% 
% TFSF_dist_tri_R = TFSF_dist_tri_R .* sin(TFSF_phi_tri_R);
% 
% TFSF_dist_tri_B = TFSF_dist_tri_B .* sin(TFSF_phi_tri_B);


%%%%%%%%%%%%%%%% MISC variables


% dist0 = sqrt((obs_x)^2+(obs_y)^2);

dist0=min(min(TFSF_dist_tri_T));

numt0 = dist0/c/deltat;

eta=sqrt(mu0/epsilon0);

phi0 = atan2(obs_y,obs_x);

%%%%%%%%%%%%%%%%%%% TFSF time lag to triangles

TFSF_lag_tri_T = TFSF_dist_tri_T/c/deltat;

% TFSF_lag_tri_L = TFSF_dist_tri_L/c/deltat;
% 
% TFSF_lag_tri_R = TFSF_dist_tri_R/c/deltat;
% 
% TFSF_lag_tri_B = TFSF_dist_tri_B/c/deltat;


%%%%%%%%%%%%%%%%%%%% TFSF time loop


for n= 1:time_tot
    
    t_temp = numt0 + n;
    
    %%%%%%%%%%%%%%%%% E field injection
    
    Ez_inc_T = sin(((2.*pi.*c)./(N_lambda_signal.*delta)).*(t_temp-TFSF_lag_tri_T-time_shift).*deltat)...
        .*BW/2.*(sinc(BW*(t_temp-TFSF_lag_tri_T-time_shift).*deltat)+1/2*sinc(BW*(t_temp-TFSF_lag_tri_T-time_shift).*deltat-1.0)+1/2*sinc(BW*(t_temp-TFSF_lag_tri_T-time_shift).*deltat+1.0));
    
    Ez_inc_12_T = sin(((2.*pi.*c)./(N_lambda_signal.*delta)).*(t_temp+0.5-TFSF_lag_tri_T-time_shift).*deltat)...
        .*BW/2.*(sinc(BW*(t_temp+0.5-TFSF_lag_tri_T-time_shift).*deltat)+1/2*sinc(BW*(t_temp+0.5-TFSF_lag_tri_T-time_shift).*deltat-1.0)+1/2*sinc(BW*(t_temp+0.5-TFSF_lag_tri_T-time_shift).*deltat+1.0));
    
    
%     Ez_inc_L = sin(((2.*pi.*c)./(N_lambda_signal.*delta)).*(t_temp-TFSF_lag_tri_L-time_shift).*deltat).*...
%         exp(-(((t_temp-TFSF_lag_tri_L-time_shift).*deltat).^2)./(0.966./BW)^2);
%     
%     Ez_inc_12_L = sin(((2.*pi.*c)./(N_lambda_signal.*delta)).*(t_temp+0.5-TFSF_lag_tri_L-time_shift).*deltat).*...
%         exp(-(((t_temp+0.5-TFSF_lag_tri_L-time_shift).*deltat).^2)./(0.966./BW)^2);
%     
%     Ez_inc_R = sin(((2.*pi.*c)./(N_lambda_signal.*delta)).*(t_temp-TFSF_lag_tri_R-time_shift).*deltat).*...
%         exp(-(((t_temp-TFSF_lag_tri_R-time_shift).*deltat).^2)./(0.966./BW)^2);
%     
%     Ez_inc_12_R = sin(((2.*pi.*c)./(N_lambda_signal.*delta)).*(t_temp+0.5-TFSF_lag_tri_R-time_shift).*deltat).*...
%         exp(-(((t_temp+0.5-TFSF_lag_tri_R-time_shift).*deltat).^2)./(0.966./BW)^2);
%     
%     Ez_inc_B = sin(((2.*pi.*c)./(N_lambda_signal.*delta)).*(t_temp-TFSF_lag_tri_B-time_shift).*deltat).*...
%         exp(-(((t_temp-TFSF_lag_tri_B-time_shift).*deltat).^2)./(0.966./BW)^2);
%     
%     Ez_inc_12_B = sin(((2.*pi.*c)./(N_lambda_signal.*delta)).*(t_temp+0.5-TFSF_lag_tri_B-time_shift).*deltat).*...
%         exp(-(((t_temp+0.5-TFSF_lag_tri_B-time_shift).*deltat).^2)./(0.966./BW)^2);
        
    Ez_inc_T = Ez_inc_T ./(BW/2)*10*sqrt(120*pi*10/4/pi)/323e3;
    
    Ez_inc_12_T = Ez_inc_12_T ./(BW/2)*10*sqrt(120*pi*10/4/pi)/323e3;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% H field injection

    H_inc_T = Ez_inc_T / eta;
%     H_inc_L = Ez_inc_L / eta;
%     H_inc_R = Ez_inc_R / eta;
%     H_inc_B = Ez_inc_B / eta;
    
    
    Hx_inc_T = - H_inc_T .* sin(phi0);
%     Hx_inc_L = - H_inc_L .* sin(phi0);
%     Hx_inc_R = - H_inc_R .* sin(phi0);
%     Hx_inc_B = - H_inc_B .* sin(phi0);
    
    Hy_inc_T = H_inc_T .* cos(phi0);
%     Hy_inc_L = H_inc_L .* cos(phi0);
%     Hy_inc_R = H_inc_R .* cos(phi0);
%     Hy_inc_B = H_inc_B .* cos(phi0);
    
    H_inc_12_T = Ez_inc_12_T / eta;
%     H_inc_12_L = Ez_inc_12_L / eta;
%     H_inc_12_R = Ez_inc_12_R / eta;
%     H_inc_12_B = Ez_inc_12_B / eta;
    
    Hx_inc_12_T = - H_inc_12_T .* sin(phi0);
%     Hx_inc_12_L = - H_inc_12_L .* sin(phi0);
%     Hx_inc_12_R = - H_inc_12_R .* sin(phi0);
%     Hx_inc_12_B = - H_inc_12_B .* sin(phi0);
    
    Hy_inc_12_T = H_inc_12_T .* cos(phi0);
%     Hy_inc_12_L = H_inc_12_L .* cos(phi0);
%     Hy_inc_12_R = H_inc_12_R .* cos(phi0);
%     Hy_inc_12_B = H_inc_12_B .* cos(phi0);
    
    
    %%%%%%%%%%%%%%%%%%%%%%% incident fields matrix
    
    Ez_INC=zeros(size(zeta));
    Hx_INC_12=zeros(size(zeta));
    Hy_INC_12=zeros(size(zeta));
    
    
    Ez_INC(TFSF_height_T,TFSF_location_T)=Ez_inc_T;
%     Ez_INC(TFSF_height_B,TFSF_location_B)=Ez_inc_B;
%     Ez_INC(TFSF_height_L,TFSF_location_L)=Ez_inc_L;
%     Ez_INC(TFSF_height_R,TFSF_location_R)=Ez_inc_R;
    
    Hx_INC_12(TFSF_height_T,TFSF_location_T)=Hx_inc_12_T;
%     Hx_INC_12(TFSF_height_B,TFSF_location_B)=Hx_inc_12_B;
%     Hx_INC_12(TFSF_height_L,TFSF_location_L)=Hx_inc_12_L;
%     Hx_INC_12(TFSF_height_R,TFSF_location_R)=Hx_inc_12_R;
    
    Hy_INC_12(TFSF_height_T,TFSF_location_T)=Hy_inc_12_T;
%     Hy_INC_12(TFSF_height_B,TFSF_location_B)=Hy_inc_12_B;
%     Hy_INC_12(TFSF_height_L,TFSF_location_L)=Hy_inc_12_L;
%     Hy_INC_12(TFSF_height_R,TFSF_location_R)=Hy_inc_12_R;
    
    %%%%%%%%%%%%%%%%%%%%%%% update equations

%     Hx = Hx.*G + G2.*flipud(ifft((-Ky.*fft((flipud(Ez)),[],1)),[],1));  
%     Hx = Hx + G2 .* Ez_INC .* flipud(ifft((Ky.*fft((flipud(zeta)),[],1)),[],1));
%     
%     
%     Hy = Hy.*A + A2.*ifft((Kx.*fft((Ez),[],2)),[],2);    
%     Hy = Hy - A2 .* Ez_INC .* ifft((Kx.*fft((zeta),[],2)),[],2);
%     
%     
%     Ezx = Ezx.*E + E2.*ifft((Kx.*fft((Hy),[],2)),[],2);
%     Ezx = Ezx - E2 .* Hy_INC_12 .* ifft((Kx.*fft((zeta),[],2)),[],2);
%     
% 
%     Ezy = Ezy.*C + C2.*flipud(ifft((-Ky.*fft((flipud(Hx)),[],1)),[],1));
%     Ezy = Ezy + C2 .* Hx_INC_12 .* flipud(ifft((Ky.*fft((flipud(zeta)),[],1)),[],1));
    
    

    HxI=HxI+Hx*deltat;
    HyI=HyI+Hy*deltat;
    EzxI=EzxI+Ezx*deltat;
    EzyI=EzyI+Ezy*deltat;
    



    
    Hx = Hx.*G + G2.*flipud(ifft((-Ky.*fft((flipud(Ez)),[],1)),[],1))./s0y;    
    Hx = Hx + G2 .* Ez_INC .* flipud(ifft((Ky.*fft((flipud(zeta)),[],1)),[],1))./s0y;
%     Hx = Hx - conductivity1.*sigma_starx./mu.*G2.*HxI;
    
    
    Hy = Hy.*A + A2.*ifft((Kx.*fft((Ez),[],2)),[],2)./s0x;   
    Hy = Hy - A2 .* Ez_INC .* ifft((Kx.*fft((zeta),[],2)),[],2)./s0x;
%     Hy = Hy - conductivity1.*sigma_stary./mu.*A2.*HyI;
    
    
    Ezx = Ezx.*E + E2.*ifft((Kx.*fft((Hy),[],2)),[],2)./s0x;
    Ezx = Ezx - E2 .* Hy_INC_12 .* ifft((Kx.*fft((zeta),[],2)),[],2)./s0x;
    Ezx = Ezx - conductivity1.*sigmax_orig./epsilon1.*E2.*EzxI;
    

    Ezy = Ezy.*C + C2.*flipud(ifft((-Ky.*fft((flipud(Hx)),[],1)),[],1))./s0y;
    Ezy = Ezy + C2 .* Hx_INC_12 .* flipud(ifft((Ky.*fft((flipud(zeta)),[],1)),[],1))./s0y;
    Ezy = Ezy - conductivity1.*sigmay_orig./epsilon1.*C2.*EzyI;
    
    
        
    
    Ez = Ezx + Ezy;
    
    
    Ezt = Ez;
    Hxt = Hx;
    Hyt = Hy;
    
    
    % Collect data in time domain at specified point
    Ez_data(n,:) = Ezt(surface_height,surface_location);
%     Ez_data_L(n,:) = Ezt(surface_height_L,surface_location_L).';
%     Ez_data_R(n,:) = Ezt(surface_height_R,surface_location_R).';
%     Ez_data_B(n,:) = Ezt(surface_height_B,surface_location_B);
    
    Hx_data(n,:) = Hxt(surface_height,surface_location);
%     Hy_data_L(n,:) = Hyt(surface_height_L,surface_location_L).';
%     Hy_data_R(n,:) = Hyt(surface_height_R,surface_location_R).';
%     Hx_data_B(n,:) = Hxt(surface_height_B,surface_location_B);
    
    % Display video and write to file if necessary
%     set(h,'cdata',real(Ezt)); %getframe(ax(3));
%     set(h,'cdata',abs(Ezt).'); getframe(f1);
    set(h,'cdata',20*log10(abs(real(Ezt)))); %getframe(f1);
    title(['time: ' num2str(((n)+numt0)*deltat*1e6) ' \mus'],'FontSize',20)
    F(n) = getframe(gcf);
%     display([num2str(n) '/' num2str(time_tot)])
%     toc

end

% fig = figure;
% movie(fig,F,1)

F=F(1:10:end);

myVideo = VideoWriter(nameflag,'MPEG-4');
open(myVideo);
% writeVideo(myVideo, F);
writeVideo(myVideo, F);
close(myVideo);

% load(['./1km_N2F/DATASET_' nameflag])