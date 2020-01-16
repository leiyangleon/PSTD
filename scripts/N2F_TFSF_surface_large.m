%%%%%%%%%%%         N2F

dist0r = sqrt((0*obs_x_r_sub)^2+(obs_y_r_sub)^2)-(source_height-surface_height);

numt0r = dist0r/c/deltat;






Jz = -Hx_data;
Mx = -Ez_data;

surface_x = (surface_location - center_location) * delta;
surface_y = -(surface_height - center_height) * delta * ones(size(surface_x));
% surface_y_L = -(surface_height_L - center_height) * delta;
% surface_y_R = -(surface_height_R - center_height) * delta;
% surface_x_L = (surface_location_L - center_location) * delta * ones(size(surface_y_L));
% surface_x_R = (surface_location_R - center_location) * delta * ones(size(surface_y_R));
% surface_x_B = (surface_location_B - center_location) * delta;
% surface_y_B = -(surface_height_B - center_height) * delta * ones(size(surface_x_B));


single_bound=round(400/2/delta);

num_single_bound=floor(24500/400/2*1);

complex_waveform=0;

% figure;
% kk=1;

for ii=-num_single_bound:num_single_bound
    if ii==0
        display(['aaa']);
    end
%     ii
    low_bound=center_location-single_bound+ii*single_bound*2;
    up_bound=center_location+single_bound+ii*single_bound*2;
    cen_ind=center_location+ii*single_bound*2;surface_x_cen=surface_x(cen_ind);
    surface_x=surface_x(low_bound:up_bound);surface_y=surface_y(low_bound:up_bound);Jz=Jz(:,low_bound:up_bound);Mx=Mx(:,low_bound:up_bound);
%     surface_x=surface_x-ii*single_bound*2*delta;
    surface_x=surface_x-surface_x_cen;
    [m,n]=size(Jz);



    % fftn=2^16;
    fftn=2^(round(log2(time_tot*10)));


    FJz=fft(Jz,fftn,1);
    % FEz_L=fft(Ez_data_L,fftn,1);
    % FEz_R=fft(Ez_data_R,fftn,1);
    % FEz_B=fft(Ez_data_B,fftn,1);

    FMx=fft(Mx,fftn,1);
    % FHy_L=fft(Hy_data_L,fftn,1);
    % FHy_R=fft(Hy_data_R,fftn,1);
    % FHx_B=fft(Hx_data_B,fftn,1);


    dkt = 2*pi/deltat/fftn;

    clear kt;


    for i=1:floor(fftn/2)
        kt(i)=(i-1)*dkt; 
        kt(floor(fftn/2)+i)=-pi/deltat+(i-1)*dkt; 
    end



    kt = kt * sqrt(mu0*epsilon0);

%     phi=atan2(obs_y_r_sub,-(obs_x_r_sub+ii*single_bound*2*delta));
    phi=atan2(obs_y_r_sub,-(surface_x_cen));
    
    
%     win_fun=blackmanharris(size(FJz,2))';
% 
% %     win_fun=tukeywin(size(FJz,2),1)';
% 
%     ww = ones(size(FJz,1),1)*(win_fun);
% 
%     FJz=FJz.*ww;
%     FMx=FMx.*ww;
    
    

    FNe=sum((-1)*delta*FJz.*exp(1j*kt'*(surface_x*cos(phi)+surface_y*sin(phi))),2);
    % FNe_L=sum(delta*FHy_L.*exp(1j*kt'*(surface_x_L*cos(phi)+surface_y_L*sin(phi))),2);
    % FNe_R=sum((-1)*delta*FHy_R.*exp(1j*kt'*(surface_x_R*cos(phi)+surface_y_R*sin(phi))),2);
    % FNe_B=sum((-1)*delta*FHx_B.*exp(1j*kt'*(surface_x_B*cos(phi)+surface_y_B*sin(phi))),2);

    FLi=sum(-sin(phi)*delta*FMx.*exp(1j*kt'*(surface_x*cos(phi)+surface_y*sin(phi))),2);
    % FLi_L=sum(-cos(phi)*delta*FEz_L.*exp(1j*kt'*(surface_x_L*cos(phi)+surface_y_L*sin(phi))),2);
    % FLi_R=sum(cos(phi)*delta*FEz_R.*exp(1j*kt'*(surface_x_R*cos(phi)+surface_y_R*sin(phi))),2);
    % FLi_B=sum(-sin(phi)*delta*FEz_B.*exp(1j*kt'*(surface_x_B*cos(phi)+surface_y_B*sin(phi))),2);


    % FNe=FNe+FNe_L+FNe_R;
    % FLi=FLi+FLi_L+FLi_R;

    % FNe=FNe+FNe_L+FNe_R+FNe_B;
    % FLi=FLi+FLi_L+FLi_R+FLi_B;

%     r=sqrt((obs_x_r_sub+ii*single_bound*2*delta)^2+obs_y_r_sub^2);
    r=sqrt((surface_x_cen)^2+obs_y_r_sub^2);

    eta=sqrt(mu0/epsilon0);

    FEe=-sqrt(1j/(2*4*pi*r)*kt').*exp(-1j*kt'*r).*(FLi+eta*FNe);

    FEe=FEe.*exp(1j*kt'*dist0r);


    FHi=FEe/eta;

    Ee=ifft(FEe,fftn,1);

    Hi=ifft(FHi,fftn,1);



%     figure;
%     axx(1)=subplot(2,2,1);plot(((1:fftn)+numt0+numt0r)*deltat*1e6,real(-Ee))
%     title('Far field','FontSize',20);xlabel('Time / \mus','FontSize',20);grid on;set(gca,'fontsize',20)

%     axx(2)=subplot(2,2,2);plot(((1:time_tot)+numt0)*deltat*1e6,real(Ez_data(:,round(size(Ez_data,2)/2))))
%     title('Near field','FontSize',20);xlabel('Time / \mus','FontSize',20);grid on;set(gca,'fontsize',20)

    %%%%%%%%%%%%%%%%%%%%%%%%%%  quadrature demodulation
    FFTt = ((1:fftn)+numt0+numt0r)*deltat;
    FFTdk = 2*pi/(FFTt(end)-FFTt(1));
    FFTdknum = round(2*pi*fcenter/FFTdk);
    cosfun=cos(((2.*pi.*c)./(N_lambda_signal.*delta)).*(FFTt));
    sinfun=sin(((2.*pi.*c)./(N_lambda_signal.*delta)).*(FFTt));
    output = real(-Ee)';
    FFTcos=fft(output.*cosfun);
    FFTcos(FFTdknum:length(FFTt)-FFTdknum)=0;
    FFTsin=fft(-output.*sinfun);
    FFTsin(FFTdknum:length(FFTt)-FFTdknum)=0;
    envelope = 2*abs(ifft(FFTcos)+1j*ifft(FFTsin));
%     axx(3)=subplot(2,2,3);plot(FFTt*1e6,envelope);
%     title('Amplitude','FontSize',20);xlabel('Time / \mus','FontSize',20);grid on;set(gca,'fontsize',20)
    pha = angle(ifft(FFTcos)+1j*ifft(FFTsin));
%     axx(4)=subplot(2,2,4);plot(FFTt*1e6,pha);
%     title('Phase','FontSize',20);xlabel('Time / \mus','FontSize',20);grid on;set(gca,'fontsize',20)
    % linkaxes(axx,'x');

    fy1=abs(FEe).';
    FEz_source=fft(Ez_source,fftn);
    fy2=abs(FEz_source);
    ind=round(fcenter/(FFTdk/2/pi)+1);
    flag=1;
    if flag==1
        freq_bound=0;
    elseif flag==2
        freq_bound=10;
    else
        freq_bound=30;
    end
    peak=sum(fy1(ind-freq_bound:ind+freq_bound))/sum(fy2(ind-freq_bound:ind+freq_bound));
    
%     subplot(7,9,kk);plot(FFTt*1e6,envelope);axis([2000 2002 0 1e-5])%axis([2000 2002 -pi pi]);%
%     kk=kk+1;
    
    
    complex_waveform=complex_waveform+envelope.*exp(1j*pha);
    
    
    Jz = -Hx_data;
    Mx = -Ez_data;

    surface_x = (surface_location - center_location) * delta;
    surface_y = -(surface_height - center_height) * delta * ones(size(surface_x));
end

% time_cutoff=2*r/3e8*1e6+time_shift*deltat*1e6;
% 
% n_cutoff=round((time_cutoff-FFTt(1)*1e6)/(FFTt(2)*1e6-FFTt(1)*1e6)+1);

% complex_waveform(n_cutoff:end)=0;

% figure;plot(FFTt*1e6,abs(complex_waveform));

envelope=abs(complex_waveform);
pha=angle(complex_waveform);
