dep=ydim+50;
wid=xdim;

dd=-((1:ydim+50)-center_height)*delta;
ww=((1:xdim)-center_location)*delta;

tomo=zeros(dep,wid);
% tomo_abs2=zeros(dep,wid);

% fftn=2^16;
fftn=2^(round(log2(time_tot*50)));

ENVELOPE=[];
PHA=[];
FY1=[];
FY2=[];
PEAK=[];
NUMT0R=[];
NUMT0=[];

if length(obs_y_r)==1
    if focus=='-'
        delta_l=c/2/BW;
    else
        delta_l=focus;
    end
    r0=obs_y_r;
    L=free_space_wavelength*r0/delta_l/2;
    dl=free_space_wavelength*r0/wid/delta/2;
    Nl=round(L/dl);
    if mod(Nl,2)==0
        Nl=Nl+1;
    end
    tomo_data=zeros(Nl,fftn);
    numt0_data=zeros(1,Nl);
    numt0r_data=zeros(1,Nl);
    lx=(-(Nl-1)/2:(Nl-1)/2)*dl;
    ly=ones(size(lx))*r0;
    for lind=1:Nl
        obs_y=ly(lind);obs_x=lx(lind);
        if sigType==0
            TFSF_surface_CR_GPML;
            Ez=zeros(ydim,xdim);Ezx=zeros(ydim,xdim);Ezy=zeros(ydim,xdim);Hy=zeros(ydim,xdim);Hx=zeros(ydim,xdim);
        else
            error('Waveform NOT defined!!!')
        end
        obs_y_r_sub=obs_y;obs_x_r_sub=obs_x;
        N2F_TFSF_surface_noplot;
    %     close all;
        tomo_data(lind,:)=envelope.*exp(1j*pha);
        numt0_data(lind)=numt0;
        numt0r_data(lind)=numt0r;
        ENVELOPE=[ENVELOPE;envelope];
        PHA=[PHA;pha];
        FY1=[FY1;fy1];
        FY2=[FY2;fy2];
        PEAK=[PEAK;peak];
        NUMT0R=[NUMT0R;numt0r];
        NUMT0=[NUMT0;numt0];
    end
else
    Nl=length(obs_y_r);
    tomo_data=zeros(Nl,fftn);
    numt0_data=zeros(1,Nl);
    numt0r_data=zeros(1,Nl);
    lx=obs_x_r;ly=obs_y_r;
    for lind=1:Nl
        obs_y=ly(lind);obs_x=lx(lind);
        if sigType==0
            TFSF_surface_CR_GPML;
            Ez=zeros(ydim,xdim);Ezx=zeros(ydim,xdim);Ezy=zeros(ydim,xdim);Hy=zeros(ydim,xdim);Hx=zeros(ydim,xdim);
        else
            error('Waveform NOT defined!!!')
        end
        obs_y_r_sub=obs_y;obs_x_r_sub=obs_x;
        N2F_TFSF_surface_noplot;
    %     close all;
        tomo_data(lind,:)=envelope.*exp(1j*pha);
        numt0_data(lind)=numt0;
        numt0r_data(lind)=numt0r;
        ENVELOPE=[ENVELOPE;envelope];
        PHA=[PHA;pha];
        FY1=[FY1;fy1];
        FY2=[FY2;fy2];
        PEAK=[PEAK;peak];
        NUMT0R=[NUMT0R;numt0r];
        NUMT0=[NUMT0;numt0];
    end
end
    

envelope=ENVELOPE;
pha=PHA;
fy1=FY1;
fy2=FY2;
peak=PEAK;
numt0r=NUMT0R;
numt0=NUMT0;


for dep_ind=1:dep
    for wid_ind=1:wid
        for lind=1:Nl
%             dep_ind, wid_ind, lind
            temp_delay=2*sqrt((lx(lind)-ww(wid_ind))^2+(ly(lind)-dd(dep_ind))^2)/c;
            temp_index=round(temp_delay/deltat-numt0_data(lind)+time_shift-numt0r_data(lind));
            if (temp_index<=fftn)&&(temp_index>=1)
                tomo(dep_ind,wid_ind)=tomo(dep_ind,wid_ind)+tomo_data(lind,temp_index)*exp(1j*omega*temp_delay);
%                 tomo_abs2(dep_ind,wid_ind)=tomo_abs2(dep_ind,wid_ind)+abs(tomo_data(lind,temp_index)*exp(1j*omega*temp_delay)).^2;
            end
        end
    end
end

figure;imagesc(ww+x_dist/2,-dd+source_height*delta,abs(tomo))
% figure;imagesc(ww+x_dist/2,-dd+source_height*delta,abs(tomo)./sqrt(tomo_abs2))         