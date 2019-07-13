

dep=ydim+50;
wid=xdim;


dd=-((1:ydim+50)-center_height)*delta;
ww=((1:xdim)-center_location)*delta;

tomo=zeros(dep,wid);
% tomo_abs=zeros(dep,wid);




% fftn=2^16;
fftn=2^(round(log2(time_tot*50)));


ENVELOPE=[];
PHA=[];
FY1=[];
FY2=[];
PEAK=[];
NUMT0R=[];


if length(obs_y_r)==1
    if focus=='-'
        delta_l=c/2/BW;
    else
        delta_l=focus;
    end
    r0=obs_y_r;
    L=free_space_wavelength*r0/delta_l;
    dl=free_space_wavelength*r0/wid/delta;
    Nl=round(L/dl);
    if mod(Nl,2)==0
        Nl=Nl+1;
    end
    tomo_data=zeros(Nl,fftn);
    numt0r_data=zeros(1,Nl);
    lx=(-(Nl-1)/2:(Nl-1)/2)*dl;
    ly=ones(size(lx))*r0;
    for lind=1:Nl
        obs_y_r_sub=ly(lind);obs_x_r_sub=lx(lind);
        N2F_TFSF_surface_noplot;
%         close all;
        tomo_data(lind,:)=envelope.*exp(1j*pha);
        numt0r_data(lind)=numt0r;
        ENVELOPE=[ENVELOPE;envelope];
        PHA=[PHA;pha];
        FY1=[FY1;fy1];
        FY2=[FY2;fy2];
        PEAK=[PEAK;peak];
        NUMT0R=[NUMT0R;numt0r];
    end
else
    Nl=length(obs_y_r);
    tomo_data=zeros(Nl,fftn);
    numt0r_data=zeros(1,Nl);
    lx=obs_x_r;ly=obs_y_r;
    for lind=1:Nl
        obs_y_r_sub=ly(lind);obs_x_r_sub=lx(lind);
        N2F_TFSF_surface_noplot;
%         close all;
        tomo_data(lind,:)=envelope.*exp(1j*pha);
        numt0r_data(lind)=numt0r;
        ENVELOPE=[ENVELOPE;envelope];
        PHA=[PHA;pha];
        FY1=[FY1;fy1];
        FY2=[FY2;fy2];
        PEAK=[PEAK;peak];
        NUMT0R=[NUMT0R;numt0r];
    end
end

envelope=ENVELOPE;
pha=PHA;
fy1=FY1;
fy2=FY2;
peak=PEAK;
numt0r=NUMT0R;


for dep_ind=1:dep
    for wid_ind=1:wid
        for lind=1:Nl
%             dep_ind, wid_ind, lind
            temp_delay=sqrt((lx(lind)-ww(wid_ind))^2+(ly(lind)-dd(dep_ind))^2)/c+sqrt((ly(lind)-dd(dep_ind))^2)/c;
            temp_index=round(temp_delay/deltat-numt0+time_shift-numt0r_data(lind));
            if (temp_index<=fftn)&&(temp_index>=1)
                tomo(dep_ind,wid_ind)=tomo(dep_ind,wid_ind)+tomo_data(lind,temp_index)*exp(1j*omega*temp_delay);
%                 tomo_abs(dep_ind,wid_ind)=tomo_abs(dep_ind,wid_ind)+abs(tomo_data(lind,temp_index)*exp(1j*omega*temp_delay));
            end
        end
    end
end

figure;imagesc(ww+x_dist/2,-dd+source_height*delta,abs(tomo))
% figure;imagesc(ww+x_dist/2,-dd+source_height*delta,abs(tomo)./tomo_abs)            