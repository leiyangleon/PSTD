%  In this script, we set up an example run of active focusing along with refraction correction using output after running PSTD.m.
%  Here, PSTD.m was run by using the array of user-specified receiver locations ("obs_x_r" and "obs_y_r" are user-defined vectors).
%  One can follow this setup strategy to obtain an example run for the passive focusing case.




%%  INPUT for the PSTD run

x_dist=200;
y_dist=200;

xdim=round(x_dist/delta);
ydim=round(y_dist/delta);

source_height=round(60/delta);
center_height=source_height;
center_location=round(x_dist/2/delta);

time_shift=round(0.3e-6/deltat);
time_tot=round(1.5e-6/deltat);

ref_corr=[-40 2.3];

omega=2*pi*20e6;
c=3e8;

%%% in addition, must declare "obs_x_r" and "obs_y_r" in MATLAB workspace as already done for running PSTD.m %%%





%%  INPUT (based on OUTPUT after running PSTD.m) needed to set up the example run

numt0_data=numt0;numt0r_data=numt0r;tomo_data=envelope.*exp(1j*pha);





%%  Direct copy specific lines from the active focusing script: active_focus_result.m

dep=ydim+50;
wid=xdim;

dd=-((1:ydim+50)-center_height)*delta;
ww=((1:xdim)-center_location)*delta;

tomo=zeros(dep,wid);

fftn=2^(round(log2(time_tot*50)));

lx=obs_x_r;ly=obs_y_r;Nl=length(obs_y_r);

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




%%  Direct copy specific lines from the active focusing with refraction correction script: active_focus_corr_result.m

tomo=zeros(dep,wid);

SUR=ref_corr(1);
n_rf=sqrt(ref_corr(2));

for dep_ind=1:dep
    for wid_ind=1:wid
        for lind=1:Nl
            if dd(dep_ind)>SUR
%             dep_ind, wid_ind, lind
                temp_delay=2*sqrt((lx(lind)-ww(wid_ind))^2+(ly(lind)-dd(dep_ind))^2)/c;
                temp_index=round(temp_delay/deltat-numt0_data(lind)+time_shift-numt0r_data(lind));
                if (temp_index<=fftn)&&(temp_index>=1)
                    tomo(dep_ind,wid_ind)=tomo(dep_ind,wid_ind)+tomo_data(lind,temp_index)*exp(1j*omega*temp_delay);
                end
            else
%                 theta_rf=(0:.001:90)/180*pi;
%                 check_rf=(ly(lind)-SUR)*tan(theta_rf)+(SUR-dd(dep_ind))*sin(theta_rf)./sqrt(n_rf^2-sin(theta_rf).^2)-abs(lx(lind)-ww(wid_ind));
%                 [min_rf,mid_rf]=min(abs(check_rf));
%                 theta_rf=theta_rf(mid_rf);
                [x_iter, l1_iter, l2_iter, rho_iter, theta_rf, tta_iter] = solveRefractionPointFlat2(ly(lind)-SUR,abs(lx(lind)-ww(wid_ind)),SUR-dd(dep_ind),n_rf^2,20);
                temp_delay=(SUR-dd(dep_ind))/cos(asin(sin(theta_rf)/n_rf))*n_rf+(ly(lind)-SUR)/cos(theta_rf);
                temp_delay=2*temp_delay/c;
                temp_index=round(temp_delay/deltat-numt0_data(lind)+time_shift-numt0r_data(lind));
                if (temp_index<=fftn)&&(temp_index>=1)
                    tomo(dep_ind,wid_ind)=tomo(dep_ind,wid_ind)+tomo_data(lind,temp_index)*exp(1j*omega*temp_delay);
                end
            end
        end
    end
end

figure;imagesc(ww+x_dist/2,-dd+source_height*delta,abs(tomo))