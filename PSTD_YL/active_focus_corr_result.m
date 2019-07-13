
% dep=ydim+50-0*source_height-0*bound_width_x;
% wid=xdim-0*bound_width_x;
 
% dd=-((0*source_height+1:ydim+50-0*bound_width_x)-center_height)*delta;
% ww=((0*bound_width_x+1:xdim-0*bound_width_x)-center_location)*delta;
% 
tomo=zeros(dep,wid);


% delta_l=15;
% r0=obs_y;
% L=free_space_wavelength*r0/delta_l;
% dl=free_space_wavelength*r0/wid/delta;
% 
% Nl=round(L/dl);
% if mod(Nl,2)==0
%     Nl=Nl+1;
% end
% 
% lx=(-(Nl-1)/2:(Nl-1)/2)*dl;
% ly=ones(size(lx))*r0;
% 
% fftn=2^16;
% tomo_data=zeros(Nl,fftn);
% 
% for lind=1:Nl
%     obs_y_r=ly(lind);obs_x_r=lx(lind);
%     N2F_TFSF_surface;
%     close all;
%     tomo_data(lind,:)=envelope.*exp(1j*pha);
% end


% % SUR=-40;
% SUR=-30;
% n_rf=sqrt(2);
% % n_rf=sqrt(2.3);

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
