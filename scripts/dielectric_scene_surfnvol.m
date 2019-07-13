function [epsilon1,conductivity1] = dielectric_scene_surfnvol(surface_position,rmsh_1D,lc_1D,type_1D,std_eps,lc_2D,type_2D,relative_permittivity,sigma,...
        epsilon1,conductivity1,delta,epsilon0,omega,wavelength);

[ydim,xdim]=size(epsilon1);

x_ext=1:xdim;

y_ext=1:ydim;

rmsh_1D=wavelength*(rmsh_1D);

lc_1D=wavelength*(lc_1D);

rho_curv=lc_1D^2/2/sqrt(3)/rmsh_1D*(1+2*rmsh_1D^2/lc_1D^2)^(3/2);

eps_r=relative_permittivity;

eps_real=real(eps_r);

eps_imag=-imag(eps_r);

mean_eps=eps_real;


%%%%%%%%%%%%%%%%%%  2D volume generation over the entire domain
h_2D=rough_volume(length(x_ext),length(y_ext),length(x_ext)*delta,length(y_ext)*delta,std_eps,lc_2D,lc_2D,type_2D);
h_2D=real(h_2D);
eps_block=mean_eps*(1+h_2D)*epsilon0;


%%%%%%%%%%%%%%%%%%  1D surface generation across the entire domain
h_1D=rough_surface(length(x_ext),length(x_ext)*delta,rmsh_1D,lc_1D,type_1D);
h_1D=round(h_1D/delta);
surface_position=surface_position+h_1D;



for col=1:length(h_1D)
    eps_block(1:surface_position(col),col)=epsilon0;
end


cond_block=(eps_block/epsilon0>1)*(sigma+(eps_imag)*epsilon0*omega);

epsilon1(y_ext,x_ext)=eps_block;

conductivity1(y_ext,x_ext)=cond_block;