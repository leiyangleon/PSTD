function [epsilon1,conductivity1] = dielectric_scene_cylinder(target_height,target_location,radius,relative_permittivity,sigma,...
    epsilon1,conductivity1,delta,epsilon0,omega)

x0=target_location;

y0=target_height;

ext=round(radius/delta);

xx=(x0-ext-10):(x0+ext+10);

yy=(y0-ext-10):(y0+ext+10);

eps_r=relative_permittivity;

eps_real=real(eps_r);

eps_imag=-imag(eps_r);

for xxi=1:length(xx)
    for yyi=1:length(yy)
        if sqrt((xx(xxi)-x0).^2+(yy(yyi)-y0).^2)<ext
            epsilon1(yy(yyi),xx(xxi))=(eps_real)*epsilon0;
            conductivity1(yy(yyi),xx(xxi))=sigma+(eps_imag)*epsilon0*omega;
        end
    end
end


