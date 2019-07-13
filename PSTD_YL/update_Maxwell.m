% Calculate maxwells constants at each grid point

%Magnetic conductivity matrix obtained by Perfectly Matched Layer condition
%This is also split into x and y directions in Berenger's model
sigma_starx=(sigmax.*mu)./epsilon1;
sigma_stary=(sigmay.*mu)./epsilon1;

sigmax=sigmax+conductivity1;
sigmay=sigmay+conductivity1;
epsilon=epsilon1;

%Multiplication factor matrices for H matrix update to avoid being calculated many times
%in the time update loop so as to increase computation speed
 G=((mu-0.5*deltat*sigma_stary)./(mu+0.5*deltat*sigma_stary));
 G2=deltat*ones(size(sigma_stary))./(mu+sigma_stary*deltat/2);
 A=((mu-0.5*deltat*sigma_starx)./(mu+0.5*deltat*sigma_starx));
 A2=deltat*ones(size(sigma_starx))./(mu+sigma_starx*deltat/2);
 
%Multiplication factor matrices for E matrix update to avoid being calculated many times
%in the time update loop so as to increase computation speed
 C=((epsilon-0.5*deltat*sigmay)./(epsilon+0.5*deltat*sigmay));
 C2=deltat*ones(size(sigmay))./(epsilon+sigmay*deltat/2);
 E=((epsilon-0.5*deltat*sigmax)./(epsilon+0.5*deltat*sigmax));
 E2=deltat*ones(size(sigmax))./(epsilon+sigmax*deltat/2);