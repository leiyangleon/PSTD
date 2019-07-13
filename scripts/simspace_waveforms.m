%% Required EM constants
epsilon0=(1/(36*pi))*1e-9;
c=3e+8;
mu0=4*pi*1e-7;

%% Courant stability factor
S=sqrt(2)/pi;
% S=S/3;

%% Define spatial and temporal discretization requirements
f_solution = 1*fcenter;
free_space_wavelength = c/f_solution;
delta = (free_space_wavelength/cellsperwavelength);
deltat=S*delta/c;

time_tot=round(time_tot/deltat);
time_shift=round(time_shift/deltat);

%% Define waveforms based on center frequency
wavelength = c/fcenter;
N_lambda_signal=wavelength/delta;
n_use = 1:time_tot;

%% Waveform - Use sine modulated gaussian waveform, with gaussian width as a function of center frequency
if sigType==0
    Ez_source = sin(((2.*pi.*c)./(N_lambda_signal.*delta)).*(n_use-time_shift).*deltat).*exp(-(((n_use-time_shift).*deltat).^2)./(0.966./BW)^2);
elseif sigType==1
    fc=fcenter;
    T=1.55/fc;
    a1=-0.488;
    a2=0.145;
    a3=-0.01022222;
    for i_use=1:length(n_use)
        if (((n_use(i_use)-time_shift)*deltat)<T)&&(((n_use(i_use)-time_shift)*deltat)>0)
            Ez_source(i_use)=a1*pi/T*sin(2*pi*((n_use(i_use)-time_shift).*deltat)/T)+a2*2*pi/T*sin(2*2*pi*((n_use(i_use)-time_shift).*deltat)/T)+a3*3*pi/T*sin(2*3*pi*((n_use(i_use)-time_shift).*deltat)/T);
            Ez_source(i_use)=-Ez_source(i_use)/2.757e7;
        else
            Ez_source(i_use)=0;
        end
    end
elseif sigType==2
    Ez_source = sin(((2.*pi.*c)./(N_lambda_signal.*delta)).*(n_use-time_shift).*deltat).*sinc(BW*(n_use-time_shift).*deltat);
else
    error('Waveform NOT defined!!!')
end