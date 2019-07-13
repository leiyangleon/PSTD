% Spatial grid 
xdim = round(x_dist/delta);
ydim = round(y_dist/delta);

% Check for odd number and add 1 to avoid errors
if mod(xdim,2) == 1
    xdim = xdim + 1;
end
if mod(ydim,2) == 1
    ydim = ydim + 1;
end

% Define original waveform - a cosine modulated gaussian pulse
omega = 2*pi*fcenter;

% Define epsilon and conductivity and mu
mu = mu0*ones(ydim,xdim);
epsilon = epsilon0*ones(ydim,xdim);
conductivity = zeros(ydim,xdim);

% % PML thickness as a function of number of wavelengths
if PML_type==0
    bound_width_x=20;
elseif PML_type==1
    bound_width_x=40;
else
    error('PML type NOT defined!!!')
end
bound_width_y = bound_width_x;

% Initialization of field matrices
Ez=zeros(ydim,xdim);
Ezx=zeros(ydim,xdim);
Ezy=zeros(ydim,xdim);
Hy=zeros(ydim,xdim);
Hx=zeros(ydim,xdim);
