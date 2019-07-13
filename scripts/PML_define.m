% Define PML

% Order of polynomial on which sigma is modeled
gradingorder = 3; 
% Required reflection co-efficient
refl_coeff=1e-16; 


% 2D FDTD update for PML boundary as used in previous program
% Initializing electric conductivity matrices in x and y directions
sigmax=conductivity;
sigmay=conductivity;

% Perfectly matched layer boundary design

% Polynomial model for sigma
sigmamax_y=(-log(refl_coeff)*(gradingorder+1)*epsilon0*c)/(2*bound_width_y*delta); % check log vs log10
% boundfact1=((epsilon(bound_width_y,xdim/2)/epsilon0)*sigmamax_y)/((bound_width_y^gradingorder)*(gradingorder+1));
% boundfact2=((epsilon(ydim-bound_width_y,xdim/2)/epsilon0)*sigmamax_y)/((bound_width_y^gradingorder)*(gradingorder+1));

sigmamax_x=(-log(refl_coeff)*(gradingorder+1)*epsilon0*c)/(2*bound_width_x*delta); % check log vs log10
% boundfact3=((epsilon(ydim/2,bound_width_x)/epsilon0)*sigmamax_x)/((bound_width_x^gradingorder)*(gradingorder+1));
% boundfact4=((epsilon(ydim/2,xdim-bound_width_x)/epsilon0)*sigmamax_x)/((bound_width_x^gradingorder)*(gradingorder+1));

y=0:1:bound_width_y;
for i=1:1:xdim
    boundfact1=((epsilon1(bound_width_y+1,i)/epsilon0)*sigmamax_y)/((bound_width_y^gradingorder)*(gradingorder+1));
    boundfact2=((epsilon1(ydim-bound_width_y,i)/epsilon0)*sigmamax_y)/((bound_width_y^gradingorder)*(gradingorder+1));
    sigmay(bound_width_y+1:-1:1,i)=boundfact1*((y+0.5*ones(1,bound_width_y+1)).^(gradingorder+1)-(y-0.5*[0 ones(1,bound_width_y)]).^(gradingorder+1));
    sigmay(ydim-bound_width_y:1:ydim,i)=boundfact2*((y+0.5*ones(1,bound_width_y+1)).^(gradingorder+1)-(y-0.5*[0 ones(1,bound_width_y)]).^(gradingorder+1));
    epsilon1(bound_width_y+1:-1:1,i)=epsilon1(bound_width_y+1,i);
    epsilon1(ydim-bound_width_y:1:ydim,i)=epsilon1(ydim-bound_width_y,i);
end

x=0:1:bound_width_x;
for i=1:1:ydim
    boundfact3=((epsilon1(i,bound_width_x)/epsilon0)*sigmamax_x)/((bound_width_x^gradingorder)*(gradingorder+1));
    boundfact4=((epsilon1(i,xdim-bound_width_x)/epsilon0)*sigmamax_x)/((bound_width_x^gradingorder)*(gradingorder+1));
    sigmax(i,bound_width_x+1:-1:1)=boundfact3*((x+0.5*ones(1,bound_width_x+1)).^(gradingorder+1)-(x-0.5*[0 ones(1,bound_width_x)]).^(gradingorder+1))';
    sigmax(i,xdim-bound_width_x:1:xdim)=boundfact4*((x+0.5*ones(1,bound_width_x+1)).^(gradingorder+1)-(x-0.5*[0 ones(1,bound_width_x)]).^(gradingorder+1))';
    epsilon1(i,bound_width_x+1:-1:1)=epsilon1(i,bound_width_x+1);
    epsilon1(i,xdim-bound_width_x:1:xdim)=epsilon1(i,xdim-bound_width_x);
end
