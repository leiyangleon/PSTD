function h=rough_volume(Nx,Ny,Lx,Ly,rmsh,lcx,lcy,type)
Knx=2*pi/Lx*[0:(Nx/2) -fliplr(1:(Nx/2-1))];
Kny=2*pi/Ly*[0:(Ny/2) -fliplr(1:(Ny/2-1))];

[KnX,KnY]=meshgrid(Knx,Kny);

if strcmp(type,'exp')
    %%%%%%%%%%%%%%%%    exponential
    W=rmsh^2/2/pi*lcx*lcy*ones(size(KnX))./(1+lcx^2*KnX.^2+lcy^2*KnY.^2).^(3/2);
elseif strcmp(type,'norm')
    %%%%%%%%%%%%%%%%    gaussian
    W=rmsh^2/4/pi*lcx*lcy*exp(-lcx^2*KnX.^2/4-lcy^2*KnY.^2/4);
else
    h=0;
    error('Bad Func Type!!!')
end

R1=randn(Ny/2-1,Nx/2-1)+1j*randn(Ny/2-1,Nx/2-1);

R2=randn(Ny/2-1,Nx/2-1)+1j*randn(Ny/2-1,Nx/2-1);

R=randn(Ny,Nx);

Bn=R.*sqrt(2*pi*Lx*2*pi*Ly*W);

Bn(2:Ny/2,2:Nx/2)=R1.*sqrt(2*pi*Lx*2*pi*Ly*W(2:Ny/2,2:Nx/2)/2);

Bn(Ny/2+2:end,Nx/2+2:end)=fliplr(flipud(conj(Bn(2:Ny/2,2:Nx/2))));


Bn(Ny/2+2:end,2:Nx/2)=R2.*sqrt(2*pi*Lx*2*pi*Ly*W(Ny/2+2:end,2:Nx/2)/2);

Bn(2:Ny/2,Nx/2+2:end)=fliplr(flipud(conj(Bn(Ny/2+2:end,2:Nx/2))));


h=Nx*Ny/Lx/Ly*ifft2(Bn);
% figure;imagesc(real(h),[-8 8]);