% Definitions (wave-vector and constants) for PSTD
kmax = pi/delta;
dk = kmax/(xdim/2);


for i=1:xdim/2
    kx(i)=(i-1)*dk; 
    kx(xdim/2+i)=-kmax+(i-1)*dk; 
end

kx=sqrt(-1)*kx;
Kx = ones(ydim,1)*kx;

kmax = pi/delta;
dk = kmax/(ydim/2);


for i=1:ydim/2
    ky(i)=(i-1)*dk; 
    ky(ydim/2+i)=-kmax+(i-1)*dk; 
end

ky=sqrt(-1)*ky;
Ky = ky.'*ones(1,xdim);

disp('Starting PSTD...');
