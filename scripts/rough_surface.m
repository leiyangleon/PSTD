function h=rough_surface(N,L,rmsh,lc,type)
Kn=2*pi/L*[0:(N/2) -fliplr(1:(N/2-1))];
R=randn(2,N/2);

if strcmp(type,'exp')
    %%%%%%%%%%%%%%%%    exponential
    W=rmsh^2/pi*lc*ones(size(Kn))./(1+lc^2*Kn.^2);
elseif strcmp(type,'norm')
    %%%%%%%%%%%%%%%%    gaussian
    W=rmsh^2/sqrt(pi)/2*lc*exp(-lc^2*Kn.^2/4);
else
    h=0;
    error('Bad Func Type!!!')
end



bn=R(1,:)+1j*R(2,:);
bn=[real(bn(1)) bn(2:end) imag(bn(1))];
bn=[sqrt(2*pi*L*W(1))*bn(1) sqrt(pi*L*W(2:(N/2))).*bn(2:end-1) sqrt(2*pi*L*W(N/2+1))*bn(end)];

bn=[bn fliplr(conj(bn(2:end-1)))];

h=N/L*ifft(bn);