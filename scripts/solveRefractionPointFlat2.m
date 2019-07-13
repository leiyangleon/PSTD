function [x, l1, l2, rho, tia, tta] = solveRefractionPointFlat2(h,L,d,er,N)


tia = atan2(L,(d+h));
for n=1:N,
    tta = asin(sin(tia)./sqrt(er));
    x = L - d.*tan(tta);
    tia = atan2(x,h);
end

l1 = sqrt(x.^2 + h.^2);
l2 = sqrt((L-x).^2 + d.^2);
rho = l1 + sqrt(er).*l2;

