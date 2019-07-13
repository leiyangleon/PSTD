

[mm,nn]=size(Ez_data);



dep=ydim+50-1*source_height;
wid=xdim;


radar=zeros(dep,nn);

dd=(1:dep)*delta;
tt=(obs_y+(source_height-surface_height)*delta+2*dd)/c+time_shift*deltat;
tt_index=round(tt/deltat-numt0);

tt_index(find(tt_index>mm))=mm;tt_index(find(tt_index<1))=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%  quadrature demodulation
FFTt = ((1:time_tot)+numt0)*deltat;

FFTdk = 2*pi/(FFTt(end)-FFTt(1));

FFTdknum = round(2*pi*fcenter/FFTdk);

cosfun=cos(((2.*pi.*c)./(N_lambda_signal.*delta)).*(FFTt));
sinfun=sin(((2.*pi.*c)./(N_lambda_signal.*delta)).*(FFTt));
output=Ez_data;

for ntemp=1:nn
    FFTcos=fft(Ez_data(:,ntemp)'.*cosfun);
    FFTcos(FFTdknum:length(FFTt)-FFTdknum)=0;
    FFTsin=fft(-Ez_data(:,ntemp)'.*sinfun);
    FFTsin(FFTdknum:length(FFTt)-FFTdknum)=0;
    envelope = 2*(ifft(FFTcos)+1j*ifft(FFTsin));
    output(:,ntemp)=envelope.';
    radar(:,ntemp)=output(tt_index,ntemp);
end

figure;imagesc(surface_location*delta,dd+source_height*delta,abs(radar))