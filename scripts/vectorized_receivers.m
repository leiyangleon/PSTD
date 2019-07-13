



% fftn=2^16;
fftn=2^(round(log2(time_tot*50)));


ENVELOPE=[];
PHA=[];
FY1=[];
FY2=[];
PEAK=[];
NUMT0R=[];



Nl=length(obs_y_r);


for lind=1:Nl
    obs_y_r_sub=obs_y_r(lind);obs_x_r_sub=obs_x_r(lind);
    N2F_TFSF;
    ENVELOPE=[ENVELOPE;envelope];
    PHA=[PHA;pha];
    FY1=[FY1;fy1];
    FY2=[FY2;fy2];
    PEAK=[PEAK;peak];
    NUMT0R=[NUMT0R;numt0r];
end


envelope=ENVELOPE;
pha=PHA;
fy1=FY1;
fy2=FY2;
peak=PEAK;
numt0r=NUMT0R;

         