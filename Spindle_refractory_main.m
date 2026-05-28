 clear all; clc

T = 60*1000 ;%
dt = 1e-1;
time_array = 0:dt:T;
vec_len = length(time_array);

%%
Cm=1;
taot=20;taor=20; %ms
Qt_max=400e-3;Qr_max=400e-3;
theta=-58.5;
sigma_t=6;sigma_r=6;
gammae=70e-3;gammag=58.6e-3;gammar=100e-3;
v=120e-3;
%Npp=120;Nip=72;Npi=90;Nii=90;Ntp=2.6;Nrp=2.6;Nrt=3;Ntr=5;Nrr=19;Npt=2.5;Nit=2.5;
Npp=115;Nip=72;Npi=90;Nii=90;Ntp=2.6;Nrp=2.6;Nrt=3;Ntr=5;Nrr=25;Npt=2.5;Nit=2.5;
WAMPA=1;WGABA=1;
g_T_t=3;g_T_r=2.3;
E_L_p=-64;E_L_i=-64;E_L_t=-70;E_L_r=-70;
E_K=-100;E_Ca=120;EAMPA=0;EGABA=-70;
alphaNa=2;tao_Na=1.7;R_pump=0.09;Na0=9.5;alphaCa=-51.8e-6;tao_Ca=10;
Ca0=2.4e-4;
k2=4e-4;k3=1e-1;k4=1e-3;
np=4;ginc=2;

%
k1=2.5e7;
E_h=-40;
C=(pi./sqrt(3));

g_L_K=0.019;
%g_h=0.062;
HA_max=0.03;
max_ha=10;%
min_ha=2.1;
EC_50=0.02;
b=1.6;
sigma_h=3;
g_h1=0.061;
K_TC_h=0.04;
%g_AHP=6;

alphaAHP=48;
bataAHP=0.09;


%% 
%phi_T_sd=10e-3;
% phi_T_sd=10e-3;
% 
% noise3 = normrnd(0,phi_T_sd,size(time_array)); 

%noise3 = gammae^2*sqrt(phi_T_sd*dt)*randn(1,vec_len);

%alpha=1.2;beta=0.6;D1=0.05;gamma1=D1^(1/alpha);delta=0;D2=0.005;gamma2=D2^(1/alpha);
%noise3 = levy(alpha,beta,gamma2,delta,vec_len);

noise3 =zeros(1,vec_len);%

%% 
Vp = zeros(1,vec_len);
Vi = zeros(1,vec_len);
Vt = zeros(1,vec_len);
Vr = zeros(1,vec_len);
Sep = zeros(1,vec_len);
Sepdot = zeros(1,vec_len);
Sei = zeros(1,vec_len);
Seidot = zeros(1,vec_len);
Set = zeros(1,vec_len);
Setdot = zeros(1,vec_len);
Ser = zeros(1,vec_len);
Serdot = zeros(1,vec_len);
Sgp = zeros(1,vec_len);
Sgpdot = zeros(1,vec_len);
Sgi = zeros(1,vec_len);
Sgidot = zeros(1,vec_len);
Srt = zeros(1,vec_len);
Srtdot = zeros(1,vec_len);
Srr = zeros(1,vec_len);
Srrdot = zeros(1,vec_len);
phip = zeros(1,vec_len);
phipdot = zeros(1,vec_len);
phit = zeros(1,vec_len);
phitdot = zeros(1,vec_len);
h_T_t = zeros(1,vec_len);
h_T_r = zeros(1,vec_len);
m_h1 = zeros(1,vec_len);
m_h2 = zeros(1,vec_len);
concentration_Ca = zeros(1,vec_len);
concentration_Na = zeros(1,vec_len);
R_AHP=zeros(1,vec_len);
HA=zeros(1,vec_len);
R_h=zeros(1,vec_len);
s_AHP=zeros(1,vec_len);
I=zeros(1,vec_len);
g=zeros(1,vec_len);


Vp(1) = -64;
Vi(1) = -64;
Vt(1) = -70;
Vr(1) = -70;
concentration_Ca(1) = 2.4E-4;
concentration_Na(1) =Na0;

%% 

for t = 2:vec_len
    
    m_limit_t=1/(1+exp(-(Vt(t-1)+59)/6.2));
    m_limit_r=1/(1+exp(-(Vr(t-1)+52)/7.4));
    h_limit_t=1/(1+exp((Vt(t-1)+81)/4));
    h_limit_r=1/(1+exp((Vr(t-1)+80)/5));
    
    tao_h_t=(30.8+(211.4+exp((Vt(t-1)+115.2)/5))/(1+exp((Vt(t-1)+86)/3.2)))/3^(1.2);
    tao_h_r=(85+1./(exp((Vr(t-1)+48)/4)+exp(-(Vr(t-1)+407)/50)))/3^(1.2);
    m_limit_h=1./(1+exp((Vt(t-1)+75)/5.5));
    tao_m_h=(20+1000./(exp((Vt(t-1)+71.5)/14.2)+exp(-(Vt(t-1)+89)/11.6)));
    Ph=k1*(concentration_Ca(t-1))^np./(k1*(concentration_Ca(t-1))^np+k2);
    
    Na_pump=R_pump*((concentration_Na(t-1)^3./(concentration_Na(t-1)^3+3375))-(Na0^3./(Na0^3+3375)));
   
    Qt=Qt_max./(1+exp(-C*(Vt(t-1)-theta)/(sigma_t)));
    Qr=Qr_max./(1+exp(-C*(Vr(t-1)-theta)/(sigma_r)));
    
    J_L_p=(Vp(t-1)-E_L_p);
    J_L_i=(Vi(t-1)-E_L_i);
    J_L_t=(Vt(t-1)-E_L_t);
    J_L_r=(Vr(t-1)-E_L_r);
    JAMPA_Sep=WAMPA*Sep(t-1)*(Vp(t-1)-EAMPA);
    JAMPA_Sei=WAMPA*Sei(t-1)*(Vi(t-1)-EAMPA);
    JAMPA_Set=WAMPA*Set(t-1)*(Vt(t-1)-EAMPA);
    JAMPA_Ser=WAMPA*Ser(t-1)*(Vr(t-1)-EAMPA);
    JGABA_Sgp=WGABA*Sgp(t-1)*(Vp(t-1)-EGABA);
    JGABA_Sgi=WGABA*Sgi(t-1)*(Vi(t-1)-EGABA);
    JGABA_Srt=WGABA*Srt(t-1)*(Vt(t-1)-EGABA);
    JGABA_Srr=WGABA*Srr(t-1)*(Vr(t-1)-EGABA);
    
    I_LK_t=g_L_K*(Vt(t-1)-E_K);
    I_LK_r=g_L_K*(Vr(t-1)-E_K);
    I_T_t=g_T_t*m_limit_t*m_limit_t*h_T_t(t-1).*(Vt(t-1)-E_Ca);
    I_T_r=g_T_r*m_limit_r*m_limit_r*h_T_r(t-1)*(Vr(t-1)-E_Ca);
    %I_h=g_h*(m_h1(t-1)+ginc*m_h2(t-1))*(Vt(t-1)-E_h);
    concentration_Ca(t)=concentration_Ca(t-1)+dt*(alphaCa*I_T_t-(concentration_Ca(t-1)-Ca0)./tao_Ca);
   
    %gAHP
    HA=HA_max/(1+exp(-(Vt(t-1)+80)/sigma_h));
    g_AHP=((max_ha-min_ha)/(1+(HA/EC_50))^b)+min_ha;
    g_h=g_h1+K_TC_h*HA;
    g(t)=g_h;
    %%%
    s_AHP(t)=s_AHP(t-1)+dt*(alphaAHP*(concentration_Ca(t)^2)*(1-s_AHP(t-1))-bataAHP*s_AHP(t-1));
    I_AHP=g_AHP*s_AHP(t-1)*(Vt(t-1)-E_K);

    %I_AHP=g_AHP*(concentration_Ca(t)./(30+concentration_Ca(t)))*(Vt(t-1)-E_K);
    I_h=g_h*(m_h1(t-1)+ginc*m_h2(t-1))*(Vt(t-1)-E_h);
    I(t)=I_h;

    %Vt(t)=Vt(t-1)+dt*((-J_L_t-JAMPA_Set-JGABA_Srt-Cm*taot*(I_LK_t+I_T_t+I_h+I_AHP))./taot);
    %Vr(t)=Vr(t-1)+dt*((-J_L_r-JAMPA_Ser-JGABA_Srr-Cm*taor*(I_LK_r+I_T_r))./taor);
    Vt(t)=Vt(t-1)+dt*((-J_L_t-JAMPA_Set-JGABA_Srt-Cm*taot*(I_LK_t+I_T_t+I_h))./taot);
    Vr(t)=Vr(t-1)+dt*((-J_L_r-JAMPA_Ser-JGABA_Srr-Cm*taor*(I_LK_r+I_T_r+I_AHP))./taor);
    R_AHP(t)=I_AHP;
      R_h(t)=I_h;


    Set(t) = Set(t-1) + Setdot(t-1)*dt;
    Setdot(t) = Setdot(t-1) + dt*(gammae^2*(noise3(t-1)-Set(t-1))-2*gammae*Setdot(t-1) );
    %Setdot(t) = Setdot(t-1) + dt*(gammae^2*(-Set(t-1))-2*gammae*Setdot(t-1) )+noise3(t-1);
    Ser(t) = Ser(t-1) + Serdot(t-1)*dt;
    Serdot(t) = Serdot(t-1) + dt*(gammae^2*(Nrt*Qt-Ser(t-1))-2*gammae*Serdot(t-1) );
    Srt(t) = Srt(t-1) +Srtdot(t-1)*dt;
    Srtdot(t) = Srtdot(t-1) + dt*(gammar^2*(Ntr*Qr-Srt(t-1))-2*gammar*Srtdot(t-1) );
    
    Srr(t) = Srr(t-1) +Srrdot(t-1)*dt;
    Srrdot(t) = Srrdot(t-1) + dt*(gammar^2*(Nrr*Qr-Srr(t-1))-2*gammar*Srrdot(t-1) );
    
    h_T_t(t)=h_T_t(t-1)+dt*((h_limit_t-h_T_t(t-1))./tao_h_t);
    h_T_r(t)=h_T_r(t-1)+dt*((h_limit_r-h_T_r(t-1))./tao_h_r);

    m_h1(t)=m_h1(t-1)+dt*((m_limit_h*(1-m_h2(t-1))-m_h1(t-1))./tao_m_h-k3*Ph*m_h1(t-1)+k4*m_h2(t-1));
    m_h2(t)=m_h2(t-1)+dt*(k3*Ph*m_h1(t-1)-k4*m_h2(t-1));
  
end
% =========Figure =========
figure
subplot(3,1,1)
plot(time_array/1000,Vt)
%plot(time_array(500000:end)/1000,Vt(500000:end))
ylim([-90 -30])
xlim([10 50])%原代码数据[5 30][10 55]
%%xlabel(['g_L_K=' num2str(g_L_K)])
xlabel('Time(s)')
ylabel('Vt(mV)')

subplot(3,1,2)
plot(time_array/1000,Vr)
%%plot(time_array(500000:end)/1000,Vr(500000:end))
ylim([-90 -30])
%%ylim([-100 0])
xlim([10 50])%原代码数据[5 30]
%%figure
%%plot(time_array,noise1)
xlabel('Time(s)')
ylabel('Vr(mV)')

subplot(3,1,3)
plot(time_array/1000,R_AHP)
%plot(time_array/1000,g)
%ylim([-100 0])
xlim([10 50])%原代码数据[5 30]
%figure
%plot(time_array,noise1)
xlabel('Time(s)')
ylabel('R_AHP')
%

[pks, locs] = findpeaks(Vt(15000:end));%
peakInterval = diff(locs)/10000;
fd=1./peakInterval;%
x=1:length(fd);%

% %% detect spindle number， output is spindle_number.
V1=Vt(15000:500000-1);%%615000-1   500000-1
V=downsample(V1(1:end),100);

load('N2_data_last.mat');
qqq=N2_data_last(3,1:1992000);
jjj=downsample(qqq(1:end),2);
hh=cell(5,1);
hh{1}=jjj;hh{2}=jjj;hh{3}=jjj;hh{4}=jjj;hh{5}=jjj;


MMM=((V-mean(V)).*5-5);

anss=wamsley_spindle_detection(hh,MMM,100); 

[m,n]=find(anss==1);
index=m;
pp=fix(index);
qq=length(pp);

spindle=find(diff(pp)>2);  
spindle_number=length(spindle)+1;
% 
% % ========= Computing tau_ref =========

M1=find(anss==1);
count=0;
for i=1:length(M1)-1
if M1(i+1)~=M1(i)+1
count=count+1;
K(count)=(M1(i+1)-M1(i)+1)*0.01;
end
end
figure
plot(K)
