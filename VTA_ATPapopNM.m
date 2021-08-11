%% VTA model - soma+terminal+ATP+Apoptosis
function VTA_ATPapopNM(dur,gl,mt,filename)

%% CREDITS
% Modified by Group B RISE 2021 based on SNc model from
% Vignayanandam R. Muddapu (Ph.D. scholar)
% C/o Prof. V. Srinivasa Chakravarthy
% Indian Institute of Technology Madras
% India

% Based in part on https://senselab.med.yale.edu/ModelDB/ShowModel?model=266419&file=/modelDB_NicandACh/Both_NicandAChDAmodel.cpp#tabs-2


%% CODE

%%
tic;
% clear;clc;
%dur=1000;
dt=0.1;
% ATP=0.411;
%gl=1;
%mt=1;
% gl=1;
% mt=1;
%filename=strcat('ATP_');

% time1=clock;
% Time parameters & Random seeding
tspan=dt:dt:dur;
% size(tspan)
Ttime=numel(tspan);
% load('randseed.mat')
% srnd = rng;

cai_array=zeros(1,Ttime);atpused_array=zeros(1,Ttime);apop_array=zeros(1,Ttime);
eda_array=zeros(1,Ttime);V_snc_array=zeros(1,Ttime);ros_mit_array=zeros(1,Ttime);
cda_array=zeros(1,Ttime);vda_array=zeros(1,Ttime);phi_er=zeros(1,Ttime);
calb_array=zeros(1,Ttime);cam_array=zeros(1,Ttime);phi_mt=zeros(1,Ttime);
caer_array=zeros(1,Ttime);camt_array=zeros(1,Ttime);LDOPA_array=zeros(1,Ttime);
ATP_array=zeros(1,Ttime);LAC_array=zeros(1,Ttime);PYR_array=zeros(1,Ttime);
GAP_array=zeros(1,Ttime);GSH_array=zeros(1,Ttime);F6P_array=zeros(1,Ttime);
F26P_array=zeros(1,Ttime);PCr_array=zeros(1,Ttime);NADPH_array=zeros(1,Ttime);
ROS_array=zeros(1,Ttime);ASYN_array=zeros(1,Ttime);ASYNA_array=zeros(1,Ttime);
ASYNT_array=zeros(1,Ttime);ASYNG_array=zeros(1,Ttime);LB_array=zeros(1,Ttime);
nai_array=zeros(1,Ttime);ki_array=zeros(1,Ttime);

V_id=zeros(1,Ttime);
V_dp=zeros(1,Ttime);
V_er=zeros(1,Ttime);
V_rel=zeros(1,Ttime);
V_pro=zeros(1,Ttime);

%%%%%%%%%%%%%%%% Initializing the differential equations %%%%%%%%%%%%%%%%%%

% TODO: Modify some of the parameters over the next hundred lines to alter
% the ephys of the cell. You can refer to the model code in
% https://senselab.med.yale.edu/ModelDB/ShowModel?model=266419&file=/modelDB_NicandACh/Both_NicandAChDAmodel.cpp#tabs-2
% for some ideas of parameter values to try. NOTE: this will be complex
% because that code was written by another author who probably uses
% different conventions for variable names. The same biological parameter
% probably has different names here and in that linked code.
%
% check for either membrane capacitance (R, M, rm, or other nomenclature) or leak channel conductance (gL)
% gL = 1/mem_capacitance
% This is a good parameter to alter as part of your parameter search or
% trying to match a VTA neuron
%
% another one would be voltage activation / voltage dependence of Na
% channels (lower priority but also interesting is same parameter for K, Ca
% channels)
% Na parameter might be something like Na_v, Na v shift, Na m ....
%
% And resting membrane potential
%
% Internal concentrations of ions - if not given, solve using the reversal
% potential and the Nernst equation 

V_sncinit = -56; % RMP
Ca_iinit = 0.00006; % internal conc of element mM
Na_iinit = 11.5931; % internal conc of element 
K_iinit = 125.1474; % internal conc of element 
Calbinit = 0.0026; % calcium buffering or binding parameter
Caminit = 0.0222; %
m_calinit = 0.006271; % voltage activation inactivation of this channel type
m_nainit = 0.0952; % voltage activation inactivation of this channel type
h_nainit = 0.1848; % voltage activation inactivation of this channel type
O_hcninit = 0.003; % voltage activation inactivation of this channel type
m_kdrinit = 0.0932; % voltage activation inactivation of this channel type
y_pcinit = 0.483;
y_nkinit = 0.6213;
Ca_erinit = 1.0*0.001; %mM
Ca_mtinit = 0.1*0.001; % mM % 0.4e-3
cdainit = 1e-4; %mM%1e-4
vdainit = 500; %mM 500
edainit = 4e-6; %mM
Iextinit=0; %Iext
ATPusedinit=0;
calinit=1; %mM
cai_calinit=0; %mM
cal_actinit=0; %mM
casp12init=1; %mM
cal_act_casp12init=0; %mM
casp12_actinit=0; %mM
casp9init=1; %mM
casp12_act_casp9init=0; %mM
casp9_actinit=0; %mM
casp3init=1; %mM
casp9_act_casp3init=0; %mM
casp3_actinit=0; %mM
apopinit=0; %mM
ROS_mitinit=0;
PTP_mit_actinit=0;
Cytc_mitinit=1;
Cytcinit=0;
Cytc_casp9init=0;
IAPinit=0;
casp9_act_IAPinit=0;
casp3_act_IAPinit=0;
NADPHinit=250*0.001;%mM
GSHinit=2500*0.001;%mM
F6Pinit = 0.175883476634895;%0.2
F26Pinit = 0.002191750879602;%0.001
GAPinit = 0.082507126186107;%0.0405
PYRinit = 0.123910489378719;%0.1
LACinit = 0.598605032933119;%0.5
ATPinit = 2.395615876085214;%2.402
PCrinit = 18.044071098085976;%18.14
ROSinit=1*0.001;%mM
ASYNinit=100*0.001;%mM
ASYNAinit=1*0.001;%mM
ASYNTinit=0.01*0.001;%mM
ASYNGinit=0*0.001;%mM
LBinit=0*0.001;%mM
LDOPAinit=36e-5;%mM

V_snc=V_sncinit;m_cal=m_calinit;m_na=m_nainit;
h_na=h_nainit;O_hcn=O_hcninit;Calb=Calbinit;
Cam=Caminit;y_nk=y_nkinit;y_pc=y_pcinit;m_kdr=m_kdrinit;
K_i=K_iinit;Na_i=Na_iinit;Ca_i=Ca_iinit;Ca_er=Ca_erinit;Ca_mt=Ca_mtinit;

cda=cdainit;vda=vdainit;
eda=edainit;Iexts=Iextinit;ATPused=ATPusedinit;
cal=calinit;cai_cal=cai_calinit;cal_act=cal_actinit;casp12=casp12init;
cal_act_casp12=cal_act_casp12init;casp12_act=casp12_actinit;casp9=casp9init;
casp12_act_casp9=casp12_act_casp9init;casp9_act=casp9_actinit;casp3=casp3init;
casp9_act_casp3=casp9_act_casp3init;casp3_act=casp3_actinit;apop=apopinit;

ROS_mit=ROS_mitinit;
PTP_mit_act=PTP_mit_actinit;
Cytc_mit=Cytc_mitinit;
Cytc=Cytcinit;
Cytc_casp9=Cytc_casp9init;
IAP=IAPinit;
casp9_act_IAP=casp9_act_IAPinit;
casp3_act_IAP=casp3_act_IAPinit;
NADPH=NADPHinit;
GSH=GSHinit;
F6P=F6Pinit;
F26P=F26Pinit;
GAP=GAPinit;
PYR=PYRinit;
LAC=LACinit;
ATP=ATPinit;
PCr=PCrinit;
ROS=ROSinit;
ASYN=ASYNinit;
ASYNA=ASYNAinit;
ASYNT=ASYNTinit;
ASYNG=ASYNGinit;
LB=LBinit;
LDOPA=LDOPAinit;

sim_mM=1e-3;
sim_mM_msec=1e-3/3.6e6;
sim_msec=1/3.6e6;
sim_msec_mM=1/((1e-3)*(3.6e6));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants
%Soma
R = 8314.472; %mJ/mol. K
T = 310.15; %K
F = 96485.30929; %coulomb/mol.
Ca_o = 1.8; %mM
Na_o = 137; %mM
K_o = 5.4; %mM
vol_pmu = 5; %pl
fr_cyt = 0.5;
C_sp = 1e6; %pF/cm2
SVR_pmu = 1.6667e4; %1/cm
% ATP = 0.411; %mM 0.411 - bursting
Calbtot = 0.005; % mM
Camtot = 0.0235; % mM
kcal_1 = 10; %1/mM.ms
kcal_2 = 2e-3; %1/ms
kcam_cd = 0.003; %1/ms
kcam_nd = 3; %1/ms
g_cal = 2101.2; %pA/mM
g_na = 907.68; %pA/mM
A_mna = 1.9651; %1/ms
B_mna = 0.0424; %1/ms
A_hna = 9.566e-5; %1/ms
B_hna = 0.5296; %1/ms
za_mna = 1.7127;
zb_mna = 1.5581;
za_hna = -2.4317;
zb_hna = -1.1868;
g_nalk = 0.0053; %pA/mM
g_nahcn = 51.1; %pA/mM
cAMP = 1e-5; %mM
g_ksk = 2.2515; %pA/mM
g_kdr = 31.237; %nS
g_kir = 13.816; %nS
k_2pc = 0.001; %1/ms
k_3pc = 0.001; %1/ms
k_4pc = 1; %1/ms
K_pco = 2; %mM
k_pmca = 2.233;
dell = 0.35;
k_xm = 0.0166; %pA
k_2nk = 0.04; %1/ms
k_3nk = 0.01; %1/ms
k_4nk = 0.165; %1/ms
K_nknai = 4.05; %mM
K_nknao = 69.8; %mM
K_nkki = 32.88; %mM
K_nkko = 0.258; %mM
k_nk = 1085.7; %pA
V_tau = (R*T)/F;
vol_cyt = fr_cyt*vol_pmu;
P_c = 1.00000/(1.00000+cAMP/0.00116300);
P_o = 1.00000/(1.00000+cAMP/1.45000e-05);
P_E2Spc = 1.00000/(1.00000+K_pco/Ca_o);
A_pmu = (SVR_pmu*vol_pmu*0.00100000*0.00100000*0.00100000)/1.00000;
P_E2pc = 1.00000-P_E2Spc;
beta_pc = k_2pc*P_E2Spc+k_4pc*P_E2pc;

%CICR
% ER
rho_er = 0.01; % rho_er
beta_er = 0.0025; %beta_er
k_pump = 20.0/1000; % 1/ms
k_ch = 3000.0/1000; %1/ms
K1 = 5.0*0.001; %mM
k_leak = 0.05/1000; %1/ms

% Mito
rho_mt = 0.01; % rho_mt
beta_mt = 0.0025; % beta_mt
k_in = 0.0055*300*0.001/1000; % mM/ms
K2 = 0.8*0.001; % mM
k_out = 125.0/1000; % 1/ms
k_m = 0.00625/1000; % 1/ms
K3 = 5.0*0.001; % mM

% DA terminal (Tello-Bravo (2012))
krel = 0.031; %(mM) 0.055 % latest-0.063
psi = 17.4391793; %(mM/ms)
nRRP = 5; % :RANGE ~ 10-30
Veda_max = 1e-6; %(mM/ms)
Keda = 3e-5; %(mM)
kcomt = 0.0083511; %(1/ms)
%vda = 500; %(mM)
vdao = 500; %(mM)
vdas = 1e-2; %(mM)
dara = 5e-5; %(mM)
dars = 1e-2; %(mM)
Vsynt_max = 250e-5;%(mM/ms)%30.2e-6 %25e-6 latest-50e-5
Ksynt = 35e-4; %(mM)
Ktyr = 46e-3; %(mM)
TYR = 126e-3; %(mM)
Kicda = 11e-2; %(mM)
Kieda = 46e-3; %(mM)
Vcda_max = 0.2*133.33e-6; %(mM/ms)%133.33e-6*0.03 latest-0.009*133.33e-6
Kcda = 238e-4; %(mM)
kmao = 0.00016; %(1/ms)

% Apoptosis pathway (Hong et.al., (2012))
k3f=1; % (muM*sec)-1
k3b=1/1e3; % (sec)-1
k4f=1/1e3; % (sec)-1
k5f=1; % (muM*sec)-1
k5b=1/1e3; % (sec)-1
k6f=1/1e3; % (sec)-1
k7f=10; % (muM*sec)-1
k7b=0.5/1e3; % (sec)-1
k8f=1/1e3; % (sec)-1
k9f=10; % (muM*sec)-1
k9b=0.5/1e3; % (sec)-1
k10f=0.1/1e3; % (sec)-1
k11f=1; % (muM*sec)-1

k29f=0.5; % (mM*msec)-1
k30f=0.5; % (mM*msec)-1
k31f=1; % (mM*msec)-1
k27f=1; % (mM*msec)-1
k27b=1/1e3; % (msec)-1
k28f=1/1e3; % (msec)-1
k12f=5; % (mM*msec)-1
k12b=0.0035/1e3; % (msec)-1
k13f=5; % (mM*msec)-1
k13b=0.0035/1e3; % (msec)-1

Mit=1;
Sig_ers=0;%0.0001;
Sig_mts=0;%0.0001;
PTP_mit=1;

% Energy Metabolism
GLCe=1;%mM
Vmax_hk = 2.5/1000;%mM/ms
Km_ATP_hk = 0.5;%mM
KI_F6P = 0.068;%mM
Vmax_pfk = 3.85/1000;%mM/ms
Km_ATP_pfk = 0.05;%mM
Km_F6P_pfk = 0.18;%mM
Km_F26P_pfk = 0.01;%mM
Vmaxf_pfk2 = 2e-04/1000;%mM/ms
Vmaxr_pfk2 = 1.036e-04/1000;%mM/ms
Km_ATP_pfk2 = 0.05;%mM
Km_F6P_pfk2 = 0.01;%mM
Km_F26P_pfk2 = 0.0001;%mM
Vmax_pk = 5.0/1000;%mM/ms
Km_ADP_pk = 0.005;%mM
Km_GAP_pk = 0.4;%mM
Vmax_op = 1.0/1000;%mM/ms
Km_ADP_op = 0.005;%mM
Km_PYR_op = 0.5;%mM
kf_ldh = 12.5/1000;%1/ms
kr_ldh = 2.5355/1000;%1/ms
kf_ck = 3.0/1000;%1/mM.ms
kr_ck = 1.26/1000;%1/mM.ms
PCr_tot = 20.0;%mM
Vmax_ATPase = 0.9355/1000;%mM/ms
Km_ATP = 0.5;%mM
Vlac_0 = 0.355/1000;%mM/ms
K_lac_eff = 0.71/1000;%1/ms
K_lac = 0.641;
ANP = 2.51;%mM
Q_adk = 0.92;
nATP = 0.4;
KI_ATP = 1.0;%mM
nAMP = 0.5;
Ka_AMP = 0.05;%mM
Kamp_pfk2 = 0.005;%mM
nh_amp = 2;
beta_ldh_ros=0.25;
Kldh_ros=10*sim_mM;%muM

snc_firings=[];

kf_gr=0.65*sim_msec_mM;%1/mM.ms
kr_gr=1.25e-3*sim_msec_mM;%1/mM.ms
GSH_tot=2500*sim_mM;%mM
NADPH_tot=250*sim_mM;%mM
Vmax_ppp=1.43e6*sim_mM_msec;%mM/ms
Ki_nadph=20;

% PD pathology pathways
eta_op_max=0.995;
beta_op_asyn=0.08;
Kasyn=8.5*sim_mM; %mM
Kros_cat=235*1*sim_msec; %1/ms
Vros_ex=0;%0.1
Kros_dopa=1500*sim_msec_mM;%1/mM.ms
Kros_dox=0.27*sim_msec;%1/ms
Kasyn_syn=50*sim_mM_msec;%mM/ms
Kasyn_ox=7e-5*sim_msec_mM;%1/mM.ms
Kasyn_to=0.5*sim_msec;%1/ms
Krasyn_agg=7.5e-4*sim_msec;%1/ms
Kasyn_agg=7.5*sim_mM;%mM
Ub_tot=10.5*sim_mM;%mM
Kasyn_tag=2.75e-7*sim_msec_mM;%1/mM.ms
Krasyn_prt=7.5e-4*sim_msec;%1/ms
Kasyn_prt=5*sim_mM;%mM
beta_asyn_prt=0.25;
Kasyn_lyso=7.5e-5*sim_msec;%1/ms
Krasyn_lb=7.5e-5*sim_msec;%1/ms
Kasyn_lb=5*sim_mM;%

% LDOPA uptake
%Reed et al(2012)
Vaadc_max = 2.78e-6;% (mM/ms)
Kaadc = 0.13;% (mM)
Vtran_max = 5.11e-7;% (mM/ms)5.94e-8
% sLD = 36e-3;% (mM)
% sLD = 0;% (mM)
sTYR = 63e-3;% (mM)
sTRP = 82e-3;% (mM)
Ksld = 32e-3;% (mM)
Kstyr = 64e-3;% (mM)
Kstrp = 15e-3;% (mM)
sLD=3.63685e-3;%35.39059803e-3; %LDOPA=36e-5 Org-36e-3 %mM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Stimulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Istim=0.000; % Phasic bursting Istim=0.0001
del=5000/dt;
dur=200/dt;
sigg1=0;sigg2=0;phier=0;phimt=0;ada=1;
counttt=1;lam=0.0001;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:Ttime
    
    if(k > 5000/dt)
%         enfatra=0.2;%0.2
%         enfamit=0.1;
        
        enfatra=gl;%0.2
        enfamit=mt;
        
%         Renfatra=1.*exp(-counttt.*lam);
%         Renfamit=1.*exp(-counttt.*lam);
%         counttt=counttt+1;
%         
%         enfatra=Renfatra;
%         enfamit=Renfamit;
%         
%         enfatra(enfatra<0.2)=0.2;
%         enfamit(enfamit<0.1)=0.1;
        
    else
        enfatra=1;
        enfamit=1;
    end
    
    if Ca_mt>0.019
        sigg1=1;
    end
    
    if phier>0.0
        sigg2=1;
    end
    
    if sigg1==1
        Sig_mts=0.01;
    end
    
    if sigg2==1
        Sig_ers=0.01;
    end

    k_1pc = 1.00000/(1.00000+0.100000/ATP);
    k_1nk = 0.370000/(1.00000+0.0940000/ATP);
    
    %Soma
    %Membrane potential
    VD = V_snc/V_tau;
    
    % HCN current
    kf_free = 0.00600000/(1.00000+exp((V_snc+87.7000)/6.45000));
    kf_bnd = 0.0268000/(1.00000+exp((V_snc+94.2000)/13.3000));
    kf_hcn = kf_free*P_c+kf_bnd*(1.00000-P_c);
    kr_free = 0.0800000/(1.00000+exp(-(V_snc+51.7000)/7.00000));
    kr_bnd = 0.0800000/(1.00000+exp(-(V_snc+35.5000)/7.00000));
    kr_hcn = kr_free*P_o+kr_bnd*(1.00000-P_o);
    
    % Calcium binding proteins
    CaCalb = Calbtot-Calb;
    J_calb = kcal_1*Calb*Ca_i-kcal_2*CaCalb;
    
    CaCam = Camtot-Cam;
    kcam_cb = 12000.0*(power(Ca_i, 2.00000));
    kcam_nb = 3.70000e+06*(power(Ca_i, 2.00000));
    alpha_cam = kcam_cb*kcam_nb*(1.00000/(kcam_cb+kcam_nd)+1.00000/(kcam_cd+kcam_nd));
    beta_cam = kcam_cd*kcam_nd*(1.00000/(kcam_cb+kcam_nd)+1.00000/(kcam_cd+kcam_nd));
    J_cam = alpha_cam*Cam-beta_cam*CaCam;
    
    
    K_pci = (173.600/(1.00000+CaCam/5.00000e-05)+6.40000)*1.00000e-05;
    P_E1Spc = 1.00000/(1.00000+K_pci/Ca_i);
    P_E1pc = 1.00000-P_E1Spc;
    alpha_pc = k_1pc*P_E1Spc+k_3pc*P_E1pc;
    
    V_Ca = 0.500000*log(Ca_o/Ca_i);
    h_cal = 0.000450000/(0.000450000+Ca_i);
    I_CaL = (g_cal*m_cal*h_cal*(power(Ca_i*Ca_o, 1.0/2))*sinh(VD-V_Ca))/(sinh(VD)/VD);
    K_pmca = k_pmca*((10.5600*CaCam)/(CaCam+5.00000e-05)+1.20000);
    I_pmca = K_pmca*(k_1pc*P_E1Spc*y_pc-k_2pc*P_E2Spc*(1.00000-y_pc))*1.00000;
    Dr = (1.00000+0.00100000*((power(Na_i, 3.00000))*Ca_o+(power(Na_o, 3.00000))*Ca_i))*(1.00000+Ca_i/0.00690000);
    I_xm = (k_xm*((power(Na_i, 3.00000))*Ca_o*exp(dell*VD)-(power(Na_o, 3.00000))*Ca_i*exp((dell-1.00000)*VD)))/Dr;
    J_ca = (-1.00000/(2.00000*F*vol_cyt))*((I_CaL+2.00000*I_pmca)-2.00000*I_xm);
    
    V_Na = log(Na_o/Na_i);
    O_na = (power(m_na, 3.00000))*h_na;
    I_Na = (g_na*O_na*(power(Na_i*Na_o, 1.0/2))*sinh(0.500000*(VD-V_Na)))/(sinh(0.500000*VD)/(0.500000*VD));
    I_Nalk = (g_nalk*(power(Na_i*Na_o, 1.0/2))*sinh(0.500000*(VD-V_Na)))/(sinh(0.500000*VD)/(0.500000*VD));
    I_NaHCN = (g_nahcn*O_hcn*(power(Na_i*Na_o, 1.0/2))*sinh(0.500000*(VD-V_Na)))/(sinh(0.500000*VD)/(0.500000*VD));
    P_E1Snk = 1.00000/(1.00000+(K_nknai/Na_i)*(1.00000+K_i/K_nkki));
    Na_eff = Na_o*exp(-0.820000*VD);
    P_E2Snk = 1.00000/(1.00000+(K_nknao/Na_eff)*(1.00000+K_o/K_nkko));
    I_nk = k_nk*(k_1nk*P_E1Snk*y_nk-k_2nk*P_E2Snk*(1.00000-y_nk))*1.00000;
    J_Na = (-1.00000/(F*vol_cyt))*(3.00000*I_nk+3.00000*I_xm+I_Na+I_Nalk+I_NaHCN);
    
    P_E1Dnk = 1.00000/(1.00000+(K_nkki/K_i)*(1.00000+Na_i/K_nknai));
    alpha_nk = k_1nk*P_E1Snk+k_3nk*P_E1Dnk;
    P_E2Dnk = 1.00000/(1.00000+(K_nkko/K_o)*(1.00000+Na_eff/K_nknao));
    beta_nk = k_2nk*P_E2Snk+k_4nk*P_E2Dnk;
    
    V_K = log(K_o/K_i);
    O_sk = (power(Ca_i, 4.20000))/(power(0.000350000, 4.20000)+power(Ca_i, 4.20000));
    I_Ksk = (g_ksk*O_sk*(power(K_i*K_o, 1.0/2))*sinh(0.500000*(VD-V_K)))/(sinh(0.500000*VD)/(0.500000*VD));
    O_kdr = power(m_kdr, 3.00000);
    I_Kdr = g_kdr*O_kdr*(V_snc-V_K*V_tau);
    O_kir = 1.00000/(1.00000+exp((V_snc+85.0000)/12.1000));
    I_Kir = g_kir*O_kir*(V_snc-V_K*V_tau);
    I_K = I_Ksk+I_Kdr+I_Kir;
    J_K = (-1.00000/(F*vol_cyt))*(I_K-2.00000*I_nk);
    
    % ER
    J_pump = k_pump*Ca_i*ATP; %J_pump
    J_ch = k_ch*((power(Ca_i, 2.00000))/(power(K1, 2.00000)+power(Ca_i, 2.00000)))*(Ca_er-Ca_i); %J_ch
    J_leak = k_leak*(Ca_er-Ca_i); %J_leak
    
    % Mito
    J_out = (k_out*((power(Ca_i, 2.00000))/(power(K3, 2.00000)+power(Ca_i, 2.00000)))+k_m)*Ca_mt; % J_out
    J_in = k_in*((power(Ca_i, 8.00000))/(power(K2, 8.00000)+power(Ca_i, 8.00000)));%*ATP; % J_in
    
    % Calcium dynamics
    J_Ca = J_ca-1*(J_calb+4.00000*J_cam)-J_pump+J_ch+J_leak-J_in+J_out;
    
    adca=0;
    %Terminal
    Vsynt = Vsynt_max/(((Ksynt/(adca+Ca_i))^4)+1);
    jsynt = (Vsynt/(1+((Ktyr/TYR)*(1+(cda/Kicda)+(eda/Kieda)))));
    
    jvmat = ada*Vcda_max * MM_kin(cda,Kcda,1);%*(ATP);
    
    jida = kmao * cda;
    
    %     nRRP = (40/((1+exp(-(vda-vdao)/vdas))*(1+exp((eda-dara)/dars))));
    %     nRRP = 1;
    
    % ATP-dependent DA packing
%     ada=RescaleRange(ATP,0.2,2.3,0.001,1);
    ada=0.001*(exp(3*ATP));
    
    % ATP-dependent vescile recycling
%     nRRP=5;
    nRRP=1*(exp(0.7*ATP));
    %     if k>3000/dt && k<5000/dt
    % %         nRRP=1;
    %         ada=0.2;
    %     else
    %         ada=1;
    % %         nRRP=10;
    %     end
    
    
    
    prob = 0.14 * MM_kin((adca+Ca_i),krel,4);
    jrel = psi * nRRP * prob;
    
%     Vrel(k)=jrel;
%     Vmat(k)=jvmat;
    
    jdat = Veda_max * MM_kin(eda,Keda,1);
    
    jeda = kcomt * eda;
    
    jldopa = Vaadc_max * MM_kin(LDOPA,Kaadc,1);
    
    % Energy metabolism
    % Energy consumed by active pumps
    % V_pumps=0;
    V_pumps1 = 1*(1.00000/(F*vol_cyt))*(I_nk+I_pmca);
    V_pumps2 = 1*(jvmat);
    V_pumps3 = 100*jrel;
    v_stim=0;
    % v_stim1=(1.00000/(F*vol_cyt))*(I_nk+I_pmca);
    v_stim1=0.0*(I_nk+I_pmca);
    J_er=(beta_er/rho_er)*(J_pump);
    
    V_id(k)=V_pumps1;
    V_dp(k)=V_pumps2;
    V_er(k)=J_er;
    V_rel(k)=V_pumps3;
    
    V_pumps=V_id(k)+V_dp(k)+V_er(k)+V_rel(k);
    
    uADP = power(Q_adk, 2.00000)+4.00000*Q_adk*(ANP/ATP-1.00000);
    ADP = (ATP/2.00000)*(-Q_adk+power(uADP, 1.0/2));
    Cr = PCr_tot-PCr;
    V_ck = 0*(kf_ck*PCr*ADP-kr_ck*Cr*ATP);
    
    ATP_inh = power((1.00000+nATP*(ATP/KI_ATP))/(1.00000+ATP/KI_ATP), 4.00000);
    V_pk = Vmax_pk*(GAP/(GAP+Km_GAP_pk))*(ADP/(ADP+Km_ADP_pk))*ATP_inh;
    pa=1;
    V_op = enfamit*Vmax_op*((pa*PYR)/((pa*PYR)+Km_PYR_op))*(ADP/(ADP+Km_ADP_op))*(1.00000/(1.00000+0.100000*(ATP/ADP)));
    
    AMP = ANP-(ATP+ADP);
    AMP_act = power((1.00000+AMP/Ka_AMP)/(1.00000+nAMP*(AMP/Ka_AMP)), 4.00000);
    V_pfk = Vmax_pfk*(F6P/(F6P+Km_F6P_pfk))*(ATP/(ATP+Km_ATP_pfk))*(F26P/(F26P+Km_F26P_pfk))*ATP_inh*AMP_act;
    
    eta_ldh=1-beta_ldh_ros*((ROS^4)/((ROS^4)+(Kldh_ros^4)));
    V_ldh = 1*eta_ldh*(kf_ldh*PYR-kr_ldh*LAC);
    V_lac = Vlac_0*(1.00000+v_stim1*K_lac)-K_lac_eff*LAC;
    
    V_hk = enfatra*Vmax_hk*(ATP/(ATP+Km_ATP_hk))*(power(1.00000+power(F6P/KI_F6P, 4.00000), -1.00000))*GLCe;
    AMP_pfk2 = (power(AMP/Kamp_pfk2, nh_amp))/(1.00000+power(AMP/Kamp_pfk2, nh_amp));
    V_pfk2 = Vmaxf_pfk2*(ATP/(ATP+Km_ATP_pfk2))*(F6P/(F6P+Km_F6P_pfk2))*AMP_pfk2-Vmaxr_pfk2*(F26P/(F26P+Km_F26P_pfk2));
    
    V_ATPase = Vmax_ATPase*(ATP/(ATP+Km_ATP))*(1.00000+v_stim);
    dAMP_dATP = -1.00000+Q_adk/2.00000+-(0.500000*(power(uADP, 1.0/2)))+Q_adk*(ANP/(ATP*(power(uADP, 1.0/2))));
    
    GSSG=(GSH_tot-GSH)/2;
    NADP=NADPH_tot-NADPH;
    Vppp = Vmax_ppp*(F6P/(F6P+Km_F6P_pfk))/(1+((NADPH/NADP)/Ki_nadph));
    Vgr = kf_gr*GSSG*NADPH-kr_gr*GSH*NADP;
    
    % PD pathology pathways
    eta_op=eta_op_max-beta_op_asyn*(((ASYNA^4)/((ASYNA^4)+(Kasyn^4))));
    
    Vros_leak=(0.5282/ATP)*(1-eta_op)*V_op;%0.221 or (0.5282/ATP)
    %     Vros_leak=0.221*(1-eta_op)*V_op;%0.221 or (0.5282/ATP)
    
    Vros_cat=Kros_cat*ROS;
    
    Vros_dopa=0*Kros_dopa*(((ASYNA^4)/((ASYNA^4)+(Kasyn^4))));
    
    Vros_dox=Kros_dox*GSH*ROS;
    
    Vasyn_syn=Kasyn_syn;
    
    Vasyn_ox=Kasyn_ox*ROS*ASYN;
    
    Vasyn_to=Kasyn_to*ASYN;
    
    Vasyn_agg=Krasyn_agg*ASYNA*(((ASYNA^6)/((ASYNA^6)+(Kasyn_agg^6))));
    
    Ub=Ub_tot-ASYNT;
    Vasyn_tag=Kasyn_tag*ASYNA*Ub*ATP;
    
    Vasyn_prt=Krasyn_prt*ASYNT*ATP*(1-beta_asyn_prt*(((ASYNG^4)/((ASYNG^4)+(Kasyn_prt^4)))));
    
    Vasyn_lyso=Kasyn_lyso*ASYNG*ATP;
    
    Vasyn_lb=Krasyn_lb*ASYNG*(((ASYNG^6)/((ASYNG^6)+(Kasyn_lb^6))));
    
    V_pro(k)=25*Vasyn_prt+1*(3*Vasyn_tag+10*Vasyn_lyso);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Stimulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Ibg=0;
    if(k < del + dur && k > del)
        Iapp = Ibg +Istim;
    else
        Iapp = Ibg;
    end
    Iext=Iapp;
    
    if k>3000/dt && k<7500/dt
%         sLD=36e-3;
        sLD=3.63685e-3;
    else
        sLD=3.63685e-3;
    end
    
    
    
    %     if ATP<0.2
    %         ATP=0.2;
    %     end
    
    %%%%%%%%%%%%%%%%%%%%%% Differential equations %%%%%%%%%%%%%%%%%%%%%%%%%
    
    V_sncnxt = V_snc + (((F*vol_cyt)/(C_sp*A_pmu))*(J_Na+J_K+2.00000*J_Ca+Iext))*dt;
    Ca_inxt = Ca_i + (J_Ca)*dt;
    Na_inxt = Na_i + (J_Na)*dt;
    K_inxt = K_i + (J_K)*dt;
    Calbnxt = Calb + (-J_calb)*dt;
    Camnxt = Cam + (-J_cam)*dt;
    m_calnxt = m_cal + ((1.00000/(1.00000+exp(-(V_snc+15.0000)/7.00000))-m_cal)/(7.68000*exp(-(power((V_snc+65.0000)/17.3300, 2.00000)))+0.723100))*dt;
    m_nanxt = m_na + (A_mna*exp(za_mna*VD)*(1.00000-m_na)-B_mna*exp(-zb_mna*VD)*m_na)*dt;
    h_nanxt = h_na + (A_hna*exp(za_hna*VD)*(1.00000-h_na)-B_hna*exp(-zb_hna*VD)*h_na)*dt;
    O_hcnnxt = O_hcn + (kf_hcn*(1.00000-O_hcn)-kr_hcn*O_hcn)*dt;
    m_kdrnxt = m_kdr + ((1.00000/(1.00000+exp(-(V_snc+25.0000)/12.0000))-m_kdr)/(18.0000/(1.00000+exp((V_snc+39.0000)/8.00000))+1.00000))*dt;
    y_pcnxt = y_pc + (beta_pc*(1.00000-y_pc)-alpha_pc*y_pc)*dt;
    y_nknxt = y_nk + (beta_nk*(1.00000-y_nk)-alpha_nk*y_nk)*dt;
    ATPusednxt = ATPused+(-ATPused+(1.00000/(F*vol_cyt))*(I_nk+I_pmca))*dt;%ATPused
    Ca_ernxt = Ca_er + ((beta_er/rho_er)*(J_pump-(J_ch+J_leak)))*dt;
    Ca_mtnxt = Ca_mt + ((beta_mt/rho_mt)*(J_in-J_out))*dt;
    cdanxt = cda+(jsynt + jdat - jvmat - jida + jldopa)*dt;%cda
    vdanxt = vda+(jvmat - jrel)*dt;%vda
    edanxt = eda+(jrel - jdat - jeda)*dt;%eda
    calnxt = cal+(-k3f*(Sig_ers*cal)+k3b*(cai_cal))*dt;%cal
    cai_calnxt = cai_cal+(k3f*(Sig_ers*cal)-k3b*(cai_cal)-k4f*(cai_cal))*dt;%cai_cal
    cal_actnxt = cal_act+(k4f*(cai_cal)-k5f*(cal_act*casp12)+k5b*(cal_act_casp12))*dt;%cal_act
    casp12nxt = casp12+(-k5f*(cal_act*casp12)+k5b*(cal_act_casp12))*dt;%casp12
    cal_act_casp12nxt = cal_act_casp12+(k5f*(cal_act*casp12)-k5b*(cal_act_casp12)-k6f*(cal_act_casp12))*dt;%cal_act_casp12
    casp12_actnxt = casp12_act+(k6f*(cal_act_casp12)-k7f*(casp12_act*casp9)+k7b*(casp12_act_casp9))*dt;%casp12_act
    casp9nxt = casp9+(-k7f*(casp12_act*casp9)+k7b*(casp12_act_casp9))*dt;%casp9
    casp12_act_casp9nxt = casp12_act_casp9+(k7f*(casp12_act*casp9)-k7b*(casp12_act_casp9)-k8f*(casp12_act_casp9))*dt;%casp12_act_casp9
    casp9_actnxt = casp9_act+(k8f*(1*casp12_act_casp9)+k9b*(casp9_act_casp3)-k9f*(casp9_act*casp3)+1*k28f*Cytc_casp9-k12f*casp9_act*IAP+k12b*casp9_act_IAP)*dt;%casp9_act
    casp3nxt = casp3+(-k9f*(casp9_act*casp3)+k9b*(casp9_act_casp3))*dt;%casp3
    casp9_act_casp3nxt = casp9_act_casp3+(-k10f*(casp9_act_casp3)-k9b*(casp9_act_casp3)+k9f*(casp9_act*casp3))*dt;%casp9_act_casp3
    casp3_actnxt = casp3_act+(k10f*(casp9_act_casp3)-k11f*(casp9_act*casp3_act)-k13f*casp3_act*IAP+k13b*casp3_act_IAP)*dt;%casp3_act
    apopnxt =  apop+(k11f*(casp9_act*casp3_act))*dt;%apop
    
    F6Pnxt = F6P+(V_hk-(V_pfk-V_pfk2)-Vppp*(1/6))*dt;%F6P
    F26Pnxt = F26P+(V_pfk2)*dt;%F26P
    GAPnxt = GAP+(2.00000*V_pfk-V_pk)*dt;%GAP
    PYRnxt = PYR+(V_pk-(V_op+V_ldh))*dt;%PYR
    LACnxt = LAC+(2.25000*V_ldh+V_lac)*dt;%LAC
    ATPnxt = ATP+(((1*(1*2.00000*V_pk+15.0000*eta_op*V_op+V_ck))-(V_hk+V_pfk+V_pfk2+V_ATPase+V_pumps+25*Vasyn_prt+1*(3*Vasyn_tag+10*Vasyn_lyso)))*(power(1.00000-dAMP_dATP, -1.00000)))*dt;%ATP
    PCrnxt = PCr+(-V_ck)*dt;%PCr
    ROSnxt = ROS+(Vros_leak + Vros_ex - Vros_cat + Vros_dopa - Vros_dox)*dt;%ROS
    ASYNnxt = ASYN+(Vasyn_syn - Vasyn_ox - Vasyn_to)*dt;%ASYN
    ASYNAnxt = ASYNA+(Vasyn_ox - Vasyn_agg - Vasyn_tag)*dt;%ASYNA
    ASYNTnxt = ASYNT+(Vasyn_tag - Vasyn_prt)*dt;%ASYNT
    ASYNGnxt = ASYNG+(Vasyn_agg - Vasyn_lyso - Vasyn_lb)*dt;%ASYNG
    LBnxt = LB+(Vasyn_lb)*dt;%LB
    
    ROS_mitnxt = ROS_mit+(k29f*Sig_mts*Mit)*dt;%ROS_mit
    PTP_mit_actnxt = PTP_mit_act+(k30f*ROS_mit*PTP_mit)*dt;%PTP_mit_act
    Cytc_mitnxt = Cytc_mit+(-k31f*PTP_mit_act*Cytc_mit)*dt;%Cytc_mit
    Cytcnxt = Cytc+(-k27f*Cytc*casp9+k27b*Cytc_casp9+k31f*PTP_mit_act*Cytc_mit)*dt;%Cytc
    Cytc_casp9nxt = Cytc_casp9+(k27f*Cytc*casp9-k27b*Cytc_casp9-k28f*Cytc_casp9)*dt;%Cytc_casp9
    IAPnxt = IAP+(-k12f*casp9_act*IAP+k12b*casp9_act_IAP-k13f*casp3_act*IAP+k13b*casp3_act_IAP)*dt;%IAP
    casp9_act_IAPnxt = casp9_act_IAP+(k12f*casp9_act*IAP-k12b*casp9_act_IAP)*dt;%casp9_act_IAP
    casp3_act_IAPnxt = casp3_act_IAP+(k13f*casp3_act*IAP-k13b*casp3_act_IAP)*dt;%casp3_act_IAP
    
    NADPHnxt = NADPH+(2*Vppp - Vgr)*dt;%NADPH
    GSHnxt = GSH+(2*Vgr - 2*Vros_dox)*dt;%
    
    LDOPAnxt = LDOPA+(((Vtran_max*sLD)/(Ksld*(1+(sTYR/Kstyr)+(sTRP/Kstrp))+sLD)) - jldopa)*dt;%ldopa
    
    
    V_snc=V_sncnxt;m_cal=m_calnxt;m_kdr=m_kdrnxt;m_na=m_nanxt;
    h_na=h_nanxt;O_hcn=O_hcnnxt;Calb=Calbnxt;
    Cam=Camnxt;y_nk=y_nknxt;y_pc=y_pcnxt;
    K_i=K_inxt;Na_i=Na_inxt;Ca_i=Ca_inxt;
    Ca_er=Ca_ernxt;Ca_mt=Ca_mtnxt;
    cda=cdanxt;vda=vdanxt;eda=edanxt;ATPused=ATPusednxt;
    cal=calnxt;cai_cal=cai_calnxt;cal_act=cal_actnxt;
    casp12=casp12nxt;cal_act_casp12=cal_act_casp12nxt;casp12_act=casp12_actnxt;
    casp9=casp9nxt;casp12_act_casp9=casp12_act_casp9nxt;casp9_act=casp9_actnxt;
    casp3=casp3nxt;casp9_act_casp3=casp9_act_casp3nxt;casp3_act=casp3_actnxt;
    apop=apopnxt;
    ROS_mit=ROS_mitnxt;
    PTP_mit_act=PTP_mit_actnxt;
    Cytc_mit=Cytc_mitnxt;
    Cytc=Cytcnxt;
    Cytc_casp9=Cytc_casp9nxt;
    IAP=IAPnxt;
    casp9_act_IAP=casp9_act_IAPnxt;
    casp3_act_IAP=casp3_act_IAPnxt;
    NADPH=NADPHnxt;
    GSH=GSHnxt;
    F6P=F6Pnxt;
    F26P=F26Pnxt;
    GAP=GAPnxt;
    PYR=PYRnxt;
    LAC=LACnxt;
    ATP=ATPnxt;
    PCr=PCrnxt;
    ROS=ROSnxt;
    ASYN=ASYNnxt;
    ASYNA=ASYNAnxt;
    ASYNT=ASYNTnxt;
    ASYNG=ASYNGnxt;
    LB=LBnxt;
    LDOPA=LDOPAnxt;
    
    V_snc_array(k)=V_snc;
    
    %     inds=find(V_snc_array(k) >= -20 && V_snc_array(k) > V_snc_array(k-1) && V_snc_array(k) > V_snc_array(k+1));
    inds=find(V_snc <=80 & V_snc >=-20);
    snc_firings=[snc_firings; k+0*inds,inds+0*inds];
    
    nai_array(k)=Na_i;ki_array(k)=K_i;
    
    cai_array(k)=Ca_i;atpused_array(k)=ATPused;apop_array(k)=apop;
    eda_array(k)=eda;cda_array(k)=cda;vda_array(k)=vda;ros_mit_array(k)=ROS_mit;
    
    calb_array(k)=Calb;cam_array(k)=Cam;caer_array(k)=Ca_er;camt_array(k)=Ca_mt;
    
    ATP_array(k)=ATP;LAC_array(k)=LAC;PYR_array(k)=PYR;GAP_array(k)=GAP;GSH_array(k)=GSH;
    F6P_array(k)=F6P;F26P_array(k)=F26P;PCr_array(k)=PCr;NADPH_array(k)=NADPH;
    
    ROS_array(k)=ROS;ASYN_array(k)=ASYN;ASYNA_array(k)=ASYNA;ASYNT_array(k)=ASYNT;
    ASYNG_array(k)=ASYNG;LB_array(k)=LB;
    LDOPA_array(k)=LDOPA;
    
%     phier=Ca_i-Ca_er;
%     phimt=Ca_i-Ca_mt;
%     
%     phi_er(k)=phier;
%     phi_mt(k)=phimt;
    
    disp(k*dt)
end

phi_er=log(cai_array./caer_array);
phi_mt=log(cai_array./camt_array);

%[snc_firings1]=ConvertAPtoST(snc_firings,1);
snc_firings1=snc_firings;

base1=1/2;
sncfrequency=size(snc_firings1,1)/(2*base1.*dt.*(Ttime).*1e-3);

%% plot
sec=0.001;
fig1=figure(1);
set(fig1, 'Position', [5, 50, 1920, 955]);
sizz=15;
subplot(311)
set(gca,'fontsize',sizz);
%plot(sec*dt*(snc_firings1(:,1)),snc_firings1(:,2),'k.','MarkerSize',10);
% xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('# of neurons','fontweight','bold')
tit=strcat('SNc firings (Freq = ',num2str(sncfrequency),' Hz)');
title(tit,'fontsize',sizz,'fontweight','bold')
xlim([0 sec*dt*Ttime]);
% ylim([-90 90]);
subplot(312)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(V_snc_array)),V_snc_array,'b')
% xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('V (mV)','fontweight','bold')
title('Membrane potential','fontsize',sizz,'fontweight','bold')
subplot(313)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(cai_array)),cai_array,'r')
xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('Ca^{2+} conc. (mM)','fontweight','bold')
title('Ca^{2+} conc.','fontsize',sizz,'fontweight','bold')
% f1=strcat('SNc_',filename);
% saveas(fig1,f1,'png');

fig2=figure(2);
set(fig2, 'Position', [5, 50, 1920, 955]);
sizz=15;
subplot(411)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(ros_mit_array)),ros_mit_array,'r')
% xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('ROS_{mit}','fontweight','bold')
title('ROS_{mit}','fontsize',sizz,'fontweight','bold')
subplot(412)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(apop_array)),apop_array,'r')
% xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('Apoptosis signal','fontweight','bold')
title('Apoptosis signal','fontsize',sizz,'fontweight','bold')
subplot(413)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(atpused_array)),atpused_array,'k')
% xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('ATPused (mM)','fontweight','bold')
title('ATPused','fontsize',sizz,'fontweight','bold')
subplot(414)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(eda_array)),eda_array,'b')
xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('eDA conc. (mM)','fontweight','bold')
title('eDA conc.','fontsize',sizz,'fontweight','bold')
% f2=strcat('Apop_',filename);
% saveas(fig2,f2,'png');
% 
fig3=figure(3);
set(fig3, 'Position', [5, 50, 1920, 955]);
sizz=15;
subplot(411)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(calb_array)),calb_array,'r')
% xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('Calb conc. (mM)','fontweight','bold')
title('Calb','fontsize',sizz,'fontweight','bold')
subplot(412)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(cai_array)),cai_array,'k')
% xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('Ca_i conc. (mM)','fontweight','bold')
title('Ca_i','fontsize',sizz,'fontweight','bold')
subplot(413)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(caer_array)),caer_array,'k')
% xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('ER Ca^{2+} conc. (mM)','fontweight','bold')
title('ER Calcium','fontsize',sizz,'fontweight','bold')
refline([0 mean(caer_array)]);
fh=strcat(num2str(mean(caer_array)));legend(fh);
subplot(414)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(camt_array)),camt_array,'k')
xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('MT Ca^{2+} conc. (mM)','fontweight','bold')
title('MT Calcium','fontsize',sizz,'fontweight','bold')
refline([0 mean(camt_array)]);
fh=strcat(num2str(mean(camt_array)));legend(fh);
% f3=strcat('Ca_',filename);
% saveas(fig3,f3,'png');
%
fig4=figure(4);
set(fig4, 'Position', [5, 50, 1920, 955]);
sizz=10;
subplot(311)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(F6P_array)),F6P_array,'r')
% xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('F6P conc. (mM)','fontweight','bold')
title('F6P conc.','fontsize',sizz,'fontweight','bold')
% ylim([0.09 1.11])
subplot(312)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(F26P_array)),F26P_array,'r')
% xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('F26P conc. (mM)','fontweight','bold')
title('F26P conc.','fontsize',sizz,'fontweight','bold')
% ylim([2.3 2.5])
subplot(313)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(GAP_array)),GAP_array,'r')
% ylim([0.38 0.410])
xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('GAP conc. (mM)','fontweight','bold')
title('GAP conc.','fontsize',sizz,'fontweight','bold')
% f4=strcat('EM1_',filename);
% saveas(fig4,f4,'png');

fig5=figure(5);
set(fig5, 'Position', [5, 50, 1920, 955]);
sizz=10;
subplot(411)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(PYR_array)),PYR_array,'r')
% ylim([0.8 1.1]*1e-3)
% xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('PYR conc. (mM)','fontweight','bold')
title('PYR','fontsize',sizz,'fontweight','bold')
subplot(412)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(LAC_array)),LAC_array,'r')
% ylim([0.9 1.1]*1e-3)
xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('LAC conc. (mM)','fontweight','bold')
title('LAC','fontsize',sizz,'fontweight','bold')
subplot(413)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(ATP_array)),ATP_array,'r')
% xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('ATP conc. (mM)','fontweight','bold')
title('ATP conc.','fontsize',sizz,'fontweight','bold')
subplot(414)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(PCr_array)),PCr_array,'r')
xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('PCr conc. (mM)','fontweight','bold')
title('PCr conc.','fontsize',sizz,'fontweight','bold')
% f5=strcat('EM2_',filename);
% saveas(fig5,f5,'png');

fig6=figure(6);
set(fig6, 'Position', [5, 50, 1920, 955]);
sizz=10;
subplot(611)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(ROS_array)),ROS_array,'r')
% xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('ROS conc. (mM)','fontweight','bold')
title('ROS conc.','fontsize',sizz,'fontweight','bold')
% ylim([0.09 1.11])
subplot(612)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(ASYN_array)),ASYN_array,'r')
% xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('\alpha-syn conc. (mM)','fontweight','bold')
title('\alpha-syn conc.','fontsize',sizz,'fontweight','bold')
subplot(613)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(ASYNA_array)),ASYNA_array,'r')
% xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('\alpha-syn* conc. (mM)','fontweight','bold')
title('\alpha-syn* conc.','fontsize',sizz,'fontweight','bold')
subplot(614)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(ASYNT_array)),ASYNT_array,'r')
% xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('\alpha-syn-tag conc. (mM)','fontweight','bold')
title('\alpha-syn-tag conc.','fontsize',sizz,'fontweight','bold')
subplot(615)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(ASYNG_array)),ASYNG_array,'r')
% xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('\alpha-syn-agg conc. (mM)','fontweight','bold')
title('\alpha-syn-agg conc.','fontsize',sizz,'fontweight','bold')
subplot(616)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(LB_array)),LB_array,'r')
xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('LB conc. (mM)','fontweight','bold')
title('LB conc.','fontsize',sizz,'fontweight','bold')
% f6=strcat('PDP_',filename);
% saveas(fig6,f6,'png');

fig7=figure(7);
set(fig7, 'Position', [5, 50, 1920, 955]);
sizz=10;
subplot(211)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(NADPH_array)),NADPH_array,'r')
% xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('NADPH conc. (mM)','fontweight','bold')
title('NADPH conc.','fontsize',sizz,'fontweight','bold')
subplot(212)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(GSH_array)),GSH_array,'r')
xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('GSH conc. (mM)','fontweight','bold')
title('GSH conc.','fontsize',sizz,'fontweight','bold')
% f7=strcat('Add_',filename);
% saveas(fig7,f7,'png');
% %
fig8=figure(8);
set(fig8, 'Position', [5, 50, 1920, 955]);
sec=0.001;
sizz=10;
subplot(411)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(cda_array)),cda_array,'r')
% xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('cDA conc. (mM)','fontweight','bold')
title('cDA conc.','fontsize',sizz,'fontweight','bold')
% ylim([0.09 1.11])
refline([0 mean(cda_array)]);
fh=strcat(num2str(mean(cda_array)));legend(fh);
subplot(412)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(vda_array)),vda_array,'r')
% xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('vDA conc. (mM)','fontweight','bold')
title('vDA conc.','fontsize',sizz,'fontweight','bold')
% ylim([2.3 2.5])
refline([0 mean(vda_array)]);
fh=strcat(num2str(mean(vda_array)));legend(fh);
subplot(413)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(eda_array)),eda_array,'r')
% ylim([0.38 0.410])
xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('eDA conc. (mM)','fontweight','bold')
title('eDA conc.','fontsize',sizz,'fontweight','bold')
refline([0 mean(eda_array)]);
fh=strcat(num2str(mean(eda_array)));legend(fh);
subplot(414)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(LDOPA_array)),LDOPA_array,'r')
% ylim([0.38 0.410])
xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('LDOPA conc. (mM)','fontweight','bold')
title('LDOPA conc.','fontsize',sizz,'fontweight','bold')
refline([0 mean(LDOPA_array)]);
fh=strcat(num2str(mean(LDOPA_array)));legend(fh);
% f8=strcat('DA_',filename);
% saveas(fig8,f8,'png');

fig10=figure(10);
set(fig10, 'Position', [5, 50, 1920, 955]);
sec=0.001;
sizz=10;
%subtitle('Energy consumption in different cellular processes')
subplot(411)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(V_id)),V_id,'r')
% xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('E_{id}','fontweight','bold')
title('Ion dynamics.','fontsize',sizz,'fontweight','bold')
% ylim([0.09 1.11])
refline([0 mean(V_id)]);
fh=strcat(num2str(mean(V_id)));legend(fh);
subplot(412)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(V_dp)),V_dp,'r')
% xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('E_{dp}','fontweight','bold')
title('Dopamine packing','fontsize',sizz,'fontweight','bold')
% ylim([2.3 2.5])
refline([0 mean(V_dp)]);
fh=strcat(num2str(mean(V_dp)));legend(fh);
subplot(413)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(V_rel)),V_rel,'r')
% xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('E_{vr}','fontweight','bold')
title('Vesicle recycling','fontsize',sizz,'fontweight','bold')
% ylim([2.3 2.5])
refline([0 mean(V_rel)]);
fh=strcat(num2str(mean(V_rel)));legend(fh);
subplot(414)
set(gca,'fontsize',sizz);
plot(sec*dt*(1:numel(V_er)),V_er,'r')
% ylim([0.38 0.410])
xlabel('Time (sec)','fontsize',sizz,'fontweight','bold')
ylabel('E_{er}','fontweight','bold')
title('ER','fontsize',sizz,'fontweight','bold')
refline([0 mean(V_er)]);
fh=strcat(num2str(mean(V_er)));legend(fh);
% f10=strcat('ATPothers_',filename);
% saveas(fig10,f10,'png');

% save(filename)

toc
% end