% Esteban Vazquez-Hidalgo
% last update 02 .27.2021
% rates, initial conditions, constants, etc...
rows = 10; % number of rows
cols = rows; % number of columms
actlen = 3000; % actin length in nm                                         
nactin = rows*cols; % number of actin filaments                            
nmotors = 60; % number of myosin heads                                     
seconds = 100;
realtime = seconds*1000;% converted to  ms
delta_t = 1;
runtime = realtime/delta_t; % simulations time in ms
mappingMatrix % call to get matrix of neighbors for every element
%ms = 1e-3; % seconds to miliseconds conversion
ms = 1e-3;
ks0 = 0.1; % unloaded off-rate from Bens' paper... find original source
ks0 = ks0*ms;       
Fb = 5; % pN, calculated from kBt/epsilon 
act_ka = 23*ms;  
act_kd = 0.01*ms; 
sf_kon = 0.65*ms; % vinculin binding rate limiting
sf_koff = 8*ms;              % filament deattachment rate s^-1           
int_bind = 1e3*ms;                 % integrin-ligand attachment rate s^-1   
vincb = 0.65*ms; % s^-1, vinculin->talin binding rate 
trates = [0.018 99 2.5e-5 4.2e-6 1.7e-8]*ms; % ms^-1, converted to  talin unfolding rates DOI:10.1038/ncomms11966
FitA = 3309; 			% slip ms  
FitB = 0.0000003942; % catch ms^-1 
FitC = 0.05819;	% slip ms^-1 
k12 = 0.00014;% in ms, rate of transition from state 1 to state 2, 
k23 = 0.007;% in ms, rate of transition from state 2 to state 3
k32 = 0.011;% in ms, rate of transition from state 3 to state 2
k34 = 0.00016;% in ms, rate of transition from state 3 to state 4
k41 = 0.028;% in ms, rate of transition from state 4 to state 1
tao = 1; % previous time stamp
epsilon = 0.74; % value for dissociation probability                       % Ben
km = 4; % pN/nm, stiffness of motor stalk
F = 0; % applied force in pN
t = 0; % time ms
y = 5.3; % myosin step size nm
bsd = 2.7;  % binding site distance on actin nm
toler = 1.017;%brownian motion in nm
kb = 1.3806e-2;  % Boltzmanm Constant
T = 300;% Kelvin
kbT = kb*T; % JK^1
mu = 14.3; % myosin spacing (nm)
sigma = 0; % standard deviation
drag = 6*10^-4; % pN*ms/nm
len_scale = 0.1; % um
k_spring_vals = [1 2 10 20 100 200 1000 10000]; % substrate stiffness kPa% experimental ratios
k_spring_vals = k_spring_vals*len_scale;
WB = [0.1 0.2 1 2 10];
ks152 = .01;% phosphorylated to dephosphorylated
% -------------------------------------------------------------------------
%                    constant probabilities
%--------------------------------------------------------------------------
vincprob = 1-(exp(-vincb)); % probability of vinculin binding
taluprob = 1-(exp(-trates)); % talin unfolding probablility
actproba = 1-exp(-act_ka*delta_t); % probability of integrin activating
actprobd = 1-exp(-act_kd*delta_t); % probability of integrin deactivating
attprob = 1-exp(-int_bind*delta_t); % probability of integrin attaching to ligand
sfproba = 1-exp(-sf_kon*delta_t); % probability of sf attaching
sfprobd = 1-exp(-sf_koff*delta_t); % probability of sf detaching
%--------------------------------------------------------------------------
%                vectors and matrices for storage
%--------------------------------------------------------------------------
pprob = zeros(nactin,runtime,1);
r1 = zeros(nactin,runtime,1); % regions of talin. max 1 vbs
r1max = 1*2;
r1max = 8;
r2 = zeros(nactin,runtime,1); % max 4 vbs 
r2max = 4*2;
r3 = zeros(nactin,runtime,1); % max 1 ibs or 5 vbs. setting this as an ibs 
r3max = 1;
r4 = zeros(nactin,runtime,1); % max 2 vbs
r4max = 2*2;
int_act = zeros(nactin,runtime ,1); % activated integrin matrix
int_act(1) = 1;
int_att = zeros(nactin, runtime,1); % attached integrin matrix
int_att(1) = 1;
sf_att = zeros(nactin, runtime,1); % stress fiber attached matrix
int_up = zeros(nactin,runtime,1); % filament update matrix
sf_up = zeros(nactin,runtime,1); % filament update matrix
sf_recruit = zeros(nactin, runtime,1); % stress fiber recruited matrix
sf_branch = zeros(nactin, runtime,1); % stress fiber branched matrix
mixmat2 = zeros(nactin, runtime, 1); % matrix for filament recruitment
mixmat3 = zeros(nactin, runtime, 1); % matrix for integrin branching
Force = zeros(nactin, runtime,1);
n12=zeros(nactin,runtime,1); % count of 1 to 2 transitions per run
n34=zeros(nactin,runtime,1); 
N = zeros(nactin,runtime,1);% count of state 1 + state 4
nn = zeros(nactin,runtime,1); % count of state 1 + state 4 + state 3
myo_pos = zeros(nmotors,nactin,1); % zero matrix for myosin position matrix
%  motorhead spacing exponentially distributed spacing of motor heads
%  first motor at (0,0)
for i=1:nmotors
    for j = 1:nactin
        if i==1
            myo_pos(i,j)=0;
        else
            myo_pos(i,j)=myo_pos(i-1,1)+exprnd(mu);
        end
    end
end
myo_pos = myo_pos';
steps = zeros(nactin,nmotors,1);% matrix of zeros to record number of steps (full cycle) for each motor
change = zeros(nactin,nmotors,1);% matrix of zeros to record number of state changes for each motor
tao = ones(nactin,nmotors);% matrix of ones for previous time stamp
state = ones(nactin,nmotors)*15;     
init_on_motors = randperm(nmotors)';
state(init_on_motors(1:ceil(1*nmotors)),:)=2;
a= state';
[m,nn] = size(a) ;
b = a ;
for i = 1:m
    idx = randperm(nn);
    b(i,idx) = a(i,:);
end
state = b';
clear a b i idx
prev_state = state; % previous state of motor
delta = zeros(nactin,runtime,1); % total distance each motor moved for the given runtime
xx = zeros(nactin,runtime,1);% distance the motor moved at a specific time point
strain = zeros(nactin,nmotors,1); % strain created by attached motors
 