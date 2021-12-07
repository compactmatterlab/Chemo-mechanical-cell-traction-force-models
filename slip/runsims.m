% Esteban Vazquez-Hidalgo
% last update 07.27.2021
clear all; %close all
datetime
tic
params
WBratio = WB(1);% select which ratio to use
k215 = (WBratio^-1)*ks152;% dephosporylated to phosphorylated 
k_spring = k_spring_vals(1); % choose value for substrate stiffness
tractionForce
save(sprintf('a1'))
toc; 
% 
datetime
tic
params
WBratio = WB(2);% select which ratio to use
k215 = (WBratio^-1)*ks152;% dephosporylated to phosphorylated 
k_spring = k_spring_vals(1); % choose value for substrate stiffness
tractionForce
save(sprintf('a2'))
toc; 
% 
datetime
tic
params
WBratio = WB(3);% select which ratio to use
k215 = (WBratio^-1)*ks152;% dephosporylated to phosphorylated 
k_spring = k_spring_vals(1); % choose value for substrate stiffness
tractionForce
save(sprintf('a3'))
toc; 

datetime
tic
params
WBratio = WB(4);% select which ratio to use
k215 = (WBratio^-1)*ks152;% dephosporylated to phosphorylated 
k_spring = k_spring_vals(1); % choose value for substrate stiffness
tractionForce
save(sprintf('a4'))
toc; 

datetime
tic
params
WBratio = WB(5);% select which ratio to use
k215 = (WBratio^-1)*ks152;% dephosporylated to phosphorylated 
k_spring = k_spring_vals(1); % choose value for substrate stiffness
tractionForce
save(sprintf('a5'))
toc; 