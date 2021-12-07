% Esteban Vazquez-Hidalgo
% last update 07.27.2021
params
WBratio = WB(1);% select which ratio to use
k215 = (WBratio^-1)*ks152;% dephosporylated to phosphorylated 
k_spring = k_spring_vals(1); % choose value for substrate stiffness
tractionForce