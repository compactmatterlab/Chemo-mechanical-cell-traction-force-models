% Esteban Vazquez-Hidalgo
% last updated 07.13.2021
% tractionForce.m calls:
% 1. updateStatus.m to update the status of allintegrins and filaments
% 2. myosin.m to update the state of all motors for all filaments
% 3. calculatexx.m to find the total displacement by myosin on a filament
options = optimset('Display','off');

while t < runtime
    t = t+1;
    
    for jj = 1:nactin
        updateStatus
    end
%     
    for jj = 1:nactin
        myosin 
    end
    
    for jj = 1:nactin
        calculatexx
    end
end

figure
tiledlayout(2,2)
nexttile
plot(Force')
hold on
plot(mean(Force),'r','LineWidth',2)
title('Slip - kPa 10, WB 10')

nexttile
title('Slip')
imagesc(Force)
colorbar

nexttile
title('Slip')
imagesc(sf_att)

nexttile
title('Slip')
plot(sum(sf_att))

close all


