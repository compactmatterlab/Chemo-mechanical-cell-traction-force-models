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
    
    for jj = 1:nactin
        myosin 
    end
    
    for jj = 1:nactin
        calculatexx
    end
end


figure(1)
tiledlayout(2,1)
nexttile
imagesc(int_att)
colorbar
title('int-att')
% figure(2)
nexttile
imagesc(sf_att)
colorbar
title('sf-att')
colorbar

figure(2)
tiledlayout(2,1)
nexttile
imagesc(delta)
title('delta')
colorbar
% figure(6)
nexttile
imagesc(Force)
title('Force')
colorbar

figure(3)
tiledlayout(2,1)
nexttile
imagesc(r1)
colorbar
title('r1')
nexttile
imagesc(mixmat2)
colorbar
title('mixmat2')

figure(4)
tiledlayout(2,1)
nexttile
imagesc(r3)
colorbar
title('r3')
nexttile
imagesc(mixmat3)
colorbar
title('mixmat3')