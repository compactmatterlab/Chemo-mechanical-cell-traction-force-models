% Esteban Vazquez-Hidalgo
% last updated 07.13.2021
% updateStatus.m
% update status integrin activation, integrin attachment, and filament
% status to time (t) based on status at (t-1). we update to time (t) so
% that downstream we can determine which filaments are generating force and
% which aren't.
% update status (t) based on status (t-1)
% updated status used to calculate delta (t) based on status (t)
if t == 1
    int_act(jj,t) = int_act(jj,t);
    int_att(jj,t) = int_att(jj,t);
    sf_att(jj,t) = sf_att(jj,t);
else
    tempvec = mapps(jj,:);
    nbors = [...
        int_act(jj,t-1), int_act(jj,t);
        int_att(jj,t-1), int_att(jj,t);
        sf_att(jj,t-1), sf_att(jj,t)];
    
    if int_up(jj,t) == 0
        int_act(jj,t) = int_act(jj,t-1);
        int_att(jj,t) = int_att(jj,t-1);
%         mixmat3(jj,t) = mixmat3(jj,t-1);
    end
    
    if sf_up(jj,t) == 0
        int_act(jj,t) = int_act(jj,t-1);
        sf_att(jj,t) = sf_att(jj,t-1);
%         mixmat2(jj,t) = mixmat2(jj,t-1);
    end
%     r1(jj,t) = r1(jj,t-1);
%     r3(jj,t) = r3(jj,t-1);
    % %     mixmat3(jj,t) = mixmat3(jj,t-1);
    % if status at t == 0, integrin/filament need to update. if
    % int_act(jj,t), (5), and (6) == 0, activate current filament by neighbors
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nbors(1)==1
        for kk = 1:8
            % if neighboring integrin is t-1==0 and t==0, then it is
            % avaialbe to be activated. if t-1==1, then it is not
            % avaialable to be activated, and anything that happens to
            % that neighbor should happen when it is jj
            if rand <= actproba && int_act(tempvec(kk),t-1) == 0 && int_up(tempvec(kk),t) == 0
                int_act(tempvec(kk),t) = 1;
                int_up(tempvec(kk),t) = 1;
            end
        end
    end
    
    if int_up(jj,t) == 0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if nbors(1) == 0
            for kk = 1:8
                if int_act(jj,t) == 0
                    if rand <= actproba && int_act(tempvec(kk),t-1) == 1
                        int_act(jj,t) = 1;
                        int_up(jj,t) = 1;
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif nbors(1) == 1 && nbors(2) == 0
            if rand <= attprob
                int_att(jj,t) = 1;
                int_up(jj,t) = 1;
            else
                if rand <= actprobd
                    int_act(jj,t) = 0;
                    int_up(jj,t) = 1;
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif nbors(1) == 1 && nbors(2) == 1 && nbors(3) == 0
            if rand <= sfproba && sf_up(jj,t)==0
                sf_att(jj,t) = 1;
                sf_up(jj,t) = 1;
                int_up(jj,t) = 1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif nbors(1) == 1 && nbors(2) == 1 && nbors(3) == 1
            Force(jj,t-1) = delta(jj,t-1) * k_spring;
            koff(jj,t-1) = ks0*exp(Force(jj,t-1)/Fb);
%             koff(jj,t-1) = (FitA*exp((-Force(jj,t-1)*epsilon)/(kbT))+(FitB*exp(Force(jj,t-1)*epsilon/(kbT))+FitC*exp((-Force(jj,t-1)*epsilon)/(kbT)))^(-1))^(-1);
            pprob(jj,t-1) = 1-exp(-koff(jj,t-1) * delta_t);% probability of detaching
            if rand <= pprob(jj,t-1)
                sf_att(jj,t) = 0;
                int_att(jj,t) = 0;
                sf_up(jj,t) = 1;
                int_up(jj,t) = 1;
            end
        end
    end
end
% end
% 
%             else
%                 if (Force(jj,t-1) > 1 && Force(jj,t-1) <= 15) || Force(jj,t-1) > 21
%                     for mm = 1:8
%                         if rand <= taluprob(1) && sf_att(tempvec(mm),t-1) == 0 && sf_up(tempvec(mm),t) == 0 ...
%                                 && r1(jj,t-1) < r1max && r1(jj,t) < r1max
%                             
%                             disp('rec 1')
%                             sf_att(jj,t) = 2;
%                             jj;
%                             sf_att(tempvec(mm),t) = 2;
%                             tempvec(mm);
%                             int_up(jj,t) = 1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             sf_up(jj,t) = 1;
%                             sf_up(tempvec(mm),t) = 1;
%                             mixmat2(jj,t) = jj;
%                             mixmat2(tempvec(mm),t) = jj;
%                             r1(jj,t) = r1(jj,t) + 1;
%                         end
%                     end
%                 elseif Force(jj,t-1) > 15 && Force(jj,t-1) <= 21
%                     for mm = 1:8
%                         if rand <= taluprob(2) && int_att(tempvec(mm),t-1) == 0 && int_up(tempvec(mm),t) == 0 ...
%                                 && r3(jj,t-1) < r3max && r3(jj,t) < r3max
%                             
%                             disp('branch  3')
%                             sf_att(jj,t) = 3;
%                             int_att(jj,t) = 3;
%                             jj
%                             int_att(tempvec(mm),t) = 3;
%                             tempvec(mm)
%                             sf_up(jj,t) = 1;
%                             int_up(jj,t) = 1;
%                             int_up(tempvec(mm),t) = 1;
%                             mixmat3(jj,t) = jj;
%                             mixmat3(tempvec(mm),t) = jj;
%                             figure;plot(mixmat3(:,t));
%                             close all
%                             r3(jj,t) = r3(jj,t) + 1;
%                         end
%                     end
%                 end
%             end
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         elseif nbors(1) == 1 && nbors(2) == 1 && nbors(3) == 2
%             disp('                                         1 1 2')
%             if mixmat2(jj,t)== jj
%                 hh = find(mixmat2(:,t-1)== jj);
%                 Force(jj,t-1) = delta(jj,t-1) * k_spring;
%                 koff(jj,t-1) = (FitA*exp((-Force(jj,t-1)*epsilon)/(kbT))+(FitB*exp(Force(jj,t-1)*epsilon/(kbT))+FitC*exp((-Force(jj,t-1)*epsilon)/(kbT)))^(-1))^(-1);
%                 pprob(jj,t-1) = 1-exp(-koff(jj,t-1) * delta_t);% probability of detaching
%                 if rand <= pprob(jj,t-1)
%                     disp('break bond')
%                     
%                     sf_att(hh(:),t) = 0;
%                     sf_up(hh(:),t) = 1;
%                     int_up(jj,t) = 1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     mixmat2(hh(:),t) = 0;
%                     int_att(jj,t) = 0;
%                     r1(jj,t) = 0;
%                     
%                 else
%                     if (Force(jj,t-1) > 1 && Force(jj,t-1) <= 15) || Force(jj,t-1) > 21
%                         for mm = 1:8
%                             if rand <= taluprob(1) && sf_att(tempvec(mm),t-1) == 0 && sf_up(tempvec(mm),t) == 0 ...
%                                     && r1(jj,t-1) < r1max && r1(jj,t) < r1max
%                                 disp('rec!  r1')
%                                 sf_att(jj,t) = 2;
%                                 jj;
%                                 sf_att(tempvec(mm),t) = 2;
%                                 tempvec(mm);
%                                 int_up(jj,t) = 1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 sf_up(jj,t) = 1;
%                                 sf_up(tempvec(mm),t) = 1;
%                                 mixmat2(jj,t) = jj;
%                                 mixmat2(tempvec(mm),t) = jj;
%                                 r1(jj,t) = r1(jj,t) + 1;
%                             end
%                         end
%                     end
%                 end
%             end
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         elseif nbors(1) == 1 && nbors(2) == 3 && nbors(3) == 3
%             disp('                    1 3 3')
%             Force(jj,t-1) = delta(jj,t-1) * k_spring;
%             koff(jj,t-1) = (FitA*exp((-Force(jj,t-1)*epsilon)/(kbT))+(FitB*exp(Force(jj,t-1)*epsilon/(kbT))+FitC*exp((-Force(jj,t-1)*epsilon)/(kbT)))^(-1))^(-1);
%             pprob(jj,t-1) = 1-exp(-koff(jj,t-1) * delta_t);% probability of detaching
%             hh = find(mixmat3(:,t-1)== jj);
%             gg = find(hh(:) ~= jj);
%             if rand <= pprob(jj,t-1)
%                 disp('break bond')
%                 int_att(jj,t) = 1;
%                 int_att(hh(gg),t) = 0;
%                 sf_att(jj,t) = 1;
%                 sf_up(jj,t) = 1;
%                 int_up(hh(:),t) = 1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 mixmat3(hh(:),t) = 0;
%                 r3(jj,t) = 0;
%             else
%                 int_att(hh(:),t) = int_att(hh(:),t-1);
%                 int_up(hh(:),t) = 1;
%                 sf_att(jj,t) = sf_att(jj,t-1);
%                 sf_up(jj,t) = 1;
%                 mixmat3(hh(:),t) = mixmat3(hh(:),t-1);
%                 r3(jj,t) = r3(jj,t-1);
%                 
%             end
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             child, find parent,
%         elseif nbors(1) == 1 && nbors(2) == 3 && (nbors(3) == 0 || nbors(3) == 2)
%             disp('                                         1 3 ?')
%             Force(jj,t-1) = delta(jj,t-1) * k_spring;
%             koff(jj,t-1) = (FitA*exp((-Force(jj,t-1)*epsilon)/(kbT))+(FitB*exp(Force(jj,t-1)*epsilon/(kbT))+FitC*exp((-Force(jj,t-1)*epsilon)/(kbT)))^(-1))^(-1);
%             pprob(jj,t-1) = 1-exp(-koff(jj,t-1) * delta_t);% probability of detaching
%             par = mixmat3(jj,t-1);
%             hh = find(mixmat3(:,t-1) == par);
%                             
%             if rand <= pprob(jj,t-1)
% 
%                 int_att(par,t) = 1;
%                 int_att(jj,t) = 0;
%                 sf_att(par,t) = 1;
%                 sf_up(par,t) = 1;
%                 int_up(hh(:),t) = 1;
%                 mixmat3(hh(:),t) = 0;
%                 r3(par,t) = 0;
%             else
%                 int_att(hh(:),t) = int_att(hh(:),t-1);
%                 int_up(hh(:),t) = 1;
%                 sf_att(jj,t) = sf_att(jj,t-1);
%                 sf_up(jj,t) = 1;
%                 mixmat3(hh(:),t) = mixmat3(hh(:),t-1);
%                 r3(par,t) = r3(par,t-1);
%             end
%             
%         end
%     end
% end
% 
