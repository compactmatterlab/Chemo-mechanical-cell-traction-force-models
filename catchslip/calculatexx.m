% calculatexx.m finds the distance the myosin heads stepped at time (t). 
if t == 1
    xx(jj,t) = fsolve(@(x) (-(N(jj,t)*km*y)+(km*sum(strain(jj,:)))+(nn(jj,t)*km*x)+...
        +F+(drag*x/delta_t)),0,options);
    delta(jj,t) = 0;
else
    
    if int_act(jj,t) == 1 && int_att(jj,t) == 1 && sf_att(jj,t) == 1
        xx(jj,t) = fsolve(@(x) (-(N(jj,t)*km*y)+(km*sum(strain(jj,:)))+(nn(jj,t)*km*x)+ ...
            (k_spring*delta(jj,t-1))+(k_spring*x)+F+(drag*x/delta_t)),0,options);
        delta(jj,t) = delta(jj,t-1) + xx(jj,t);
        
    elseif int_act(jj,t) == 1 && int_att(jj,t) == 1 && sf_att(jj,t) == 2
        if mixmat2(jj,t)== jj
            hh = find(mixmat2(:,t)==jj);
            tempN = sum(N(hh(:),t));
            tempnn = sum(nn(hh(:),t));
            tempstrain = sum(sum(strain(hh(:),:)));
            tempdelta = sum(delta(hh(:),t-1));
            xx(jj,t) = fsolve(@(x) (-(tempN*km*y)+(km*tempstrain)+(tempnn*km*x)+ ...
                (k_spring*tempdelta)+(k_spring*x)+F+(drag*x/delta_t)),0,options);
            delta(jj,t) = delta(jj,t-1) + xx(jj,t);
        else
            xx(jj,t) = fsolve(@(x) (-(N(jj,t)*km*y)+(km*sum(strain(jj,:)))+(nn(jj,t)*km*x)+...
                +F+(drag*x/delta_t)),0,options);
        end
        
    elseif int_act(jj,t) == 1 && int_att(jj,t) == 3 && sf_att(jj,t) == 3
        hh = find(mixmat3(:,t) == jj);
        if delta(hh(1),t) == 0 && delta(hh(2),t) == 0
            xx(jj,t) = fsolve(@(x) (-(N(jj,t)*km*y)+(km*sum(strain(jj,:)))+(nn(jj,t)*km*x)+ ...
                (2*k_spring*delta(jj,t-1))+(2*k_spring*x)+F+(drag*x/delta_t)),0,options);
            delta(hh(:),t) = delta(hh(:),t-1) + xx(jj,t);
        end
    else
        xx(jj,t) = fsolve(@(x) (-(N(jj,t)*km*y)+(km*sum(strain(jj,:)))+(nn(jj,t)*km*x)+...
            +F+(drag*x/delta_t)),0,options);
    end
end
