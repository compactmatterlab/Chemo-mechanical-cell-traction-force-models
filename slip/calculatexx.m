% calculatexx.m finds the distance the myosin heads stepped at time (t).
if t == 1
    xx(jj,t) = fsolve(@(x) (-(N(jj,t)*km*y)+(km*sum(strain(jj,:)))+(nn(jj,t)*km*x)+...
        +F+(drag*x/delta_t)),0,options);
else
    
    if int_act(jj,t) == 1 && int_att(jj,t) == 1 && sf_att(jj,t) == 1
        xx(jj,t) = fsolve(@(x) (-(N(jj,t)*km*y)+(km*sum(strain(jj,:)))+(nn(jj,t)*km*x)+ ...
            (k_spring*delta(jj,t-1))+(k_spring*x)+F+(drag*x/delta_t)),0,options);
        delta(jj,t) = delta(jj,t-1) + xx(jj,t);
    else
        xx(jj,t) = fsolve(@(x) (-(N(jj,t)*km*y)+(km*sum(strain(jj,:)))+(nn(jj,t)*km*x)+...
            +F+(drag*x/delta_t)),0,options);
    end
end