% Esteban Vazquez-Hidalgo
% Carly Farris
% Last update 07.12.2021

% myosin.m updates the state of every motor at the current filament at time
% (t) and keeps track of the number of state changes.

% initialize states
state1(jj,t) = 0;
state2(jj,t) = 0;
state3(jj,t) = 0;
state4(jj,t) = 0;

% track number of 3-4 transitions
if t == 1
    n34(jj,t) = 0;
else
    n34(jj,t) = n34(jj,t-1);
end

% loop through motors. only one state per motor per time (t)
% track number of motors in each state
for motors = 1:nmotors
    success = 0;
    % transition from 1->2 or 1->1
    if state(jj,motors) == 1 && success == 0
        strain(jj,motors) = strain(jj,motors)+xx(jj,t-1);
        % strain incorporated into the rate constant
        ks12 = k12*exp((0.5*km*(strain(jj,motors))^2)/(kb*T));
        % this is the probability that it will change state in
        % time delta_t
        a = 1-exp(-ks12*delta_t);
        % if a random variable is less than the probability that
        % it will change state, it will change
        if rand < a
            prev_state(jj,motors) = 1;
            % next state the motor will go to
            state(jj,motors) = 2;
            % stores the time when the motor passes to next step
            tao(jj,motors) = t;
            % stores the state change count
            change(jj,motors) = change(jj,motors)+1;
            % this prevents the motor from going to a new state in
            % the same time point
            success = 1;
            strain(jj,motors) = 0;
            n12(jj) = n12(jj)+1;
        end
    end
    % transition from 15-> 2 or 15->15
    if state(jj,motors) == 15 && success == 0
        e1 = (1-exp(-ks152*delta_t));
        if rand < e1
            prev_state(jj,motors) = 15;
            state(jj,motors) = 2;
            tao(jj,motors) = t;
            success = 1;
            strain(jj,motors) = 0;
        end
    end
    
    if state(jj,motors) == 2 && success == 0
        % transistion from 2->3 or 2->15 or 2->2
        ks23 = k23*exp((-0.5*km*(strain(jj,motors))^2)/(kb*T));
        % this is the probability that it will change state in time delta_t
        b1 = (1-exp(-ks23*delta_t));
        e2 = (1-exp(-k215*delta_t));
        u1 = rand;
        u2 = b1+e2;
        bb1 = b1/u2;
        ee2 = e2/u2;
        % if u1 < b1 + e2
        if u2 < 1
            if t==1
                % if a random variable is less than the probability
                % that it will change state, it will change but only
                % if actin binding site and myosin head are aligned
                % within a tolerance
                if ((rem((0+myo_pos(jj,motors)),bsd)<toler) || ...
                        ((bsd-rem((0+myo_pos(jj,motors)),bsd))<toler))...
                        && u1 < b1
                    prev_state(jj,motors) = 2;
                    % next state the motor will go to
                    state(jj,motors) = 3;
                    % stores the time when the motor passes to next step
                    tao(jj,motors) = t;
                    % stores the state change count
                    change(jj,motors) = change(jj,motors)+1;
                    % this prevents the motor from going to a new state
                    % in the same time point
                    success = 1;
                    if rem((0+myo_pos(jj,motors)),bsd)<toler
                        strain(jj,motors) = rem((0+myo_pos(jj,motors)),bsd);
                    elseif (bsd-rem((0+myo_pos(jj,motors)),bsd))<toler
                        strain(jj,motors) = -(bsd-rem((0+myo_pos(jj,motors)),bsd));
                    end
                elseif b1 < u1 && u1 < u2
                    prev_state(jj,motors) = 2;
                    state(jj,motors) = 15;
                    tao(jj,motors) = t;
                    change(jj,motors) = change(jj,motors)+1;
                    success = 1;
                elseif u1 > u2
                    prev_state(jj,motors) = 2;
                    state(jj,motors) = 2;
                    tao(jj,motors) = t;
                    success = 1;
                end
            else
                % determine state change with binding distance
                % if a random variable is less than the probability
                % that it will change state, it will change but only
                % if actin binding site and myosin head are aligned
                % within a tolerance
                if ((rem((delta(jj,t-1)+myo_pos(jj,motors)),bsd)<toler) || ...
                        ((bsd-rem((delta(jj,t-1)+myo_pos(jj,motors)),bsd))<toler)) ...
                        && u1 < b1
                    prev_state(jj,motors) = 2;
                    % next state the motor will go to
                    state(jj,motors) = 3;
                    % stores the time when the motor passes to next step
                    tao(jj,motors) = t;
                    % stores the state change count
                    change(jj,motors) = change(jj,motors)+1;
                    % this prevents the motor from going to a new state
                    % in the same time point
                    success = 1;
                    if rem((delta(jj,t-1)+myo_pos(jj,motors)),bsd)<toler
                        strain(jj,motors) = rem((delta(jj,t-1)+myo_pos(jj,motors)),bsd);
                    elseif (bsd-rem((delta(jj,t-1)+myo_pos(jj,motors)),bsd))<toler
                        strain(jj,motors) = -(bsd-rem((delta(jj,t-1)+myo_pos(jj,motors)),bsd));
                    end
                elseif b1 < u1 && u1 < u2
                    prev_state(jj,motors) = 2;
                    state(jj,motors) = 15;
                    tao(jj,motors) = t;
                    change(jj,motors) = change(jj,motors)+1;
                    success = 1;
                elseif u1 > u2
                    prev_state(jj,motors) = 2;
                    state(jj,motors) = 2;
                    tao(jj,motors) = t;
                    success = 1;
                end
            end
            % if u1 > b1 + e2
        elseif u2 > 1
            if t == 1
                % if a random variable is less than the probability
                % that it will change state, it will change but only
                % if actin binding site and myosin head are aligned
                % within a tolerance
                if ((rem((0+myo_pos(jj,motors)),bsd)<toler) || ...
                        ((bsd-rem((0+myo_pos(jj,motors)),bsd))<toler))...
                        && u1 <= bb1
                    prev_state(jj,motros) = 2;
                    % next state the motor will go to
                    state(jj,motors) = 3;
                    % stores the time when the motor passes to next step
                    tao(jj,motors) = t;
                    % stores the state change count
                    change(jj,motors) = change(jj,motors)+1;
                    % this prevents the motor from going to a new state
                    % in the same time point
                    success = 1;
                    if rem((0+myo_pos(jj,motors)),bsd)<toler
                        strain(jj,motors) = rem((0+myo_pos(jj,motors)),bsd);
                    elseif (bsd-rem((0+myo_pos(jj,motors)),bsd))<toler
                        strain(jj,motors) = -(bsd-rem((0+myo_pos(jj,motors)),bsd));
                    end
                elseif bb1 < u1 && u1 < bb1 + ee2
                    prev_state(jj,motors) = 2;
                    state(jj,motors) = 15;
                    tao(jj,motors) = t;
                    change(jj,motors) = change(jj,motors)+1;
                    success = 1;
                end
            else
                %  determine state change with binding distance
                % if a random variable is less than the probability
                % that it will change state, it will change but only
                % if actin binding site and myosin head are aligned
                % within a tolerance
                if ((rem((delta(jj,t-1)+myo_pos(jj,motors)),bsd)<toler) || ...
                        ((bsd-rem((delta(jj,t-1)+myo_pos(jj,motors)),bsd))<toler)) ...
                        && u1 <= bb1
                    prev_state(jj,motors) = 2;
                    % next state the  motor will go to
                    state(jj,motors) = 3;
                    % stores the time when the motor passes to next step
                    tao(jj,motors) = t;
                    % stores the state change count
                    change(jj,motors) = change(jj,motors)+1;
                    % this prevents the motor from going to a new state
                    % in the same time point
                    success = 1;
                    if rem((delta(jj,t-1)+myo_pos(jj,motors)),bsd)<toler  
                            strain(jj,motors) = rem((delta(jj,t-1)+myo_pos(jj,motors)),bsd);
                    elseif (bsd-rem((delta(jj,t-1)+myo_pos(jj,motors)),bsd))<toler 
                            strain(jj,motors) = -(bsd-rem((delta(jj,t-1)+myo_pos(jj,motors)),bsd));
                    end
                elseif bb1 < u1 && u1 < bb1 + ee2
                    prev_state(jj,motors) = 2;
                    state(jj,motors) = 15;
                    tao(jj,motors) = t;
                    change(jj,motors) = change(jj,motors)+1;
                    success = 1;
                end%
            end
        end
    end
    %3 -> 4 or 3 -> 2
    %  calculate strain on motors
    if state(jj,motors) == 3 && success == 0
        strain(jj,motors) = strain(jj,motors)+xx(jj,t-1);
        % strain incorporated into the rate constant
        ks34 = k34*exp((km*(strain(jj,motors))*bsd)/(kb*T));
        ks32 = k32;
        % this is the probability that it will change state in
        % time delta_t
        c = (1-exp(-ks34*delta_t));
        % this is the probability that it will change state in
        % time delta_t
        b2 = (1-exp(-ks32*delta_t));
        % random variable 1
        x1 = rand;
        x2 = c + b2;
        cc = c/x2;
        bb2 = b2/x2;
        %  determine if motors will change states
        % if both random variables are less than the probability
        % that it will change state, it will stay in the same state
        if x2 <= 1
            if x1 > x2
                prev_state(jj,motors) = 3;
                % next state the motor will go to
                state(jj,motors) = 3;
                % stores the time when the motor passes to next step
                tao(jj,motors) = t;
                % this prevents the motor from going to a new state in
                % the same time point
                success = 1;
                % if a random variable is less than the probability
                % that it will change state, it will change
            elseif c < x1 && x1 < x2
                prev_state(jj,motors) = 3;
                % next state the motor will go to
                state(jj,motors) = 2;
                % stores the time when the motor passes to next step
                tao(jj,motors) = t;
                % stores the state change count
                change(jj,motors) = change(jj,motors)+1;
                % this prevents the motor from going to a new state in
                % the same time point
                success = 1;
                strain(jj,motors) = 0;
                % if a random variable is less than the probability
                % that it will change state, it will change
            elseif x1 < c
                prev_state(jj,motors) = 3;
                % next state the motor will go to
                state(jj,motors) = 4;
                % stores the time when the motor passes to next step
                tao(jj,motors) = t;
                % stores the state change count
                change(jj,motors) = change(jj,motors)+1;
                % this prevents the motor from going to a new state
                % in the same time point
                success = 1;
                n34(jj,t) = n34(jj,t)+1;
            end%
        elseif x2 > 1
            if x1 <= cc
                prev_state(jj,motors) = 3;
                % next state the motor will go to
                state(jj,motors) = 4;
                % stores the time when the motor passes to next step
                tao(jj,motors) = t;
                % this prevents the motor from going to a new state in
                % the same time point
                success = 1;
                % if a random variable is less than the probability
                % that it will change state, it will change
            elseif cc < x1 && x1 < cc+bb2
                prev_state(jj,motors) = 3;
                % next state the motor will go to
                state(jj,motors) = 2;
                % stores the time when the motor passes to next step
                tao(jj,motors) = t;
                % stores the state change count
                change(jj,motors) = change(jj,motors)+1;
                % this prevents the motor from going to a new state in
                % the same time point
                success = 1;
                strain(jj,motors) = 0;
            end
        end
    end
    %  calculate strain for each motor
    if state(jj,motors) == 4 && success == 0
        strain(jj,motors) = strain(jj,motors)+xx(jj,t-1);
        % not a function of strain
        ks41 = k41*1;
        % this is the probability that it will change state in
        % time delta_t
        % determine if the motor will change states
        % if a random variable is less than the probability that
        % it will change state, it will change
        d = (1-exp(-ks41*delta_t));
        if rand < d
            prev_state(jj,motors)=4;
            % next state the motor will go to
            state(jj,motors) = 1;
            % stores the time when the motor passes to next step
            tao(jj,motors) = t;
            % stores the state change count
            change(jj,motors) = change(jj,motors)+1;
            % stores the step (full completion of cycle) count
            steps(jj,motors) = steps(jj,motors)+1;
            % this prevents the motor from going to a new state in
            % the same time point
            success = 1;
        end%
    end%
    % array of total steps for each motor
    motorstep(jj,motors) = steps(jj,motors);
    % array of total state changes for each motor
    motorstatechange(jj,motors) = change(jj,motors);
    %  count how many heads are in state 1,2,3,4
    if state(jj,motors) == 1
        % count of motors in state 1 for each time point
        state1(jj,t) = state1(jj,t)+1;
    elseif state(jj,motors) == 2
        % count of motors in state 2 for each time point
        state2(jj,t) = state2(jj,t)+1;
    elseif state(jj,motors) == 3
        % count of motors in state 3 for each time point
        state3(jj,t) = state3(jj,t)+1;
    elseif state(jj,motors) == 4 
        % count of motors in state 4 for each time point
        state4(jj,t) = state4(jj,t)+1;
    end
end
% number of motors in states 1, 3, and 4 (number of motors
% attached to actin)
nn(jj,t) = state3(jj,t) + state1(jj,t) + state4(jj,t);
% # of motors in state 1 and 4 (number of motors pulling actin)
N(jj,t) = state1(jj,t) + state4(jj,t);