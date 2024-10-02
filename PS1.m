% Define parameters

par.alpha = 0.36; % capital share
par.beta = 0.99; % discount factor
par.delta = 0.025;
T = [50, 100, 200];
par.eps = 1e-9;

% steady-state
kstar = ((par.delta - 1 + 1/par.beta) / par.alpha)^(1/(par.alpha - 1));
cstar = kstar^par.alpha - par.delta * kstar;

%k_0 definition
k_0 = 0.9 * kstar;


function [c,k, out] = NCGM(c_0, k_0, T, par)

    % Exogenous productivity shifter (we set A = 1 for all t)
    % Warning : for endogenous productivity shifter, RGM function should be
    % adapted to add the computation of A before the one of (k,c) in the loop 
    A = ones(1,T+1);

    % Warning : A shall be a vector of length T
    % A represents an exogenous productivity shifter
    %if ~isvector(A) 
    %    error(' A shall be a vector of length T');
    %end

    c = zeros(1, T+1);
    k = zeros(1, T+1);

    c(1) = c_0;
    k(1) = k_0;

    out = false;

    for t = 2:T+1
        %k or c should not be negative
        k(t) = (1 - par.delta)*k(t-1) + A(t-1)*k(t-1)^par.alpha - c(t-1);
        c(t) = c(t-1) * par.beta * (A(t) * par.alpha * k(t)^(par.alpha - 1) + 1 - par.delta);
        if k(t) < -par.eps || c(t) < -par.eps % par.eps insted of 0 to avoid looping indefinitely            
            out = true;
            c(t+1:end) = 0;
            k(t+1:end) = 0;
            break;
        end   
    end
end

% test over the steady state
[C, K, out] = NCGM(cstar, kstar, 500, par);
C(501) == cstar
K(501) == kstar
out == false
%it works

function [C, K] = shooting(k_0, kstar, cstar, par, T, maxiter, question)
    
    c_l = 0 ;
    c_u = 10 ;

    error = 1;
    i = 1;

    while abs(error) > par.eps && i <= maxiter
        i = i + 1;
        c_0 = (c_l + c_u)/2;
        [C, K, out] = NCGM(c_0, k_0, T, par);

        if question == 4
            if out == true
                error = -1;
                c_u = c_0;
            else
                error = K(T+1)/k_0;
                c_l = c_0;
            end
        elseif question == 5
            error = K(T+1) - kstar;
            if abs(error) < par.eps 
                    break
            elseif error < 0
                c_u = c_0;
            else
                c_l = c_0;
            end
        else 
            error("Question should either be 4 or 5");
        end 

        
    end
    if i < maxiter
        fprintf(1, "The algorithm converged after %d iterations\n", i)
    else
        fprintf(1, "The algorithm failed to converge after %d iterations\n", maxiter)
    end
end

%----------------------------------%
%%%%         Question 4         %%%% 
%----------------------------------%

% Loop over all values of T
for i = 1:length(T)
    % Run the shooting algorithm
    [C, K] = shooting(k_0, kstar, cstar, par, T(i), 2000, 4);
    
    % Plot the sequences C(t) and K(t)
    figure;
    subplot(2, 1, 1);
    plot(0:T(i)-1, C(1:end-1));
    title(['Consumption Path for T = ', num2str(T(i))]);
    xlabel('Time');
    ylabel('Consumption');
    
    subplot(2, 1, 2);
    plot(0:T(i), K(1:end));
    title(['Capital Path for T = ', num2str(T(i))]);
    xlabel('Time');
    ylabel('Capital');

    % Sauvegarder la figure
    saveas(gcf, ['qu4_figure_T_', num2str(T(i)), '.png']);
end

%----------------------------------%
%%%%         Question 5         %%%% 
%----------------------------------%
T_5 = 200;

% Run the shooting algorithm
[C, K] = shooting(k_0, kstar, cstar, par, T_5, 2000, 5);
    
% Plot the sequences C(t) and K(t)
figure;
subplot(2, 1, 1);
plot(0:T_5, C);
title(['Consumption Path for T = ', num2str(T(i))]);
xlabel('Time');
ylabel('Consumption');

subplot(2, 1, 2);
plot(0:T_5, K);
title(['Capital Path for T = ', num2str(T(i))]);
xlabel('Time');
ylabel('Capital');

saveas(gcf, 'qu5_figure_T_200.png');

