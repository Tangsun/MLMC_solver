%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FILENAME:                    MLMC_IPOC_11_03_2019.m
%  WRITTEN BY:                  Hang Yang 
%                               University of Michigan, Ann Arbor
%
%  COMMENTS ON CODE STRUCTURE:  
%
%  ORIGINAL DATE WRITTEN:       17 October 2019 - University of Michigan
%  DATE OF LAST MODIFICATION:   28 October 2019
%  MODIFICATION HISTORY:
%  DATE:                        COMMENT:
%  03 November, 2019            Evaluate a quadratic cost function rather
%                               than the states at the last time step
%                              
%  PROBLEMS/MODIFICATIONS FOR FUTURE WORK: 
%
%%%%%%%%%%%%%%%%% START OF HANG YANG's CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc

%% Parameter set-up
% Need postprocessing?
PostPro = input('Need postprocessing? (Y/N)', 's');
Yes = 'Y'; No = 'N';
if PostPro == Yes
    disp('Post processing turned on');
elseif PostPro == No
    disp('Post processing turned off');
else
    disp('Error! Please answer Y for Yes and N for No');
end

T = 1;               %prediction horizon
d = 4;               %dimensional of the state space
sigma = 0.2*eye(d);  %diffusion
X0 = [0 0 pi/8 0]';  %initial condition

N = 5000;            %number of MLMC trails
M = 2;               %refinement factor

t0 = 0;                     % initial time
dt_finest = 2^(-8);         %time step of the high-fidelity model

N_samples = [1200 70 26 8]'; %sample sizes for MLMC

Q = diag([50 0 100 0]);     %stage cost weighting matrix
np = 1;

rng 'default'

% Define the drift
ufunc = @(t, X) zeros(1, size(X, 2)); % control is zero
drift = @(t, x) dyn_IPOC(t, x, ufunc);
g = @(t, sim_out) evaluate_cost(t, sim_out, Q);

N_level = 3;
% Set up time step sizes
dt = zeros(N_level+1, 1);
Nt = zeros(N_level+1, 1);
for i_MLMC = 1:N_level+1
    dt(i_MLMC) = dt_finest*M^(N_level+1-i_MLMC);
    Nt(i_MLMC) = ceil(T/dt(i_MLMC));
end

%% Standard MC
CL = 1;    %assume that the computational cost of evaluating a simple sample at the finest level L is 1
Cost_ML = N_samples(1)/M^3 + N_samples(2)*(1/M^3+1/M^2) + N_samples(3)*(1/M^2+1/M) + N_samples(4)*(1/M+1);
N_EquivMC = ceil(Cost_ML)/CL

Pmc = zeros(N, np);
for ii = 1:N
    [t_out, MC_trajs, W] = MC_EM(X0, t0, dt(end), Nt(end)+1, drift, sigma, N_EquivMC, d);
    Pmc(ii, :) = mean(g(t_out, MC_trajs), 2)';
end

%% 3-Level MLMC
if PostPro == No
    P_MLMC = zeros(N, np);
    for i = 1:N
        P = MLMC_EM(X0, t0, dt, Nt+1, drift, sigma, N_samples, d, M, N_level, g, np);
        P_MLMC(i, :) = sum(P, 2);
    end
    P_MLMC_var = var(P_MLMC);    
    
    %plot results
    h1 = figure(1);
    set(h1,'Name','MLMC and MC Estimates');
    histogram(Pmc, 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor',[0.9290 0.6940 0.1250],'FaceAlpha',0.5, 'Normalization', 'pdf'); hold on
    histogram(P_MLMC,'FaceColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410],'FaceAlpha',0.5, 'Normalization', 'pdf'); hold on
    YL = ylim();
    plot([mean(Pmc) mean(Pmc)],YL,'Color',[0.9290 0.6940 0.1250],'LineWidth',1,'LineStyle','--'); hold on
    plot([mean(P_MLMC) mean(P_MLMC)], YL,'Color',[0 0.4470 0.7410],'LineWidth',1,'LineStyle','--'); hold off
    legend('MC','MLMC');
    xlabel('Cost','interpreter','latex','FontSize',15);
    ylabel('PDF','interpreter','latex','FontSize',15); grid on
    
elseif PostPro == Yes
    P_MLMC = zeros(N, np);
    for i = 1:N
        [P, P_var, P_level_fine, P_var_level_fine, g_out, g_out_coarse] = MLMC_EM_postproc(X0, t0, dt, Nt+1, drift, sigma, N_samples, d, M, N_level, g, np);
        P_MLMC(i, :) = sum(P, 2);
    end
    P_MLMC_var = var(P_MLMC);
    
    % Optimal MLMC Sample Allocation
    CL = 1;
    Var_target = 0.5;
    Nl_opt = Optimal_Allocation(P_var, CL, M, Var_target, N_level)
    
    % Plot Reults
    h1 = figure(1);
    set(h1,'Name','MLMC and MC Estimates');
    histogram(Pmc, 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor',[0.9290 0.6940 0.1250],'FaceAlpha',0.5, 'Normalization', 'pdf'); hold on
    histogram(P_MLMC,'FaceColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410],'FaceAlpha',0.5, 'Normalization', 'pdf'); hold on
    YL = ylim();
    plot([mean(Pmc) mean(Pmc)],YL,'Color',[0.9290 0.6940 0.1250],'LineWidth',1,'LineStyle','--'); hold on
    plot([mean(P_MLMC) mean(P_MLMC)], YL,'Color',[0 0.4470 0.7410],'LineWidth',1,'LineStyle','--'); hold off
    legend('MC','MLMC');
    xlabel('Cost','interpreter','latex','FontSize',15);
    ylabel('PDF','interpreter','latex','FontSize',15); grid on

    h2 = figure(2);
    set(h2, 'Name','Mean of Each Level in MLMC');
    plot(0:3, log10([P(1,1) P_level_fine(1,1) P_level_fine(1,2) P_level_fine(1,3)]),'--o','Color',[0 0.4470 0.7410]); hold on
    plot(1:3, log10([P(1,2) P(1,3) P(1,4)]),'-*','Color',[0 0.4470 0.7410]); hold off
    xlabel('MLMC Level');
    ylabel('Log_{10} Mean');
    legend('X_l','X_l-X_{l-1}'); 
    set(gca,'xtick', [0 1 2 3]);

    h3 = figure(3);
    set(h3, 'Name','Variance of Each Level in MLMC');
    plot(0:3, log10([P_var(1,1) P_var_level_fine(1,1) P_var_level_fine(1,2) P_var_level_fine(1,3)]),'--o','Color',[0 0.4470 0.7410]); hold on
    plot(1:3, log10([P_var(1,2) P_var(1,3) P_var(1,4)]),'-*','Color',[0 0.4470 0.7410]); hold off
    xlabel('MLMC Level');
    ylabel('Log_{10} Variance');
    legend('X_l','X_l-X_{l-1}'); 
    set(gca,'xtick', [0 1 2 3]);
end
%%%%%%%%%%%%%%%%%%%%%% END OF HANG YANG's CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
