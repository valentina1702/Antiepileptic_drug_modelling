% EN.580.640 - Final Project
% Sophia Nehs, Valentina Dsouza, Caroline Ghio, Christianne Chua, Shruthi Bare
% Levetiracetam Sim

function [out1,out2,out3,out4,out5] = LEV_sim(D0,p,DoseTimes);
%
% This function runs one simulation of repeated Levetiracetam dosing
%

% % INPUTS % %
% D0 = drug dose (amount)
% p = parameter values only including ka, kc, V, and q
% DoseTimes = the time of dose


% % OUTPUTS % %
% out1 = T = output time vector from solved ode function
% out2 = Y = output Y matrix for the concentration of Levetiracetam in blood, amount in
% gut, and amount in cleared compartment in columns 1,2,and 3 respectively.
% out3 = BalanceD = mass balance for the simulation
% out4 = AUC = area under the concentration curve, metric for drug exposure
% out5 = Ctrough = the concentration value of Levetiracetam in blood
% immediately prior to next dose.

options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);

%% RUN SOLVER FOR EACH DOSE

Ctrough = zeros(length(DoseTimes)-1,1);

for d = 1:length(DoseTimes)-1
    if d == 1 % FIRST DOSE
        % Set initial conditions
        y0 = [0 D0 0]; % y0(1) = mg/L, y0(2) and y0(3) = mg
        tspan = [DoseTimes(d):0.01:DoseTimes(d+1)];
        
        % Run solver
        [T1,Y1] = ode45(@LEV_eqns,tspan,y0,options,p);
        T = T1;
        Y = Y1;
        
        % Record Ctrough value (drug concentration immediately before next dose)
        Ctrough(d) = Y(end,1); % mg/L
        
        % Mass Balance
%         CurrentDrug(:,1) = Y(:,1)*p(3); % current drug in blood compartment (mg)
%         CurrentDrug(:,2) = Y(:,2); % current drug in gut compartment (mg)
%         CurrentDrug(:,3) = Y(:,3); % current drug in cleared compartment (mg)
        DrugIn = p(4)*T + D0*d; % cumulative drug into the system (mg)
%         DrugOut = CurrentDrug(:,3);
%         BalanceD = DrugIn - DrugOut - CurrentDrug(:,1) - CurrentDrug(:,2); 
    else % ALL SUBSEQUENT DOSES
        
        % Set initial conditions
        y0 = Y(end,:); % initialize with final values from previous dosing sim
        y0(2) = y0(2) + D0; % add a new dose to the gut compartment
        tspan = [DoseTimes(d):0.01:DoseTimes(d+1)];
        
        % Run solver
        [T2,Y2] = ode45(@LEV_eqns,tspan,y0,options,p);
        T = [T;T2(2:end)];
        Y = [Y;Y2(2:end,:)];
        
        % Record Ctrough
        Ctrough(d) = Y2(end,1); % mg/L
        
        % Mass Balance
        DrugIn2 = p(4)*T2 + D0*d;
        DrugIn = [DrugIn; DrugIn2(2:end)];
        
    end
end

%% MASS BALANCE

CurrentDrug(:,1) = Y(:,1)*p(3); % current drug in blood compartment (mg)
CurrentDrug(:,2) = Y(:,2);      % current drug in gut compartment (mg)
DrugOut = Y(:,3);

BalanceD = DrugIn - DrugOut - CurrentDrug(:,1) - CurrentDrug(:,2); %(zero = balance)

% Check if there is an imbalance
if mean(BalanceD)>1e-6
    disp('Mass imbalance possible: ');
    disp(BalanceD);
end

%% AUC CALCULATION

AUC = trapz(T,Y(:,1)); % concentrations are mg/L so AUC units = mg/L*hr

%% RETURN OUTPUTS

out1 = T;
out2 = Y;
out3 = BalanceD;
out4 = AUC;
out5 = Ctrough;
end
