% EN.580.640 - Final Project
% Sophia Nehs, Valentina Dsouza, Caroline Ghio, Christianne Chua, Shruthi Bare
% Levetiracetam Driver - Base case and population analysis

close all


%% Baseline model

% Setting up Dose, number of doses, and dose times
D0 = 1500;       % mg
num_doses = 10;  % count
dose_freq = 12;  % hrs
first_doseT = 0; % hrs
last_doseT = num_doses*dose_freq; % hrs, simulation end time (time of last dose + dose freq)
DoseTimes = linspace(first_doseT,last_doseT,(num_doses+1)); % hrs, vector for dose and total simulation timing


%%% Setting parameter values %%%

% Average values for mass, height, and CrCL from Karatza et al., 2020
mass = 79;     % kg
height = 170;  % cm
CrCL = 128.73; % mL/min - creatinine clearance


% Parameter values set by allometric relationships
BSA = sqrt((height*mass)/3600); % m^2, body surface area - (Li et al., 2021)
CrCL = (CrCL/BSA)*1.73;         % mL/min/1.73m^2, CrCL adjusted by BSA - (FDA)
CL = 3.26*(CrCL/139)^0.795;     % L/hr - (Karatza et al., 2020)
Vd = 0.6;                       % L/kg, volume of distribution - (Patsalos, 2000)
V = Vd*mass;                    % L, volume of central compartment
ka = 0.616;                     % hr^-1, absorption rate constant - (Karatza et al., 2020)
kc = CL/V;                      % hr^-1, clearance rate constant
q = 0;                          % mg/hr, 0 because no infusion

p = [ka, kc, V, q]; % parameters to input to simulation

% Run simulation for base case: avg height, weight, & CrCL 10 doses
[T,Y,BalanceD,AUC,Ctrough] = LEV_sim(D0,p,DoseTimes);

% Therapeutic range lower and upper bounds (Karatza et al., 2020)
LB = zeros(size(T))+12; % mg/L
UB = zeros(size(T))+46; % mg/L

% Visualization for drug concentration in body and mass balance
figure;
subplot(2,1,1)
plot(T,Y(:,1),'linewidth',2)
hold on
scatter(DoseTimes(2:end),Ctrough,'o')
plot(T,LB,'--r','linewidth',1)
plot(T,UB,'--r','linewidth',1)
xlabel('Time (hrs)')
ylabel('Levetiracetam (mg/L)')
legend('Concentration Profile','Ctrough','Therapeutic Range - Lower','Therapeutic Range - Upper')

subplot(2,1,2)
plot(T,BalanceD,'linewidth',2)
title('Mass Balance')
xlabel('Time (hrs)')


%% Population Analysis - Creating the virtual population

% Parameters likely to vary with population:
% - CrCL
% - weight (volume of distribution varied by weight)

%%%%% Creating a virtual population of 100 patients %%%%%%%%%%


NumberOfSubjects = 100; % size of population

% Mean & SD values from Karatza et al., 2020 and their population of 8 subjects
meanwt = 79; % kg
sdwt = 15;
cutoffwt = 30; % Minimum weight; abritrarily set to a reasonable value

meanht = 170; % cm
sdht = 8;
cutoffht = 140; % Minimum height; abritrarily set to a reasonable value

meanCrCL = 128.73; % mL/min
sdCrCL = 26.33;

% Minimum accepted bmi - assuring the randomly generated population has
% reasonable height and weight combinations
cutoffbmi = 13;  % kg/m^2 - (CDC)

% GENERATE VIRTUAL POPULATION (using random numbers)

% Initiate Random Numbers (helpful to make results reproducible)
rng(0,'twister');

% Weight, height and CrCL values for each subject randomly generated
xtemp = sdwt.*randn(NumberOfSubjects,1)+meanwt;
htemp = sdht.*randn(NumberOfSubjects,1)+meanht;
Crtemp = sdCrCL.*randn(NumberOfSubjects,1)+meanCrCL;
% Finding the BMI value for each subject using the generated height & weight
bmitemp = xtemp./((htemp/100).^2);

% This gives us a first attempt; next we must identify weights, heights and bmi
% below the cutoff and replace them
a=length(xtemp(xtemp<=cutoffwt));      % count people below x kg
b=length(htemp(htemp<=cutoffht));      % count people below x cm
c=length(bmitemp(bmitemp<=cutoffbmi)); % count people below x bmi (kg/m^2)
cycle = 1;

while a>0 | b>0 | c>0
    
    if a>0 % if there are any below cutoff, replace them
        fprintf ('series %i, cycle %i, negs %i\n',i,1,a);
        xtemp(xtemp<=cutoffwt)=sdwt.*randn(a,1)+meanwt;        
        a=length(xtemp(xtemp<=cutoffwt)); 
    end

    if b>0 % if there are any below cutoff, replace them
        fprintf ('series %i, cycle %i, negs %i\n',i,1,b);
        htemp(htemp<=cutoffht)=sdht.*randn(b,1)+meanht;        
        b=length(htemp(htemp<=cutoffht)); 
    end
    
    % recheck the bmi based on updated height and weight vals
    bmitemp = xtemp./((htemp/100).^2);
    c=length(bmitemp(bmitemp<=cutoffbmi));

    if c>0 % if there are any below cutoff, replace them
        fprintf ('series %i, cycle %i, negs %i\n',i,1,c);
        xtemp(bmitemp<=cutoffbmi)=sdwt.*randn(c,1)+meanwt; 
        htemp(bmitemp<=cutoffbmi)=sdht.*randn(c,1)+meanht; 
        c=length(bmitemp(bmitemp<=cutoffbmi)); 
    end 
    
    % Updating a, b and c values after regenerating cutoff vals
    a=length(xtemp(xtemp<=cutoffwt));
    b=length(htemp(htemp<=cutoffht)); 
    c=length(bmitemp(bmitemp<=cutoffbmi));

    cycle = cycle + 1;
end % check again for below cutoff and loop if nec.


xdist = xtemp;     % This is the final weight distribution
hdist = htemp;     % This is the final height distribution
bmidist = bmitemp; % This is the final bmi distribution

% Output the means, SDs of the simulated population(s)
simmean = mean(xdist);
hmean = mean(hdist);
Crsimmean = mean(Crtemp);
simSD = std(xdist);
hSD = std(hdist);
CrsimSD = std(Crtemp);
bmean = mean(bmidist);
bSD = std(bmidist);
fprintf('Mass - Simulated population, input mean %4.1f, simulated mean %4.1f; input SD %4.1f, simulated SD %4.1f \n',meanwt,simmean,sdwt,simSD)
fprintf('CrCL - Simulated population, input mean %4.1f, simulated mean %4.1f; input SD %4.1f, simulated SD %4.1f \n',meanCrCL,Crsimmean,sdCrCL,CrsimSD)
fprintf('Height - Simulated population, input mean %4.1f, simulated mean %4.1f; input SD %4.1f, simulated SD %4.1f \n',meanht,hmean,sdht,hSD)
fprintf('BMI - Simulated population, simulated mean %4.1f; simulated SD %4.1f \n',bmean,bSD)


patientID = (1:NumberOfSubjects)';
Weight = xdist;
Height = hdist;
CrCLp = Crtemp;
BMI = bmidist;

%% Population Analysis - Assemble populations to analyze

%%%% Population A only varying CrCL (mL/min) %%%%
wt_a = Weight;
wt_a(:,1) = 79;  % kg - set everyone to population average
ht_a = Height;
ht_a(:,1) = 170; % cm - set everyone to population average

% Calculate bmi (sanity check), BSA, BSA adjusted CrCL, and corresponding
% doses for each patient.
bmi_a = wt_a./((ht_a/100).^2);   % kg/m^2
BSA_a = sqrt((ht_a.*wt_a)/3600); % m^2
CrCL_a = (CrCLp./BSA_a)*1.73;    % mL/min/1.73m^2
D0_a = zeros(size(CrCL_a));      % mg

% Select dose for patient based on CrCL (mL/min/1.73m^2) - via FDA label
for i = 1:NumberOfSubjects
    c = CrCL_a(i);
    if c >= 80 % Normal Group
        D0_a(i) = 1500; % mg
    elseif c >= 50  % Mild Impairment Group
        D0_a(i) = 1000; % mg
    elseif c >= 30  % Moderate Impairment Group
        D0_a(i) = 750; % mg
    else % c<30, Severe Impairment Group
        D0_a(i) = 500; % mg
    end
    
end

% Create matrix for Population A containing info required for simulation
pop_a = [patientID, wt_a, ht_a, CrCLp, bmi_a, CrCL_a, D0_a];



%%%% Population B varying weight & height (CrCL held constant) %%%%
CrCL_bm = CrCLp;
CrCL_bm(:,1) = 128.73; % mL/min - set everyone to population average

% Calculate bmi (sanity check), BSA, BSA adjusted CrCL, and corresponding
% doses for each patient.
BSA_b =  sqrt((Height.*Weight)/3600); % m^2
CrCL_b = (CrCL_bm./BSA_b)*1.73;       % mL/min/1.73m^2
D0_b = zeros(size(CrCL_a));           % mg

% Select dose for patient based on CrCL (mL/min/1.73m^2) - via FDA label
for i = 1:NumberOfSubjects
    c = CrCL_b(i);
    if c >= 80 % Normal Group
        D0_b(i) = 1500; % mg
    elseif c >= 50  % Mild Impairment Group
        D0_b(i) = 1000; % mg
    elseif c >= 30  % Moderate Impairment Group
        D0_b(i) = 750; % mg
    else % c<30, Severe Impairment Group
        D0_b(i) = 500; % mg
    end
    
end

% Create matrix for Population B containing info required for simulation
pop_b = [patientID, Weight, Height, CrCL_bm, BMI, CrCL_b, D0_b];



%%%% Population C varying height, weight and CrCL %%%%
% values pulled directly from randomnly generated population

% Calculate BSA, BSA adjusted CrCL, and corresponding doses for each patient.
BSA_c = sqrt((Height.*Weight)/3600); % m^2
CrCL_c = (CrCLp./BSA_c)*1.73;        % mL/min/1.73m^2
D0_c = zeros(size(CrCL_c));          % mg

% Select dose for patient based on CrCL (mL/min/1.73m^2) - via FDA label
for i = 1:NumberOfSubjects
    c = CrCL_c(i);
    if c >= 80 % Normal Group
        D0_c(i) = 1500; % mg
    elseif c >= 50 %|| c < 80 % Mild Impairment Group
        D0_c(i) = 1000; % mg
    elseif c >= 30 %|| c < 50 % Moderate Impairment Group
        D0_c(i) = 750; % mg
    else % c<30, Severe Impairment Group
        D0_c(i) = 500; % mg
    end
    
end

% Create matrix for Population C containing info required for simulation
pop_c = [patientID, Weight, Height, CrCLp, BMI, CrCL_c, D0_c];


%%%% Population C with IIV - varying height, weight and CrCL %%%%
% Karatza et al. found IIV terms for CL, V, and ka
IIV_mean = 0;
IIV_CLstd = 0.159; % (Karatza et al., 2020)
IIV_kastd = 0.327; % (Karatza et al., 2020)
IIV_Vstd = 0.274;  % (Karatza et al., 2020)
nCL = IIV_CLstd.*randn(NumberOfSubjects,1)+IIV_mean;
nka = IIV_kastd.*randn(NumberOfSubjects,1)+IIV_mean;
nV = IIV_Vstd.*randn(NumberOfSubjects,1)+IIV_mean;


%% Population Analysis - Run the simulations for the populations

% Parameters - Baseline
num_doses = 10;  % count
dose_freq = 12;  % hrs
first_doseT = 0; % hrs
last_doseT = num_doses*dose_freq; % hrs, simulation end time (time of last dose + dose freq)
DoseTimes = linspace(first_doseT,last_doseT,(num_doses+1)); % hrs, vector for dose and total simulation timing


mass = 79;                      % kg
height = 170;                   % cm
D0 = 1500;                      % mg
Vd = 0.6;                       % L/kg
BSA = sqrt((height*mass)/3600); % m^2
CrCL = 128.73;                  % mL/min
CrCL = (CrCL/BSA)*1.73;         % mL/min/1.73m^2
CL = 3.26*(CrCL/139)^0.795;     % L/hr
V = Vd*mass;                    % L
ka = 0.616;                     % hr^-1
kc = CL/V;                      % hr^-1
q = 0;                          % mg/hr, 0 because no infusion



%%% Simulating the populations

% Initialize Output Matrices
% Key metric #1 - AUC
AUC_a = zeros(NumberOfSubjects,1);
AUC_b = zeros(NumberOfSubjects,1);
AUC_c = zeros(NumberOfSubjects,1);
AUC_cIIV = zeros(NumberOfSubjects,1);
% Key metric #2 - Ctrough
Ct_a = zeros(num_doses,NumberOfSubjects);
Ct_b = zeros(num_doses,NumberOfSubjects);
Ct_c = zeros(num_doses,NumberOfSubjects);
Ct_cIIV = zeros(num_doses,NumberOfSubjects);
% Concentration profile for each pop
Y_a = zeros(length(Y),NumberOfSubjects);
Y_b = zeros(length(Y),NumberOfSubjects);
Y_c = zeros(length(Y),NumberOfSubjects);
Y_cIIV = zeros(length(Y),NumberOfSubjects);


for p = 1:NumberOfSubjects
    
    %%% Population A - Varied CrCL (mL/min) %%%%%
    mass_a = pop_a(p,2);              % kg
    CrCL_a = pop_a(p,6);              % mL/min/1.73m^2
    CL_a = 3.26*(CrCL_a/139)^0.795;   % L/hr
    V_a = Vd*mass_a;                  % L
    kc_a = CL_a/V_a;                  % hr^-1
    D0_a = pop_a(p,7);                % mg
    param_a = [ka, kc_a, V_a, q];     % parameters to input to simulation
    
    % Run simulation, save AUC, Ctroughs and Concentration values
    [~,Ya,~,AUC_a(p),Ct_a(:,p)] = LEV_sim(D0_a,param_a,DoseTimes);
    
    Y_a(:,p) = Ya(:,1);

    
    %%% Population B - Varied Weight (kg) & Height (cm) %%%%
    mass_b = pop_b(p,2);              % kg
    CrCL_b = pop_b(p,6);              % mL/min/1.73m^2
    CL_b = 3.26*(CrCL_b/139)^0.795;   % L/hr
    V_b = Vd*mass_b;                  % L
    kc_b = CL_b/V_b;                  % hr^-1
    D0_b = pop_b(p,7);                % mg
    param_b = [ka, kc_b, V_b, q];     % parameters to input to simulation
    
     % Run simulation, save AUC, Ctroughs and Concentration values
    [~,Yb,~,AUC_b(p),Ct_b(:,p)] = LEV_sim(D0_b,param_b,DoseTimes);
    
    Y_b(:,p) = Yb(:,1);
    
    %%% Population C - Varied CrCL, Weight, and Height %%%%
    mass_c = pop_c(p,2);              % kg
    CrCL_c = pop_c(p,6);              % mL/min/1.73m^2
    CL_c = 3.26*(CrCL_c/139)^0.795;   % L/hr
    V_c = Vd*mass_c;                  % L
    kc_c = CL_c/V_c;                  % hr^-1
    D0_c = pop_c(p,7);                % mg
    param_c = [ka, kc_c, V_c, q];     % parameters to input to simulation
    
     % Run simulation, save AUC, Ctroughs and Concentration values
    [~,Yc,~,AUC_c(p),Ct_c(:,p)] = LEV_sim(D0_c,param_c,DoseTimes);
    
    Y_c(:,p) = Yc(:,1);
    
    %%% Population C with IIV - Varied CrCL, Weight, and Height %%%%
    ka_cIIV = ka*exp(nka(p));                        % hr^-1
    CL_cIIV = (3.26*(CrCL_c/139)^0.795)*exp(nCL(p)); % L/hr
    V_cIIV = (Vd*mass_c)*exp(nV(p));                 % L
    kc_cIIV = CL_cIIV/V_c;                           % hr^-1
    param_cIIV = [ka_cIIV, kc_cIIV, V_cIIV, q];      % parameters to input to simulation
    
     % Run simulation, save AUC, Ctroughs and Concentration values
    [~,YcIIV,~,AUC_cIIV(p),Ct_cIIV(:,p)] = LEV_sim(D0_c,param_cIIV,DoseTimes);
    
    Y_cIIV(:,p) = YcIIV(:,1);
    
    
end

%%% Saving the appropriate variables for visualization in R %%%%

% Save AUC values for each simulated population
save AUCpops.mat AUC_a AUC_b AUC_c AUC_cIIV

% Save Ctrough values after the 10th dose for each pop
CT_a = Ct_a(end,:)';
CT_b = Ct_b(end,:)';
CT_c = Ct_c(end,:)';
CT_cIIV = Ct_cIIV(end,:)';
save Ct10pops.mat CT_a CT_b CT_c CT_cIIV

% Save AUC and Ctrough values by dosage
pop_a = [pop_a, AUC_a, CT_a];
pop_b = [pop_b, AUC_b, CT_b];
pop_cIIV = [pop_c, AUC_cIIV, CT_cIIV];
pop_c = [pop_c, AUC_c, CT_c];

dose_metrics = [pop_a(:,6:9);pop_b(:,6:9);pop_c(:,6:9);pop_cIIV(:,6:9)];
save dose_metrics.mat dose_metrics



% Save the concentration profiles for the 4 populations
save conc_a.mat T Y_a
save conc_b.mat T Y_b
save conc_c.mat T Y_c
save conc_cIIV.mat T Y_cIIV


save pop_a.mat pop_a
save pop_b.mat pop_b
save pop_c.mat pop_c
save pop_cIIV.mat pop_cIIV

% Finding correlation between

%% Reactreating Figure 4 - Karatza paper

CrCL_k = [110, 90, 60, 45, 30, 15]; % mL/min
D0_k = [1500, 1000, 750, 500];      % mg


% Parameters - Baseline
num_doses = 10;  % count
dose_freq = 12;  % hrs
first_doseT = 0; % hrs
last_doseT = num_doses*dose_freq; % hrs, simulation end time (time of last dose + dose freq)
DoseTimes = linspace(first_doseT,last_doseT,(num_doses+1)); % hrs, vector for dose and total simulation timing


mass = 79;                      % kg
height = 170;                   % cm
Vd = 0.6;                       % L/kg
BSA = sqrt((height*mass)/3600); % m^2
V_k = Vd*mass;                  % L
ka = 0.616;                     % hr^-1
q = 0;                          % mg/hr, 0 because no infusion

% Initialize output matrices
Y_k = zeros(length(Y),(length(D0_k)*length(CrCL_k)));

% Simulations
col = 1;
for d = 1:length(D0_k)
    
    % Set the dose
    D0k = D0_k(d);
    
    for j = 1:length(CrCL_k)
        
        % set the CrCL value
        CrCL = CrCL_k(j);             % mL/min
        CrCL = (CrCL/BSA)*1.73;       % mL/min/1.73m^2
        CL_k = 3.26*(CrCL/139)^0.795; % L/hr
        kc_k = CL_k/V_k;              % hr^-1
        param_k = [ka, kc_k, V_k, q]; % parameters to input to simulation
        
        [Tk,Yk,~,~,~] = LEV_sim(D0k,param_k,DoseTimes);
        
        Y_k(:,col) = Yk(:,1);
        col = col + 1;
        
    end

    
end

Y_1500 = Y_k(:,1:6);
Y_1000 = Y_k(:,7:12);
Y_750 = Y_k(:,13:18);
Y_500 = Y_k(:,19:24);

% Save data for visualization in R
save Cprof1500.mat Tk Y_1500
save Cprof1000.mat Tk Y_1000
save Cprof750.mat Tk Y_750
save Cprof500.mat Tk Y_500

%%%% Initial visualization, final visualization in R %%%%

figure;
plot(T,Y_1500(:,1),T,Y_1500(:,2),T,Y_1500(:,3),T,Y_1500(:,4),T,Y_1500(:,5),T,Y_1500(:,6),'linewidth',2)
hold on
%scatter(DoseTimes(2:end),Ctrough,'o')
plot(T,LB,'--r','linewidth',1)
plot(T,UB,'--r','linewidth',1)
xlabel('Time (hrs)')
ylabel('Levetiracetam (mg/L)')
legend('110 mL/min','90 mL/min','60 mL/min','45 mL/mi','30 mL/min','15 mL/min','Therapeutic Range - Lower','Therapeutic Range - Upper')
title('1500mg Twice Daily')

figure;
plot(T,Y_1000(:,1),T,Y_1000(:,2),T,Y_1000(:,3),T,Y_1000(:,4),T,Y_1000(:,5),T,Y_1000(:,6),'linewidth',2)
hold on
%scatter(DoseTimes(2:end),Ctrough,'o')
plot(T,LB,'--r','linewidth',1)
plot(T,UB,'--r','linewidth',1)
xlabel('Time (hrs)')
ylabel('Levetiracetam (mg/L)')
legend('110 mL/min','90 mL/min','60 mL/min','45 mL/mi','30 mL/min','15 mL/min','Therapeutic Range - Lower','Therapeutic Range - Upper')
title('1000mg Twice Daily')

figure;
plot(T,Y_750(:,1),T,Y_750(:,2),T,Y_750(:,3),T,Y_750(:,4),T,Y_750(:,5),T,Y_750(:,6),'linewidth',2)
hold on
%scatter(DoseTimes(2:end),Ctrough,'o')
plot(T,LB,'--r','linewidth',1)
plot(T,UB,'--r','linewidth',1)
xlabel('Time (hrs)')
ylabel('Levetiracetam (mg/L)')
legend('110 mL/min','90 mL/min','60 mL/min','45 mL/mi','30 mL/min','15 mL/min','Therapeutic Range - Lower','Therapeutic Range - Upper')
title('750mg Twice Daily')

figure;
plot(T,Y_500(:,1),T,Y_500(:,2),T,Y_500(:,3),T,Y_500(:,4),T,Y_500(:,5),T,Y_500(:,6),'linewidth',2)
hold on
%scatter(DoseTimes(2:end),Ctrough,'o')
plot(T,LB,'--r','linewidth',1)
plot(T,UB,'--r','linewidth',1)
xlabel('Time (hrs)')
ylabel('Levetiracetam (mg/L)')
legend('110 mL/min','90 mL/min','60 mL/min','45 mL/mi','30 mL/min','15 mL/min','Therapeutic Range - Lower','Therapeutic Range - Upper')
title('500mg Twice Daily')
