% EN.580.640 - Final Project
% Sophia Nehs, Valentina Dsouza, Caroline Ghio, Christianne Chua, Shruthi Bare
% Missed dose analysis population variability driver

close all; clear all


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
[T,Y,BalanceD,AUC,Ctrough,peak] = LEV_sim_peak(D0,p,DoseTimes);

% Therapeutic range lower and upper bounds (Karatza et al., 2020)
LB = zeros(size(T))+12; % mg/L
UB = zeros(size(T))+46; % mg/L



%% Population Analysis

% parameters likely to vary with population
% - CrCL
% - weight (volume of distribution varied by weight)

% Population generated vayring both CrCL and weight (which affects volume
% of distribution

NumberOfSubjects = 5;
% Size of the populations - by making this a parameter, we can easily test 
%   the code with small numbers that run quickly, and then run it with 
%   larger populations once we're satisfied that it's working.

meanwt = 79;
sdwt = 15;
cutoffwt = 30; % Minumum weight; set to 0 to only remove nonpositives

meanht = 170;
sdht = 8;
cutoffht = 140; % Minumum height; set to 0 to only remove nonpositives

cutoffbmi = 13; % Minimum accepted bmi

meanCrCL = 128.73;
sdCrCL = 26.33;

% GENERATE VIRTUAL POPULATION (using random numbers)

% Initiate Random Numbers (helpful to make results reproducible)
rng(0,'twister');

xtemp = sdwt.*randn(NumberOfSubjects,1)+meanwt;
htemp = sdht.*randn(NumberOfSubjects,1)+meanht;
Crtemp = sdCrCL.*randn(NumberOfSubjects,1)+meanCrCL;
bmitemp = xtemp./((htemp/100).^2);
% This gives us a first attempt; next we must identify weights below 
% the cutoff and replace them
a=length(xtemp(xtemp<=cutoffwt)); % count people below x lb
b=length(htemp(htemp<=cutoffht)); % count people below x cm
c=length(bmitemp(bmitemp<=cutoffbmi));
cycle = 1;

while a>0 | b>0 | c>0
    
    if a>0 % if there are any nonpositives, replace them
        %fprintf ('series %i, cycle %i, negs %i\n',i,1,a);
        xtemp(xtemp<=cutoffwt)=sdwt.*randn(a,1)+meanwt;        
        a=length(xtemp(xtemp<=cutoffwt)); 
    end

    if b>0 % if there are any nonpositives, replace them
        %fprintf ('series %i, cycle %i, negs %i\n',i,1,b);
        htemp(htemp<=cutoffht)=sdht.*randn(b,1)+meanht;        
        b=length(htemp(htemp<=cutoffht)); 
    end
    
    bmitemp = xtemp./((htemp/100).^2);
    c=length(bmitemp(bmitemp<=cutoffbmi));

    if c>0 % if there are any nonpositives, replace them
        fprintf ('series %i, cycle %i, negs %i\n',i,1,c);
        xtemp(bmitemp<=cutoffbmi)=sdwt.*randn(c,1)+meanwt; 
        htemp(bmitemp<=cutoffbmi)=sdht.*randn(c,1)+meanht; 
        c=length(bmitemp(bmitemp<=cutoffbmi)); 
        %cycle = cycle + 1;
    end % check again for nonpositives and loop if nec.
    
    a=length(xtemp(xtemp<=cutoffwt));
    b=length(htemp(htemp<=cutoffht)); 
    c=length(bmitemp(bmitemp<=cutoffbmi));

    cycle = cycle + 1;
    %end % check again for nonpositives and loop if nec.
end


xdist = xtemp; % This is the final weight distribution
hdist = htemp; % This is the final height distribution
bmidist = bmitemp;

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

% population C varying weight and CrCL
BSA_c = sqrt((Height.*Weight)/3600); % m^2
CrCL_c = (CrCLp./BSA_c)*1.73; % mL/min/1.73m^2

D0_c = zeros(size(CrCL_c)); % mg
% Select dose for patient based on CrCL (mL/min/2.73m^2)

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

pop_c = [patientID, Weight, Height, CrCLp, BMI, CrCL_c, D0_c];



%%

% Parameters - Baseline
num_doses = 10;  % count
dose_freq = 12;  % hrs
first_doseT = 0; % hrs
last_doseT = num_doses*dose_freq; % hrs, simulation end time (time of last dose + dose freq)
DoseTimes = linspace(first_doseT,last_doseT,(num_doses+1)); % hrs, vector for dose and total simulation timing
missed_DoseTimes = linspace(first_doseT,last_doseT,(num_doses+1));
diff = missed_DoseTimes(2)-missed_DoseTimes(1);
missed_DoseT = zeros(10,1);
%setting up misseddosetimes
for i = 1:length(missed_DoseTimes)
    if i==6
        missed_DoseTimes1(i) = missed_DoseTimes(i) + diff/5;
        missed_DoseTimes2(i) = missed_DoseTimes(i) + 2*diff/5;
        missed_DoseTimes3(i) = missed_DoseTimes(i) + 3*diff/5;
        missed_DoseTimes4(i) = missed_DoseTimes(i) + 4*diff/5;
        missed_DoseT(i) = missed_DoseTimes(i);
    else
        missed_DoseTimes1(i) = missed_DoseTimes(i);
        missed_DoseTimes2(i) = missed_DoseTimes(i);
        missed_DoseTimes3(i) = missed_DoseTimes(i);
        missed_DoseTimes4(i) = missed_DoseTimes(i);
        missed_DoseT(i) = missed_DoseTimes(i);
    end
end
missed_DoseT(6)=[];
%

mass = 79;    % kg
height = 170; % cm
D0 = 1500;    % mg
Vd = 0.6;     % L/kg
BSA = sqrt((height*mass)/3600); % m^2
CrCL = 128.73;                  % mL/min
CrCL = (CrCL/BSA)*1.73;         % mL/min/1.73m^2
CL = 3.26*(CrCL/139)^0.795;     % L/hr
V = Vd*mass; % L
ka = 0.616;  % hr^-1
kc = CL/V;   % hr^-1
q = 0;       % mg/hr, 0 because no infusion



%%% Simulating the populations

% Initialize Output Matrices
% Key metric #1 - AUC
AUC_c = zeros(NumberOfSubjects,1);
AUC_c1 = zeros(NumberOfSubjects,1);
AUC_c2 = zeros(NumberOfSubjects,1);
AUC_c3 = zeros(NumberOfSubjects,1);
AUC_c4 = zeros(NumberOfSubjects,1);
AUC_cm = zeros(NumberOfSubjects,1);
% Key metric #2 - Ctrough
Ct_c = zeros(num_doses,NumberOfSubjects);
Ct_c1 = zeros(num_doses,NumberOfSubjects);
Ct_c2 = zeros(num_doses,NumberOfSubjects);
Ct_c3 = zeros(num_doses,NumberOfSubjects);
Ct_c4 = zeros(num_doses,NumberOfSubjects);
Ct_cm = zeros(num_doses,NumberOfSubjects);
%peak
Pk_c = zeros(num_doses,NumberOfSubjects);
Pk_c1 = zeros(num_doses,NumberOfSubjects);
Pk_c2 = zeros(num_doses,NumberOfSubjects);
Pk_c3 = zeros(num_doses,NumberOfSubjects);
Pk_c4 = zeros(num_doses,NumberOfSubjects);
Pk_cm = zeros(num_doses,NumberOfSubjects);
% Concentration profile info
Y_c = zeros(length(Y),NumberOfSubjects);
Y_c1 = zeros(length(Y),NumberOfSubjects);
Y_c2 = zeros(length(Y),NumberOfSubjects);
Y_c3 = zeros(length(Y),NumberOfSubjects);
Y_c4 = zeros(length(Y),NumberOfSubjects);
Y_cm = zeros(length(Y),NumberOfSubjects);
for p = 1:NumberOfSubjects

    %%% Population C - Varied CrCL, Weight, and Height %%%%
    mass_c = pop_c(p,2);              % kg
    %height = pop_c(p,3);             % cm
    CrCL_c = pop_c(p,6);              % mL/min/1.73m^2
    CL_c = 3.26*(CrCL_c/139)^0.795;   % L/hr
    V_c = Vd*mass_c;                  % L
    kc_c = CL_c/V_c;                  % hr^-1
    D0_c = pop_c(p,7);                % mg
    param_c = [ka, kc_c, V_c, q];     % parameters to input to simulation
    
     % Run simulation, save AUC, Ctroughs and Concentration values
    [Tc,Yc,~,AUC_c(p),Ct_c(:,p),Pk_c(:,p)] = LEV_sim_peak(D0_c,param_c,DoseTimes);
    [Tc1,Yc1,~,AUC_c1(p),Ct_c1(:,p),Pk_c1(:,p)] = LEV_sim_peak(D0_c,param_c,missed_DoseTimes1);
    [Tc2,Yc2,~,AUC_c2(p),Ct_c2(:,p),Pk_c2(:,p)] = LEV_sim_peak(D0_c,param_c,missed_DoseTimes2);
    [Tc3,Yc3,~,AUC_c3(p),Ct_c3(:,p),Pk_c3(:,p)] = LEV_sim_peak(D0_c,param_c,missed_DoseTimes3);
    [Tc4,Yc4,~,AUC_c4(p),Ct_c4(:,p),Pk_c4(:,p)] = LEV_sim_peak(D0_c,param_c,missed_DoseTimes4);
    [Tcm,Ycm,~,AUC_cm(p),Ct_c5(:,p),Pk_c5(:,p)] = LEV_sim_peak(D0_c,param_c,missed_DoseT);
    Y_c(:,p) = Yc(:,1);
    Y_c1(:,p) = Yc1(:,1);
    Y_c2(:,p) = Yc2(:,1);
    Y_c3(:,p) = Yc3(:,1);
    Y_c4(:,p) = Yc4(:,1);
    Y_cm(:,p) = Ycm(:,1);
end
%% setting up the times for peak values for missed doses
times = zeros(length(Pk_c),5);
times1 = zeros(length(Pk_c1),5);
times2 = zeros(length(Pk_c2),1);
times3 = zeros(length(Pk_c3),1);
times4 = zeros(length(Pk_c4),1);
timesm = zeros(length(Pk_c5),1);
for m=1:5
    for i=1:length(Y_c)
        for j=1:10
            for k=1:9
            if Y_c(i,m) == Pk_c(j,m)
                times(j,m)= Tc(i);
            end
            if Y_c1(i,m) == Pk_c1(j,m)
                times1(j,m)= Tc1(i);
            end
            if Y_c2(i,m) == Pk_c2(j,m)
                times2(j,m)= Tc2(i);
            end
            if Y_c3(i,m) == Pk_c3(j,m)
                times3(j,m)= Tc3(i);
            end
            if Y_c4(i,m) == Pk_c4(j,m)
                times4(j,m)= Tc4(i);
            end
            if Y_cm(i,m) == Pk_c5(k,m)
                timesm(k,m)= Tcm(i);
            end
            end
        end
    end

end
%%Only plotted for normal dose and missed by m/5, can plot more in R due to
%%more color options available there
figure;
plot(Tc,Y_c(:,1),'r',Tc,Y_c(:,2),'cyan',Tc,Y_c(:,3),'g',Tc,Y_c(:,4),'y',Tc,Y_c(:,5),'m','linewidth',2)
hold on
plot(Tc1,Y_c1(:,1),'--r',Tc1,Y_c1(:,2),'--b',Tc1,Y_c1(:,3),'--g',Tc1,Y_c1(:,4),'--y',Tc1,Y_c1(:,5),'--m','linewidth',2)
scatter(DoseTimes(2:end),Ct_c(:,1),'r*')
scatter(DoseTimes(2:end),Ct_c(:,2),'b*')
scatter(DoseTimes(2:end),Ct_c(:,3),'g*')
scatter(DoseTimes(2:end),Ct_c(:,4),'y*')
scatter(DoseTimes(2:end),Ct_c(:,5),'m*')
scatter(times(1:end,1),Pk_c(:,1),'r*')
scatter(times(1:end,2),Pk_c(:,2),'b*')
scatter(times(1:end,3),Pk_c(:,3),'g*')
scatter(times(1:end,4),Pk_c(:,4),'y*')
scatter(times(1:end,5),Pk_c(:,5),'m*')
scatter(missed_DoseTimes1(2:end),Ct_c1(:,1),'r*')
scatter(missed_DoseTimes1(2:end),Ct_c1(:,2),'b*')
scatter(missed_DoseTimes1(2:end),Ct_c1(:,3),'g*')
scatter(missed_DoseTimes1(2:end),Ct_c1(:,4),'y*')
scatter(missed_DoseTimes1(2:end),Ct_c1(:,5),'m*')
scatter(times1(1:end,1),Pk_c1(:,1),'r*')
scatter(times1(1:end,2),Pk_c1(:,2),'b*')
scatter(times1(1:end,3),Pk_c1(:,3),'g*')
scatter(times1(1:end,4),Pk_c1(:,4),'y*')
scatter(times1(1:end,5),Pk_c1(:,5),'m*')

plot(T,LB,'--r','linewidth',1)
plot(T,UB,'--r','linewidth',1)
xlabel('Time (hrs)')
ylabel('Levetiracetam (mg/L)')
%legend('110 mL/min','90 mL/min','60 mL/min','45 mL/mi','30 mL/min','15 mL/min','Therapeutic Range - Lower','Therapeutic Range - Upper')
title('Missed Dose POpk')

figure;
auc1 = [AUC_c(1) AUC_c1(1) AUC_c2(1) AUC_c3(1) AUC_c4(1) AUC_cm(1)];
auc2 = [AUC_c(2) AUC_c1(2) AUC_c2(2) AUC_c3(2) AUC_c4(2) AUC_cm(2)];
auc3 = [AUC_c(3) AUC_c1(3) AUC_c2(3) AUC_c3(3) AUC_c4(3) AUC_cm(3)];
auc4 = [AUC_c(4) AUC_c1(4) AUC_c2(4) AUC_c3(4) AUC_c4(4) AUC_cm(4)];
auc5 = [AUC_c(5) AUC_c1(5) AUC_c2(5) AUC_c3(5) AUC_c4(5) AUC_cm(5)];
change = zeros(5,6);
auc = [AUC_c AUC_c1 AUC_c2 AUC_c3 AUC_c4 AUC_cm];
for i=1:5
    for j=2:6
        change(i,j) = (auc(i,1)-auc(i,j))/auc(i,1);
    end
end
change = change(:,2:end)*100;

save change_5.mat change

plot(auc1,'r');
hold on
plot(auc2,'b');
plot(auc3,'g');
plot(auc4,'y');
plot(auc5,'m');
scatter([1 2 3 4 5 6],auc1,'r');
scatter([1 2 3 4 5 6],auc2,'b');
scatter([1 2 3 4 5 6],auc3,'g');
scatter([1 2 3 4 5 6],auc4,'y');
scatter([1 2 3 4 5 6],auc5,'m');
set(gca, 'XTick', [1 2 3 4 5 6])
set(gca, 'XTickLabel', {'AUC when dose is on time' 'AUC when dose late by m/5' 'AUC when dose late by 2m/5' 'AUC when dose late by 3m/5' 'AUC when dose late by 4m/5' 'AUC when dose is missed once'})
title('AUC')
xlabel('Different dose conditions')
y = Y_c(:,1);
y1 = Y_c1(:,1);
y2 = Y_c2(:,1);
y3 = Y_c3(:,1);
y4 = Y_c4(:,1);
ym = Y_cm(:,1);
% peaks = [peak peak1 peak2 peak3 peak4 peakm];
DoseTimes = DoseTimes(2:end)';
missed_DoseT = missed_DoseT(2:end);
missed_DoseT = [missed_DoseT missed_DoseT missed_DoseT missed_DoseT missed_DoseT];
p_times = [times1 times2 times3 times4];
ct_times = [missed_DoseTimes1(2:end)' missed_DoseTimes2(2:end)' missed_DoseTimes3(2:end)' missed_DoseTimes4(2:end)'];

save AUc_values_pop_old.mat AUC_c AUC_c1 AUC_c2 AUC_c3 AUC_c4 AUC_cm
save change_5.mat change
save c_p_values_delayed_pop_old.mat Ct_c1 Ct_c2 Ct_c3 Ct_c4 Pk_c1 Pk_c2 Pk_c3 Pk_c4 ct_times p_times
save c_p_values_missed_pop_old.mat Ct_c5 Pk_c5 missed_DoseT timesm
save conc_normal_pop_old.mat Tc Y_c 
save conc_delayed_pop_old.mat Tc1 Tc2 Tc3 Tc4 Y_c1 Y_c2 Y_c3 Y_c4  
save c_p_normal_pop_old.mat Ct_c Pk_c DoseTimes times
save conc_missed_pop_old.mat Tcm Y_cm
crcl = pop_c(:,6);
save crcl_values.mat crcl