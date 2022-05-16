close all;
clear all;

%% KEY SIMULATION PARAMETERS
flag=1;
% 1 = local univariate
% 2 = local bivariate
% 3 = local parameter vs global dose
% 4 = local parameter vs global concentration
% 5 = global sensitivity

TimeLen=12;

OutputVar = 1;

ParamDelta = 0.05; % test sensitivity to a 5% change

%% PARAMETERS 

% Model parameters for array p
mass = 80;   % kg
D0 = 1500;   % mg
Vd = 0.6;    % L/kg
CL = 3.26;   % L/hr
V = Vd*mass; % L
ka = 0.616;  % hr^-1
kc = CL/V;   % hr^-1
q = 0;       % mg/hr, 0 because no infusion
num_doses = 10;  % count
dose_freq = 12;  % hrs
first_doseT = 0; % hrs
last_doseT = num_doses*dose_freq; % hrs, simulation end time (time of last dose + dose freq)
%TimeLen = num_doses*dose_freq;        % hrs, total time for simulation
DoseTimes = linspace(first_doseT,last_doseT,(num_doses+1)); % hrs, vector for dose and total simulation timing


mapr = [linspace(0.33,1,101) linspace(.99,0,100)]; 
mapg = [linspace(0,1,101) linspace(.99,.33,100)];
mapb = [linspace(0.33,1,101) linspace(0.99,0,100)];
map =[mapr' mapg' mapb'];


%% RUN BASE CASE

p0 = [ka, kc, V, q];
p0labels = {'ka' 'kc' 'V'}';
outputvar = 1;
[t0,y0, balance0, auc0] = LEV_sim(D0,p0, DoseTimes);
y0 = y0(:,outputvar);
y=y0;
size(y)
t=t0;

%% RUN SENSITIVITY SIMULATIONS

switch flag
    case 1
        
%% 1. SENSITIVITY - LOCAL UNIVARIATE
%========================
% OUTPUT: 1 = concentration of drug 
% INPUT: all parameters

for i=1:(length(p0)-1)
    p=p0;
    p(i)=p0(i)*(1.0+ParamDelta);
    [t1, y1, balance1, auc(i)] = LEV_sim(D0,p, DoseTimes);
    t1 = t1(:,outputvar);
    y1 = y1(:,outputvar);
    %auccur(i) = auc(:,1);
    t = [t t1];
    y = [y y1];
    Sens(i) = ((auc(i)-auc0)/auc0)/((p(i)-p0(i))/p0(i)); % relative sensitivity vs parameter
    SensB(i) = ((auc(i)-auc0))/((p(i)-p0(i))); % absolute sensitivity vs parameter
    SensAbs(i) = ((auc(i)-auc0)/auc0); % relative change (not normalized to parameter)
    SensT = ((y1-y0)./y0)/((p(i)-p0(i))/p0(i)); % time-dependent relative sensitivity
    [maxS,indT] = max(abs(SensT)); 
    SensTmax(i) = t1(indT); % time to max relative sensitivity

    if i==1
        SensTarr = SensT;

    else
        SensTarr = [SensTarr SensT];
    end
end


t = t;
SensTarr = SensTarr;
Sens = Sens;
SensB = SensB;
SensTmax = SensTmax;

save 1_1.mat t SensTarr
save 1_2.mat Sens
save 1_3.mat SensB
save 1_4.mat SensTmax
save 1_5.mat t y




    case 2
%% 2. SENSITIVITY - LOCAL BIVARIATE
%========================
% OUTPUT: AUC of concentration
%         (if OutputVar = 1, then concentration in blood, etc)
% INPUT: all parameters x all parameters

for i=1:(length(p0)-1)
    p=p0;
    p(i)=p0(i)*(1.0+ParamDelta);
    [t1, y1, balance1, auc(i)] = LEV_sim(D0,p, DoseTimes);
    t1 = t1;
    y1 = y1(:,outputvar);
    t = [t t1];
    y = [y y1];
    Sens(i) = ((auc(i)-auc0)/auc0)/((p(i)-p0(i))/p0(i));
    SensB(i) = ((auc(i)-auc0))/((p(i)-p0(i)));
    SensAbs(i) = ((auc(i)-auc0)/auc0);
    SensT = ((y1-y0)./y0)/((p(i)-p0(i))/p0(i));
end

for i=1:(length(p0)-1)
for j=1:(length(p0)-1)
    p=p0;
    p(i)=p0(i)*(1.0+ParamDelta);
    p(j)=p0(j)*(1.0+ParamDelta);
    [t2, y2, balance1, auc2(i,j)] = LEV_sim(D0,p, DoseTimes);
    t2=t2(:,1);
    y2=y2(:,1);
    Sens2(i,j) = ((auc2(i,j)-auc(j))/auc(j))/((p(i)-p0(i))/p0(i));
    Sens2Abs(i,j) = ((auc2(i,j)-auc0)/auc0);
    y0=y(:,j+1);
    Sens2T = ((y1-y0)./y0)/((p(i)-p0(i))/p0(i));
    [maxS,indT] = max(abs(Sens2T));
    Sens2Tmax(i,j) = t1(indT);
end
end

for i=1:(length(p0)-1)
for j=1:(length(p0)-1)
	Sens2Syn(i,j) = Sens2Abs(i,j)/(SensAbs(i)*SensAbs(j));
end
end

save 2_1.mat Sens2
save 2_2.mat Sens2Tmax 
save 2_3.mat Sens2Syn 


    case 3
        
%% 3. SENSITIVITY - LOCAL BIVARIATE vs Dose
%========================
% OUTPUT: AUC of concentration
%         (if OutputVar = 1, then concentration, etc)
% INPUT: all parameters and range of drug doses

tic 
p1 = p0;
Doses = round(linspace(100,2500,10));
for i=1:length(Doses)
    D0 = Doses(i);
    [t0, y0, balance0, auc0] = LEV_sim(D0,p1, DoseTimes);
    t0 = t0;
    y0 = y0(:,outputvar);

for j=1:(length(p0)-1)
    p=p1;
    p(j)=p1(j)*(1.0+ParamDelta);
    [t3, y3, balance3, auc3(i,j)]= LEV_sim(D0,p, DoseTimes);
    t3 = t;
    y3 = y3(:,outputvar);
    Sens3(i,j) = ((auc3(i,j)-auc0)/auc0)/((p(j)-p1(j))/p1(j));
    Sens3T = ((y3-y0)./y0)/((p(j)-p1(j))/p1(j));
    [maxS,indT] = max(abs(Sens3T));
    Sens3Tmax(i,j) = t3(indT);
end
end

toc
%Doselabels = num2str(Doses,'%6.0f');


save 3_1.mat Sens3
save 3_2.mat Sens3Tmax

    case 4

%% 4. SENSITIVITY - LOCAL BIVARIATE vs Plasma Protein level
%========================
% OUTPUT: AUC of concentration
%         (if OutputVar = 1, then concentration, etc)
% INPUT: all parameters and range of drug doses

tic
p1 = p0;
Doses = round(linspace(100,2500,10));
for i=1:length(Doses)
    % BASE CASE
    p1(length(p1)-1) = Doses(i);
    [t0, y0, balance0, auc0] = LEV_sim(D0,p1, DoseTimes);
    t0=t0;
    y0=y0(:,1);

for j=1:length(p0)
    p=p1;
    p(j)=p1(j)*(1.0+ParamDelta);
    [t1, y1, balance, auc(i,j), ctrough] = LEV_sim(D0,p, DoseTimes);
    ctrough(i,j)=ctrough(end);
    t1 = t1(:);
    y1 = y1(:,outputvar);
    Sens3(i,j) = ((auc(i,j)-auc0)/auc0)/((p(j)-p1(j))/p1(j));
    Sens3T = ((y1-y0)./y0)/((p(j)-p1(j))/p1(j));
    [maxS,indT] = max(abs(Sens3T));
    Sens3Tmax(i,j) = t1(indT);
end
end
toc

%Doselabels = num2str(Doses,'%6.0f');
p0labels = {'1','2','3','4'};


save 4_1.mat Sens3
save 4_2.mat Sens3Tmax
save 4_3.mat ctrough



case 5
% 4. global sensitivity using values from case 2 %========================
    x=0;
    y0 = [20 0 0]'; % this is in mg/kg 
    TimeLen = 168;
    Yfinal3 = [0 0 0];
    Tfinal3 = [0];
    q=0;
    D1 = 20;
    D2 = 5;
    kA = 0;
    ydosage = [5 0 0]'; % this is in mg/kg 
    Balfinal = [0];
    Vd = 0.85;
    Doses = round(linspace(100,2500,20));
    for j=1:length(Doses)
            mass = linspace(50,150,20);
            V = Vd*mass;
            for i=1:length(mass)
                kc = CL/V(i);   % hr^-1
                p = [ka, kc, V(i), q];
                [t1, y1, balance, auc(i,j), ctrough] = LEV_sim(Doses(j),p, DoseTimes);
                ctrough(i,j)=ctrough(end);
                t1 = t1(:);
                y1 = y1(:,outputvar);
            end
        end
    Vfinal = V';
    DosesFinal = Doses';
    Final = [Vfinal DosesFinal];

    save globalAUC.mat auc
    save globalParam.mat Vfinal DosesFinal
    save globalctrough.mat ctrough


end