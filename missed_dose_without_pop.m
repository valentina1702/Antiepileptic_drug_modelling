% EN.580.640 - Final Project
% Sophia Nehs, Valentina Dsouza, Caroline Ghio, Christianne Chua, Shruthi Bare
% Missed dose analysis w/o population driver

close all; clear all


num_doses = 10;  % count
dose_freq = 12;  % hrs

first_doseT = 0; % hrs
last_doseT = num_doses*dose_freq; % hrs, simulation end time (time of last dose + dose freq)

%TimeLen = num_doses*dose_freq;        % hrs, total time for simulation

DoseTimes = linspace(first_doseT,last_doseT,(num_doses+1)); % hrs, vector for dose and total simulation timing
missed_DoseTimes = linspace(first_doseT,last_doseT,(num_doses+1));
diff = missed_DoseTimes(2)-missed_DoseTimes(1);
missed_DoseT = zeros(10,1);
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
case_no = 1;
%please edit this value to run different case_numbers
% Use case_no = 1, to have result for 500 mg
% Use case_no = 2, to have result for 1500 mg
% Use case_no = 3, to have result for 1000 mg
% Use case_no = 4, to have result for 750 mg
switch case_no
    case 1
        D0 = 500;   % mg
        [T1,Y1,BalanceD1,AUC1,Ctrough1,peak1] = LEV_sim_peak(D0,p,missed_DoseTimes1);
        [T2,Y2,BalanceD2,AUC2,Ctrough2,peak2] = LEV_sim_peak(D0,p,missed_DoseTimes2);
        [T3,Y3,BalanceD3,AUC3,Ctrough3,peak3] = LEV_sim_peak(D0,p,missed_DoseTimes3);
        [T4,Y4,BalanceD4,AUC4,Ctrough4,peak4] = LEV_sim_peak(D0,p,missed_DoseTimes4);
        [T,Y,BalanceD,AUC,Ctrough,peak] = LEV_sim_peak(D0,p,DoseTimes);
        [Tm,Ym,BalanceDm,AUCm,Ctroughm,peakm] = LEV_sim_peak(D0,p,missed_DoseT);

        times = zeros(length(peak),1);
        times1 = zeros(length(peak1),1);
        times2 = zeros(length(peak2),1);
        times3 = zeros(length(peak3),1);
        times4 = zeros(length(peak4),1);
        timesm = zeros(length(peakm),1);
        for i=1:length(Y)
            for j=1:10
                for k=1:9
                if Y(i) == peak(j,1)
                    times(j)= T(i);
                end
                if Y1(i) == peak1(j,1)
                    times1(j)= T1(i);
                end
                if Y2(i) == peak2(j,1)
                    times2(j)= T2(i);
                end
                if Y3(i) == peak3(j,1)
                    times3(j)= T3(i);
                end
                if Y4(i) == peak4(j,1)
                    times4(j)= T4(i);
                end
                if Ym(i) == peakm(k,1)
                    timesm(k)= Tm(i);
                end
                end
            end
        end




        figure;
        plot(T,Y(:,1),'k','linewidth',2)
        hold on
        plot(T1,Y1(:,1),'r','linewidth',2)
        plot(T2,Y2(:,1),'b','linewidth',2)
        plot(T3,Y3(:,1),'g','linewidth',2)
        plot(T4,Y4(:,1),'y','linewidth',2)
        plot(Tm,Ym(:,1),'c','linewidth',3)
        scatter(DoseTimes(2:end),Ctrough,'k*')
        scatter(missed_DoseTimes1(2:end),Ctrough1,'r*')
        scatter(missed_DoseTimes2(2:end),Ctrough2,'b*')
        scatter(missed_DoseTimes3(2:end),Ctrough3,'g*')
        scatter(missed_DoseTimes4(2:end),Ctrough4,'y*')
        scatter(missed_DoseT(2:end),Ctroughm,'m*')
        scatter(times,peak(:,1),'ko')
        scatter(times1,peak1(:,1),'ro')
        scatter(times2,peak2(:,1),'bo')
        scatter(times3,peak3(:,1),'go')
        scatter(times4,peak4(:,1),'yo')
        scatter(timesm,peakm(:,1),'mo')
        yline(48,'r--');
        yline(12,'r--');
        title('500mg')
        xlabel('Time (hrs)')
        ylabel('Levetiracetam (mg/L)')
        lgd = legend('concentration when dose is on time','concentration when dose late by m/5','concentration when dose late by 2m/5',...
            'concentration when dose late by 3m/5','concentration when dose late by 4m/5','concentration when dose is missed once','ctrough when dose is on time',...
            'ctrough when dose late by m/5', 'ctrough when dose late by 2m/5','ctrough when dose late by 3m/5','ctrough when dose late by 4m/5','ctrough when dose is missed once',...
            'peak conc when dose is on time','peak conc when dose late by m/5', 'peak conc when dose late by 2m/5','peak conc when dose late by 3m/5','peak conc when dose late by 4m/5',...
        'peak conc when dose is missed once','Therapeutic range-lower','Therapeutic range-upper');
        lgd.Location = 'southeastoutside';



        figure;
        auc = [AUC AUC1 AUC2 AUC3 AUC4 AUCm]';
        plot(auc);
        hold on
        scatter([1 2 3 4 5 6],auc,'r');
        set(gca, 'XTick', [1 2 3 4 5 6])
        set(gca, 'XTickLabel', {'AUC when dose is on time' 'AUC when dose late by m/5' 'AUC when dose late by 2m/5' 'AUC when dose late by 3m/5' 'AUC when dose late by 4m/5' 'AUC when dose is missed once'})
        title('AUC')
        xlabel('Different dose conditions')


        figure;
        plot(T,BalanceD,'linewidth',2)
        hold on
        plot(T1,BalanceD1,'r','linewidth',2)
        plot(T2,BalanceD2,'b','linewidth',2)
        plot(T3,BalanceD3,'g','linewidth',2)
        plot(T4,BalanceD4,'y','linewidth',2)
        plot(T4,BalanceDm,'c','linewidth',2)
        lgd = legend('Mass Balance when dose is on time','Mass Balance when dose late by m/5', 'Mass Balance when dose late by 2m/5',...
            'Mass Balance when dose late by 3m/5','Mass Balance when dose late by 4m/5','Mass Balance when dose is missed once');
        lgd.Location = 'southeastoutside';
        title('Mass Balance')
        xlabel('Time (hrs)')
        y = Y(:,1); 
        y1 = Y1(:,1); 
        y2 = Y2(:,1); 
        y3 = Y3(:,1); 
        y4 = Y4(:,1); 
        ym = Ym(:,1);
        peaks = [peak1 peak2 peak3 peak4];
        peakm = peakm(:,1);
        ctroughs = [Ctrough1 Ctrough2 Ctrough3 Ctrough4];
        ct_times = [missed_DoseTimes1(2:end)' missed_DoseTimes2(2:end)' missed_DoseTimes3(2:end)' missed_DoseTimes4(2:end)'];
        times_values = [T1 T2 T3 T4];
        p_times = [times1 times2 times3 times4];
        DoseTimes = DoseTimes(2:end)';
        missed_DoseT = missed_DoseT(2:end);
        change = zeros(1,6);
        auc = [AUC AUC1 AUC2 AUC3 AUC4 AUCm];

        for j=2:6
            change(1,j) = (auc(1)-auc(j))/auc(1);
        end

        change = change(2:end)*100;
        
        save change_means_500.mat change

        save AUc_values_500.mat auc
        save c_p_values_delayed_500.mat ctroughs peaks ct_times p_times
        save c_p_values_missed_500.mat Ctroughm peakm missed_DoseT timesm
        save conc_normal_500.mat T y 
        save conc_delayed_500.mat T1 y1 T2 y2 T3 y3 T4 y4 
        save c_p_normal_500.mat Ctrough peak DoseTimes times
        save conc_missed_500.mat Tm ym
        
        
      case 2
        D0 = 1500;   % mg
        [T1,Y1,BalanceD1,AUC1,Ctrough1,peak1] = LEV_sim_peak(D0,p,missed_DoseTimes1);
        [T2,Y2,BalanceD2,AUC2,Ctrough2,peak2] = LEV_sim_peak(D0,p,missed_DoseTimes2);
        [T3,Y3,BalanceD3,AUC3,Ctrough3,peak3] = LEV_sim_peak(D0,p,missed_DoseTimes3);
        [T4,Y4,BalanceD4,AUC4,Ctrough4,peak4] = LEV_sim_peak(D0,p,missed_DoseTimes4);
        [T,Y,BalanceD,AUC,Ctrough,peak] = LEV_sim_peak(D0,p,DoseTimes);
        [Tm,Ym,BalanceDm,AUCm,Ctroughm,peakm] = LEV_sim_peak(D0,p,missed_DoseT);

        times = zeros(length(peak),1);
        times1 = zeros(length(peak1),1);
        times2 = zeros(length(peak2),1);
        times3 = zeros(length(peak3),1);
        times4 = zeros(length(peak4),1);
        timesm = zeros(length(peakm),1);
        for i=1:length(Y)
            for j=1:10
                for k=1:9
                if Y(i) == peak(j,1)
                    times(j)= T(i);
                end
                if Y1(i) == peak1(j,1)
                    times1(j)= T1(i);
                end
                if Y2(i) == peak2(j,1)
                    times2(j)= T2(i);
                end
                if Y3(i) == peak3(j,1)
                    times3(j)= T3(i);
                end
                if Y4(i) == peak4(j,1)
                    times4(j)= T4(i);
                end
                if Ym(i) == peakm(k,1)
                    timesm(k)= Tm(i);
                end
                end
            end
        end




        figure;
        plot(T,Y(:,1),'k','linewidth',2)
        hold on
        plot(T1,Y1(:,1),'r','linewidth',2)
        plot(T2,Y2(:,1),'b','linewidth',2)
        plot(T3,Y3(:,1),'g','linewidth',2)
        plot(T4,Y4(:,1),'y','linewidth',2)
        plot(Tm,Ym(:,1),'c','linewidth',3)
        scatter(DoseTimes(2:end),Ctrough,'k*')
        scatter(missed_DoseTimes1(2:end),Ctrough1,'r*')
        scatter(missed_DoseTimes2(2:end),Ctrough2,'b*')
        scatter(missed_DoseTimes3(2:end),Ctrough3,'g*')
        scatter(missed_DoseTimes4(2:end),Ctrough4,'y*')
        scatter(missed_DoseT(2:end),Ctroughm,'m*')
        scatter(times,peak(:,1),'ko')
        scatter(times1,peak1(:,1),'ro')
        scatter(times2,peak2(:,1),'bo')
        scatter(times3,peak3(:,1),'go')
        scatter(times4,peak4(:,1),'yo')
        scatter(timesm,peakm(:,1),'mo')
        yline(48,'r--');
        yline(12,'r--');
        title('1500mg')
        xlabel('Time (hrs)')
        ylabel('Levetiracetam (mg/L)')
        lgd = legend('concentration when dose is on time','concentration when dose late by m/5','concentration when dose late by 2m/5',...
            'concentration when dose late by 3m/5','concentration when dose late by 4m/5','concentration when dose is missed once','ctrough when dose is on time',...
            'ctrough when dose late by m/5', 'ctrough when dose late by 2m/5','ctrough when dose late by 3m/5','ctrough when dose late by 4m/5','ctrough when dose is missed once',...
            'peak conc when dose is on time','peak conc when dose late by m/5', 'peak conc when dose late by 2m/5','peak conc when dose late by 3m/5','peak conc when dose late by 4m/5',...
        'peak conc when dose is missed once','Therapeutic range-lower','Therapeutic range-upper');
        lgd.Location = 'southeastoutside';



        figure;
        auc = [AUC AUC1 AUC2 AUC3 AUC4 AUCm]';
        plot(auc);
        hold on
        scatter([1 2 3 4 5 6],auc,'r');
        set(gca, 'XTick', [1 2 3 4 5 6])
        set(gca, 'XTickLabel', {'AUC when dose is on time' 'AUC when dose late by m/5' 'AUC when dose late by 2m/5' 'AUC when dose late by 3m/5' 'AUC when dose late by 4m/5' 'AUC when dose is missed once'})
        title('AUC')
        xlabel('Different dose conditions')


        figure;
        plot(T,BalanceD,'linewidth',2)
        hold on
        plot(T1,BalanceD1,'r','linewidth',2)
        plot(T2,BalanceD2,'b','linewidth',2)
        plot(T3,BalanceD3,'g','linewidth',2)
        plot(T4,BalanceD4,'y','linewidth',2)
        plot(T4,BalanceDm,'c','linewidth',2)
        lgd = legend('Mass Balance when dose is on time','Mass Balance when dose late by m/5', 'Mass Balance when dose late by 2m/5',...
            'Mass Balance when dose late by 3m/5','Mass Balance when dose late by 4m/5','Mass Balance when dose is missed once');
        lgd.Location = 'southeastoutside';
        title('Mass Balance')
        xlabel('Time (hrs)')
        y = Y(:,1); 
        y1 = Y1(:,1); 
        y2 = Y2(:,1); 
        y3 = Y3(:,1); 
        y4 = Y4(:,1); 
        ym = Ym(:,1);
        peaks = [peak1 peak2 peak3 peak4];
        peakm = peakm(:,1);
        ctroughs = [Ctrough1 Ctrough2 Ctrough3 Ctrough4];
        ct_times = [missed_DoseTimes1(2:end)' missed_DoseTimes2(2:end)' missed_DoseTimes3(2:end)' missed_DoseTimes4(2:end)'];
        times_values = [T1 T2 T3 T4];
        p_times = [times1 times2 times3 times4];
        DoseTimes = DoseTimes(2:end)';
        missed_DoseT = missed_DoseT(2:end);
        change = zeros(1,6);
        auc = [AUC AUC1 AUC2 AUC3 AUC4 AUCm];

        for j=2:6
            change(1,j) = (auc(1)-auc(j))/auc(1);
        end

        change = change(2:end)*100;
        save AUc_values_1500.mat auc
        save change_means_1500.mat change
        save c_p_values_delayed_1500.mat ctroughs peaks ct_times p_times
        save c_p_values_missed_1500.mat Ctroughm peakm missed_DoseT timesm
        save conc_normal_1500.mat T y 
        save conc_delayed_1500.mat T1 y1 T2 y2 T3 y3 T4 y4 
        save c_p_normal_1500.mat Ctrough peak DoseTimes times
        save conc_missed_1500.mat Tm ym
       case 3
        D0 = 1000;   % mg
        [T1,Y1,BalanceD1,AUC1,Ctrough1,peak1] = LEV_sim_peak(D0,p,missed_DoseTimes1);
        [T2,Y2,BalanceD2,AUC2,Ctrough2,peak2] = LEV_sim_peak(D0,p,missed_DoseTimes2);
        [T3,Y3,BalanceD3,AUC3,Ctrough3,peak3] = LEV_sim_peak(D0,p,missed_DoseTimes3);
        [T4,Y4,BalanceD4,AUC4,Ctrough4,peak4] = LEV_sim_peak(D0,p,missed_DoseTimes4);
        [T,Y,BalanceD,AUC,Ctrough,peak] = LEV_sim_peak(D0,p,DoseTimes);
        [Tm,Ym,BalanceDm,AUCm,Ctroughm,peakm] = LEV_sim_peak(D0,p,missed_DoseT);

        times = zeros(length(peak),1);
        times1 = zeros(length(peak1),1);
        times2 = zeros(length(peak2),1);
        times3 = zeros(length(peak3),1);
        times4 = zeros(length(peak4),1);
        timesm = zeros(length(peakm),1);
        for i=1:length(Y)
            for j=1:10
                for k=1:9
                if Y(i) == peak(j,1)
                    times(j)= T(i);
                end
                if Y1(i) == peak1(j,1)
                    times1(j)= T1(i);
                end
                if Y2(i) == peak2(j,1)
                    times2(j)= T2(i);
                end
                if Y3(i) == peak3(j,1)
                    times3(j)= T3(i);
                end
                if Y4(i) == peak4(j,1)
                    times4(j)= T4(i);
                end
                if Ym(i) == peakm(k,1)
                    timesm(k)= Tm(i);
                end
                end
            end
        end




        figure;
        plot(T,Y(:,1),'k','linewidth',2)
        hold on
        plot(T1,Y1(:,1),'r','linewidth',2)
        plot(T2,Y2(:,1),'b','linewidth',2)
        plot(T3,Y3(:,1),'g','linewidth',2)
        plot(T4,Y4(:,1),'y','linewidth',2)
        plot(Tm,Ym(:,1),'c','linewidth',3)
        scatter(DoseTimes(2:end),Ctrough,'k*')
        scatter(missed_DoseTimes1(2:end),Ctrough1,'r*')
        scatter(missed_DoseTimes2(2:end),Ctrough2,'b*')
        scatter(missed_DoseTimes3(2:end),Ctrough3,'g*')
        scatter(missed_DoseTimes4(2:end),Ctrough4,'y*')
        scatter(missed_DoseT(2:end),Ctroughm,'m*')
        scatter(times,peak(:,1),'ko')
        scatter(times1,peak1(:,1),'ro')
        scatter(times2,peak2(:,1),'bo')
        scatter(times3,peak3(:,1),'go')
        scatter(times4,peak4(:,1),'yo')
        scatter(timesm,peakm(:,1),'mo')
        yline(48,'r--');
        yline(12,'r--');
        title('1000mg')
        xlabel('Time (hrs)')
        ylabel('Levetiracetam (mg/L)')
        lgd = legend('concentration when dose is on time','concentration when dose late by m/5','concentration when dose late by 2m/5',...
            'concentration when dose late by 3m/5','concentration when dose late by 4m/5','concentration when dose is missed once','ctrough when dose is on time',...
            'ctrough when dose late by m/5', 'ctrough when dose late by 2m/5','ctrough when dose late by 3m/5','ctrough when dose late by 4m/5','ctrough when dose is missed once',...
            'peak conc when dose is on time','peak conc when dose late by m/5', 'peak conc when dose late by 2m/5','peak conc when dose late by 3m/5','peak conc when dose late by 4m/5',...
        'peak conc when dose is missed once','Therapeutic range-lower','Therapeutic range-upper');
        lgd.Location = 'southeastoutside';



        figure;
        auc = [AUC AUC1 AUC2 AUC3 AUC4 AUCm]';
        plot(auc);
        hold on
        scatter([1 2 3 4 5 6],auc,'r');
        set(gca, 'XTick', [1 2 3 4 5 6])
        set(gca, 'XTickLabel', {'AUC when dose is on time' 'AUC when dose late by m/5' 'AUC when dose late by 2m/5' 'AUC when dose late by 3m/5' 'AUC when dose late by 4m/5' 'AUC when dose is missed once'})
        title('AUC')
        xlabel('Different dose conditions')


        figure;
        plot(T,BalanceD,'linewidth',2)
        hold on
        plot(T1,BalanceD1,'r','linewidth',2)
        plot(T2,BalanceD2,'b','linewidth',2)
        plot(T3,BalanceD3,'g','linewidth',2)
        plot(T4,BalanceD4,'y','linewidth',2)
        plot(T4,BalanceDm,'c','linewidth',2)
        lgd = legend('Mass Balance when dose is on time','Mass Balance when dose late by m/5', 'Mass Balance when dose late by 2m/5',...
            'Mass Balance when dose late by 3m/5','Mass Balance when dose late by 4m/5','Mass Balance when dose is missed once');
        lgd.Location = 'southeastoutside';
        title('Mass Balance')
        xlabel('Time (hrs)')
        y = Y(:,1); 
        y1 = Y1(:,1); 
        y2 = Y2(:,1); 
        y3 = Y3(:,1); 
        y4 = Y4(:,1); 
        ym = Ym(:,1);
        peaks = [peak1 peak2 peak3 peak4];
        peakm = peakm(:,1);
        ctroughs = [Ctrough1 Ctrough2 Ctrough3 Ctrough4];
        ct_times = [missed_DoseTimes1(2:end)' missed_DoseTimes2(2:end)' missed_DoseTimes3(2:end)' missed_DoseTimes4(2:end)'];
        times_values = [T1 T2 T3 T4];
        p_times = [times1 times2 times3 times4];
        DoseTimes = DoseTimes(2:end)';
        missed_DoseT = missed_DoseT(2:end);
        change = zeros(1,6);
        auc = [AUC AUC1 AUC2 AUC3 AUC4 AUCm];

        for j=2:6
            change(1,j) = (auc(1)-auc(j))/auc(1);
        end

        change = change(2:end)*100;
        
        save change_means_1000.mat change

        save AUc_values_1000.mat auc
        save c_p_values_delayed_1000.mat ctroughs peaks ct_times p_times
        save c_p_values_missed_1000.mat Ctroughm peakm missed_DoseT timesm
        save conc_normal_1000.mat T y 
        save conc_delayed_1000.mat T1 y1 T2 y2 T3 y3 T4 y4 
        save c_p_normal_1000.mat Ctrough peak DoseTimes times
        save conc_missed_1000.mat Tm ym
       case 4
        D0 = 750;   % mg
        [T1,Y1,BalanceD1,AUC1,Ctrough1,peak1] = LEV_sim_peak(D0,p,missed_DoseTimes1);
        [T2,Y2,BalanceD2,AUC2,Ctrough2,peak2] = LEV_sim_peak(D0,p,missed_DoseTimes2);
        [T3,Y3,BalanceD3,AUC3,Ctrough3,peak3] = LEV_sim_peak(D0,p,missed_DoseTimes3);
        [T4,Y4,BalanceD4,AUC4,Ctrough4,peak4] = LEV_sim_peak(D0,p,missed_DoseTimes4);
        [T,Y,BalanceD,AUC,Ctrough,peak] = LEV_sim_peak(D0,p,DoseTimes);
        [Tm,Ym,BalanceDm,AUCm,Ctroughm,peakm] = LEV_sim_peak(D0,p,missed_DoseT);

        times = zeros(length(peak),1);
        times1 = zeros(length(peak1),1);
        times2 = zeros(length(peak2),1);
        times3 = zeros(length(peak3),1);
        times4 = zeros(length(peak4),1);
        timesm = zeros(length(peakm),1);
        for i=1:length(Y)
            for j=1:10
                for k=1:9
                if Y(i) == peak(j,1)
                    times(j)= T(i);
                end
                if Y1(i) == peak1(j,1)
                    times1(j)= T1(i);
                end
                if Y2(i) == peak2(j,1)
                    times2(j)= T2(i);
                end
                if Y3(i) == peak3(j,1)
                    times3(j)= T3(i);
                end
                if Y4(i) == peak4(j,1)
                    times4(j)= T4(i);
                end
                if Ym(i) == peakm(k,1)
                    timesm(k)= Tm(i);
                end
                end
            end
        end




        figure;
        plot(T,Y(:,1),'k','linewidth',2)
        hold on
        plot(T1,Y1(:,1),'r','linewidth',2)
        plot(T2,Y2(:,1),'b','linewidth',2)
        plot(T3,Y3(:,1),'g','linewidth',2)
        plot(T4,Y4(:,1),'y','linewidth',2)
        plot(Tm,Ym(:,1),'c','linewidth',3)
        scatter(DoseTimes(2:end),Ctrough,'k*')
        scatter(missed_DoseTimes1(2:end),Ctrough1,'r*')
        scatter(missed_DoseTimes2(2:end),Ctrough2,'b*')
        scatter(missed_DoseTimes3(2:end),Ctrough3,'g*')
        scatter(missed_DoseTimes4(2:end),Ctrough4,'y*')
        scatter(missed_DoseT(2:end),Ctroughm,'m*')
        scatter(times,peak(:,1),'ko')
        scatter(times1,peak1(:,1),'ro')
        scatter(times2,peak2(:,1),'bo')
        scatter(times3,peak3(:,1),'go')
        scatter(times4,peak4(:,1),'yo')
        scatter(timesm,peakm(:,1),'mo')
        yline(48,'r--');
        yline(12,'r--');
        title('750mg')
        xlabel('Time (hrs)')
        ylabel('Levetiracetam (mg/L)')
        lgd = legend('concentration when dose is on time','concentration when dose late by m/5','concentration when dose late by 2m/5',...
            'concentration when dose late by 3m/5','concentration when dose late by 4m/5','concentration when dose is missed once','ctrough when dose is on time',...
            'ctrough when dose late by m/5', 'ctrough when dose late by 2m/5','ctrough when dose late by 3m/5','ctrough when dose late by 4m/5','ctrough when dose is missed once',...
            'peak conc when dose is on time','peak conc when dose late by m/5', 'peak conc when dose late by 2m/5','peak conc when dose late by 3m/5','peak conc when dose late by 4m/5',...
        'peak conc when dose is missed once','Therapeutic range-lower','Therapeutic range-upper');
        lgd.Location = 'southeastoutside';



        figure;
        auc = [AUC AUC1 AUC2 AUC3 AUC4 AUCm]';
        plot(auc);
        hold on
        scatter([1 2 3 4 5 6],auc,'r');
        set(gca, 'XTick', [1 2 3 4 5 6])
        set(gca, 'XTickLabel', {'AUC when dose is on time' 'AUC when dose late by m/5' 'AUC when dose late by 2m/5' 'AUC when dose late by 3m/5' 'AUC when dose late by 4m/5' 'AUC when dose is missed once'})
        title('AUC')
        xlabel('Different dose conditions')


        figure;
        plot(T,BalanceD,'linewidth',2)
        hold on
        plot(T1,BalanceD1,'r','linewidth',2)
        plot(T2,BalanceD2,'b','linewidth',2)
        plot(T3,BalanceD3,'g','linewidth',2)
        plot(T4,BalanceD4,'y','linewidth',2)
        plot(T4,BalanceDm,'c','linewidth',2)
        lgd = legend('Mass Balance when dose is on time','Mass Balance when dose late by m/5', 'Mass Balance when dose late by 2m/5',...
            'Mass Balance when dose late by 3m/5','Mass Balance when dose late by 4m/5','Mass Balance when dose is missed once');
        lgd.Location = 'southeastoutside';
        title('Mass Balance')
        xlabel('Time (hrs)')
        y = Y(:,1); 
        y1 = Y1(:,1); 
        y2 = Y2(:,1); 
        y3 = Y3(:,1); 
        y4 = Y4(:,1); 
        ym = Ym(:,1);
        peaks = [peak1 peak2 peak3 peak4];
        peakm = peakm(:,1);
        ctroughs = [Ctrough1 Ctrough2 Ctrough3 Ctrough4];
        ct_times = [missed_DoseTimes1(2:end)' missed_DoseTimes2(2:end)' missed_DoseTimes3(2:end)' missed_DoseTimes4(2:end)'];
        times_values = [T1 T2 T3 T4];
        p_times = [times1 times2 times3 times4];
        DoseTimes = DoseTimes(2:end)';
        missed_DoseT = missed_DoseT(2:end);
        change = zeros(1,6);
        auc = [AUC AUC1 AUC2 AUC3 AUC4 AUCm];

        for j=2:6
            change(1,j) = (auc(1)-auc(j))/auc(1);
        end

        change = change(2:end)*100;
        
        save change_means_750.mat change

        save AUc_values_750.mat auc
        save c_p_values_delayed_750.mat ctroughs peaks ct_times p_times
        save c_p_values_missed_750.mat Ctroughm peakm missed_DoseT timesm
        save conc_normal_750.mat T y 
        save conc_delayed_750.mat T1 y1 T2 y2 T3 y3 T4 y4 
        save c_p_normal_750.mat Ctrough peak DoseTimes times
        save conc_missed_750.mat Tm ym
        
end