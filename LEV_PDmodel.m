
%% TIME TO EVENT

clear all; close all; clc;

% Dose amounts based on Karatza et al (2020)
BIDDose = [2000 1500 1000 750 500 250]; % mg
DailyDose = 2 * BIDDose;                % mg
nDosingRegimens = length(DailyDose);

% Tepop = time to event (population estimate)
Tepop = 0.736; % months

% nPatients = number of virtual patients to run
nPatients = 100;

% Tesim = time to event for each virtual patient
Tesim = exprnd(Tepop,nDosingRegimens,nPatients).*(transpose(DailyDose)/2000).^(-2.2);
Tesim_avg = mean(Tesim,2);
Tesim_med = median(Tesim,2);

% Parameters for visualization
colors = ["#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#4DBEEE"];
fontSize = 20;

%% Figure 1: Visualize time to event histograms by dosing regimen (subplots)

f1 = figure(1);
set(f1,'Position',get(0,'Screensize'));
sgtitle(sprintf("Time to Event (Te) Histograms by Dosing Regimen\n"),'FontSize',24,'FontWeight','bold');

i = 1; % 2000 mg BID
subplot(2,3,i); hold on;
histogram(Tesim(i,:),'BinWidth',1,'Normalization','probability','FaceColor',colors(i),'FaceAlpha',0.25);
xline(Tesim_avg(i),'LineWidth',2,'Color',colors(i),'Alpha',1);
text(Tesim_avg(i)+0.5,0.5,strcat('Mean = ',num2str(Tesim_avg(i))),'FontSize',14,'Color',colors(i),'FontWeight','bold');
xline(Tesim_med(i),'--','LineWidth',2,'Color',colors(i),'Alpha',1);
text(Tesim_med(i)+0.5,0.6,strcat('Median = ',num2str(Tesim_med(i))),'FontSize',14,'Color',colors(i),'FontWeight','bold');
xlim([0 20]); ylim([0 1]);
ax = gca; ax.FontSize = fontSize;
title('2000 mg BID','FontSize',fontSize);
xlabel('Te (months)','FontSize',fontSize);
ylabel('Probability','FontSize',fontSize);
hold off;

i = 2; % 1500 mg BID
subplot(2,3,i); hold on;
histogram(Tesim(i,:),'BinWidth',1,'Normalization','probability','FaceColor',colors(i),'FaceAlpha',0.25);
xline(Tesim_avg(i),'LineWidth',2,'Color',colors(i),'Alpha',1);
text(Tesim_avg(i)+0.5,0.5,strcat('Mean = ',num2str(Tesim_avg(i))),'FontSize',14,'Color',colors(i),'FontWeight','bold');
xline(Tesim_med(i),'--','LineWidth',2,'Color',colors(i),'Alpha',1);
text(Tesim_med(i)+0.5,0.6,strcat('Median = ',num2str(Tesim_med(i))),'FontSize',14,'Color',colors(i),'FontWeight','bold');
xlim([0 20]); ylim([0 1]);
ax = gca; ax.FontSize = fontSize;
title('1500 mg BID','FontSize',fontSize);
xlabel('Te (months)','FontSize',fontSize);
ylabel('Probability','FontSize',fontSize);
hold off;

i = 3; % 1000 mg BID
subplot(2,3,i); hold on;
histogram(Tesim(i,:),'BinWidth',1,'Normalization','probability','FaceColor',colors(i),'FaceAlpha',0.25);
xline(Tesim_avg(i),'LineWidth',2,'Color',colors(i),'Alpha',1);
text(Tesim_avg(i)+0.5,0.5,strcat('Mean = ',num2str(Tesim_avg(i))),'FontSize',14,'Color',colors(i),'FontWeight','bold');
xline(Tesim_med(i),'--','LineWidth',2,'Color',colors(i),'Alpha',1);
text(Tesim_med(i)+0.5,0.6,strcat('Median = ',num2str(Tesim_med(i))),'FontSize',14,'Color',colors(i),'FontWeight','bold');
xlim([0 20]); ylim([0 1]);
ax = gca; ax.FontSize = fontSize;
title('1000 mg BID','FontSize',fontSize);
xlabel('Te (months)','FontSize',fontSize);
ylabel('Probability','FontSize',fontSize);
hold off;

i = 4; % 750 mg BID
subplot(2,3,i); hold on;
histogram(Tesim(i,:),'BinWidth',1,'Normalization','probability','FaceColor',colors(i),'FaceAlpha',0.25);
xline(Tesim_avg(i),'LineWidth',2,'Color',colors(i),'Alpha',1);
text(Tesim_avg(i)+0.5,0.5,strcat('Mean = ',num2str(Tesim_avg(i))),'FontSize',14,'Color',colors(i),'FontWeight','bold');
xline(Tesim_med(i),'--','LineWidth',2,'Color',colors(i),'Alpha',1);
text(Tesim_med(i)+0.5,0.6,strcat('Median = ',num2str(Tesim_med(i))),'FontSize',14,'Color',colors(i),'FontWeight','bold');
xlim([0 20]); ylim([0 1]);
ax = gca; ax.FontSize = fontSize;
title('750 mg BID','FontSize',fontSize);
xlabel('Te (months)','FontSize',fontSize);
ylabel('Probability','FontSize',fontSize);
hold off;

i = 5; % 500 mg BID
subplot(2,3,i); hold on;
histogram(Tesim(i,:),'BinWidth',1,'Normalization','probability','FaceColor',colors(i),'FaceAlpha',0.25);
xline(Tesim_avg(i),'LineWidth',2,'Color',colors(i),'Alpha',1);
text(Tesim_avg(i)+0.5,0.5,strcat('Mean = ',num2str(Tesim_avg(i))),'FontSize',14,'Color',colors(i),'FontWeight','bold');
xline(Tesim_med(i),'--','LineWidth',2,'Color',colors(i),'Alpha',1);
text(Tesim_med(i)+0.5,0.6,strcat('Median = ',num2str(Tesim_med(i))),'FontSize',14,'Color',colors(i),'FontWeight','bold');
xlim([0 20]); ylim([0 1]);
ax = gca; ax.FontSize = fontSize;
title('500 mg BID','FontSize',fontSize);
xlabel('Te (months)','FontSize',fontSize);
ylabel('Probability','FontSize',fontSize);
hold off;

i = 6; % 250 mg BID
subplot(2,3,i); hold on;
histogram(Tesim(i,:),'BinWidth',1,'Normalization','probability','FaceColor',colors(i),'FaceAlpha',0.25);
xline(Tesim_avg(i),'LineWidth',2,'Color',colors(i),'Alpha',1);
text(Tesim_avg(i)+0.5,0.5,strcat('Mean = ',num2str(Tesim_avg(i))),'FontSize',14,'Color',colors(i),'FontWeight','bold');
xline(Tesim_med(i),'--','LineWidth',2,'Color',colors(i),'Alpha',1);
text(Tesim_med(i)+0.5,0.6,strcat('Median = ',num2str(Tesim_med(i))),'FontSize',14,'Color',colors(i),'FontWeight','bold');
xlim([0 20]); ylim([0 1]);
ax = gca; ax.FontSize = fontSize;
title('250 mg BID','FontSize',fontSize);
xlabel('Te (months)','FontSize',fontSize);
ylabel('Probability','FontSize',fontSize);
hold off;

saveas(f1,'PD_Figure_1.png');

%% Figure 2: Visualize time to event histograms by dosing regimen (stacked)

f2 = figure(2);
set(f2,'Position',get(0,'Screensize'));

hold on;
histogram(Tesim(1,:),'BinWidth',1,'Normalization','probability');
histogram(Tesim(2,:),'BinWidth',1,'Normalization','probability');
histogram(Tesim(3,:),'BinWidth',1,'Normalization','probability');
histogram(Tesim(4,:),'BinWidth',1,'Normalization','probability');
histogram(Tesim(5,:),'BinWidth',1,'Normalization','probability');
histogram(Tesim(6,:),'BinWidth',1,'Normalization','probability');
hold off;

xlim([0 14]); ylim([0 1]);
ax = gca; ax.FontSize = fontSize;
title('Time to Event (Te) Histograms by Dosing Regimen','FontSize',24,'FontWeight','bold');
xlabel('Te (months)','FontSize',fontSize);
ylabel('Probability','FontSize',fontSize);
l = legend('2000 mg BID','1500 mg BID','1000 mg BID','750 mg BID','500 mg BID','250 mg BID');
title(l,'Dosing Regimens');

saveas(f2,'PD_Figure_2.png');

%% Figure 3: Visualize Kaplan-Meier survival curves by dosing regimen (subplots)

x = 1:0.1:50;

f3 = figure(3);
set(f3,'Position',get(0,'Screensize'));
sgtitle(sprintf("Kaplan-Meier Survival Curves by Dosing Regimen\n"),'FontSize',24,'FontWeight','bold');

i = 1; % 2000 mg BID
subplot(2,3,i); hold on;
ecdf(floor(Tesim(i,:)+1),'function','survivor','Alpha',0.05,'Bounds','on');
expsurv = 1-cdf('exponential',x-0.5,Tesim_avg(i));
plot(x,expsurv,'--','LineWidth',0.5,'Color','black');
xlim([0 20]); ylim([0 1]);
ax = gca; ax.FontSize = fontSize;
title('2000 mg BID','FontSize',fontSize);
xlabel('Te (months)','FontSize',fontSize);
ylabel('Probability','FontSize',fontSize);
legend('Survival Curve','Lower CB','Upper CB','Exponential Fit');
hold off;

i = 2; % 1500 mg BID
subplot(2,3,i); hold on; plot(-1);
ecdf(floor(Tesim(i,:)+1),'function','survivor','Alpha',0.05,'Bounds','on');
expsurv = 1-cdf('exponential',x-0.5,Tesim_avg(i));
plot(x,expsurv,'--','LineWidth',0.5,'Color','black');
xlim([0 20]); ylim([0 1]);
ax = gca; ax.FontSize = fontSize;
title('1500 mg BID','FontSize',fontSize);
xlabel('Te (months)','FontSize',fontSize);
ylabel('Probability','FontSize',fontSize);
legend('Survival Curve','Lower CB','Upper CB','Exponential Fit');
hold off;

i = 3; % 1000 mg BID
subplot(2,3,i); hold on; plot(-1); plot(-1);
ecdf(floor(Tesim(i,:)+1),'function','survivor','Alpha',0.05,'Bounds','on');
expsurv = 1-cdf('exponential',x-0.5,Tesim_avg(i));
plot(x,expsurv,'--','LineWidth',0.5,'Color','black');
xlim([0 20]); ylim([0 1]);
ax = gca; ax.FontSize = fontSize;
title('1000 mg BID','FontSize',fontSize);
xlabel('Te (months)','FontSize',fontSize);
ylabel('Probability','FontSize',fontSize);
legend('Survival Curve','Lower CB','Upper CB','Exponential Fit');
hold off;

i = 4; % 750 mg BID
subplot(2,3,i); hold on; plot(-1); plot(-1); plot(-1);
ecdf(floor(Tesim(i,:)+1),'function','survivor','Alpha',0.05,'Bounds','on');
expsurv = 1-cdf('exponential',x-0.5,Tesim_avg(i));
plot(x,expsurv,'--','LineWidth',0.5,'Color','black');
xlim([0 20]); ylim([0 1]);
ax = gca; ax.FontSize = fontSize;
title('750 mg BID','FontSize',fontSize);
xlabel('Te (months)','FontSize',fontSize);
ylabel('Probability','FontSize',fontSize);
legend('Survival Curve','Lower CB','Upper CB','Exponential Fit');
hold off;

i = 5; % 500 mg BID
subplot(2,3,i); hold on; plot(-1); plot(-1); plot(-1); plot(-1);
ecdf(floor(Tesim(i,:)+1),'function','survivor','Alpha',0.05,'Bounds','on');
expsurv = 1-cdf('exponential',x-0.5,Tesim_avg(i));
plot(x,expsurv,'--','LineWidth',0.5,'Color','black');
xlim([0 20]); ylim([0 1]);
ax = gca; ax.FontSize = fontSize;
title('500 mg BID','FontSize',fontSize);
xlabel('Te (months)','FontSize',fontSize);
ylabel('Probability','FontSize',fontSize);
legend('Survival Curve','Lower CB','Upper CB','Exponential Fit');
hold off;

i = 6; % 250 mg BID
subplot(2,3,i); hold on; plot(-1); plot(-1); plot(-1); plot(-1); plot(-1);
ecdf(floor(Tesim(i,:)+1),'function','survivor','Alpha',0.05,'Bounds','on');
expsurv = 1-cdf('exponential',x-0.5,Tesim_avg(i));
plot(x,expsurv,'--','LineWidth',0.5,'Color','black');
xlim([0 20]); ylim([0 1]);
ax = gca; ax.FontSize = fontSize;
title('250 mg BID','FontSize',fontSize);
xlabel('Te (months)','FontSize',fontSize);
ylabel('Probability','FontSize',fontSize);
legend('Survival Curve','Lower CB','Upper CB','Exponential Fit');
hold off;
    
saveas(f3,'PD_Figure_3.png');

%% Figure 4: Visualize Kaplan-Meier survival curves by dosing regimen (stacked)

f4 = figure(4);
set(f4,'Position',get(0,'Screensize'));

hold on;

ecdf(floor(Tesim(1,:)+1),'function','survivor');
ecdf(floor(Tesim(2,:)+1),'function','survivor');
ecdf(floor(Tesim(3,:)+1),'function','survivor');
ecdf(floor(Tesim(4,:)+1),'function','survivor');
ecdf(floor(Tesim(5,:)+1),'function','survivor');
ecdf(floor(Tesim(6,:)+1),'function','survivor');

expsurv = 1-cdf('exponential',x-0.5,Tesim_avg(1));
plot(x,expsurv,'--','LineWidth',2,'Color',colors(1));

expsurv = 1-cdf('exponential',x-0.5,Tesim_avg(2));
plot(x,expsurv,'--','LineWidth',2,'Color',colors(2));

expsurv = 1-cdf('exponential',x-0.5,Tesim_avg(3));
plot(x,expsurv,'--','LineWidth',2,'Color',colors(3));

expsurv = 1-cdf('exponential',x-0.5,Tesim_avg(4));
plot(x,expsurv,'--','LineWidth',2,'Color',colors(4));

expsurv = 1-cdf('exponential',x-0.5,Tesim_avg(5));
plot(x,expsurv,'--','LineWidth',2,'Color',colors(5));

expsurv = 1-cdf('exponential',x-0.5,Tesim_avg(6));
plot(x,expsurv,'--','LineWidth',2,'Color',colors(6));

hold off;

xlim([0 14]); ylim([0 1]);
ax = gca; ax.FontSize = fontSize;
title('Kaplan-Meier Survival Curves by Dosing Regimen','FontSize',24,'FontWeight','bold');
xlabel('Te (months)','FontSize',fontSize);
ylabel('Probability','FontSize',fontSize);
l = legend('2000 mg BID Surv','1500 mg BID Surv','1000 mg BID Surv','750 mg BID Surv','500 mg BID Surv','250 mg BID Surv',...
    '2000 mg BID Fit','1500 mg BID Fit','1000 mg BID Fit','750 mg BID Fit','500 mg BID Fit','250 mg BID Fit');
title(l,'Dosing Regimens');

saveas(f4,'PD_Figure_4.png');

%% Export data

Tesim_T = transpose(Tesim);
Tesim_avg_T = transpose(Tesim_avg);
Tesim_med_T = transpose(Tesim_med);

save('Tesim.mat','Tesim_T');
save('Tesim_avg.mat','Tesim_avg_T');
save('Tesim_med.mat','Tesim_med_T');

survCurv1 = transpose(floor(Tesim(1,:)+1));
survCurv2 = transpose(floor(Tesim(2,:)+1));
survCurv3 = transpose(floor(Tesim(3,:)+1));
survCurv4 = transpose(floor(Tesim(4,:)+1));
survCurv5 = transpose(floor(Tesim(5,:)+1));
survCurv6 = transpose(floor(Tesim(6,:)+1));
survCurvs = [survCurv1 survCurv2 survCurv3 survCurv4 survCurv5 survCurv6];
save('survCurvs.mat','survCurvs');

survCurv1_fit = transpose(1-cdf('exponential',x-0.5,Tesim_avg(1)));
survCurv2_fit = transpose(1-cdf('exponential',x-0.5,Tesim_avg(2)));
survCurv3_fit = transpose(1-cdf('exponential',x-0.5,Tesim_avg(3)));
survCurv4_fit = transpose(1-cdf('exponential',x-0.5,Tesim_avg(4)));
survCurv5_fit = transpose(1-cdf('exponential',x-0.5,Tesim_avg(5)));
survCurv6_fit = transpose(1-cdf('exponential',x-0.5,Tesim_avg(6)));
survCurvs_fit = [transpose(x) survCurv1_fit survCurv2_fit survCurv3_fit survCurv4_fit survCurv5_fit survCurv6_fit];
save('survCurvs_fit.mat','survCurvs_fit');
