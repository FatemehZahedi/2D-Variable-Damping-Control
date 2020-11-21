%This code is for group data analysis
clc
clear
close all

load('groupave.mat');
subjects = [7 8 9 10 11 13 14 17 22 24 26 27 30 33 34 35];
%subjects = [33 8 14 7 11 35 13 22 10 9 27 30 23 32 19];% 26 20 17 21 29 25 28 34 18 24];
%subjects = [33 8 14 7 11 35 13 22 10 9 27 30 19 26 21];
for gp=1:2
    for subnum=1:length(subjects)%7:size(groupavedata,4)
        avedata(subnum,:,gp) = groupavedata(gp,:,1,subjects(subnum));
    end
end

% reordering for mean value and standard deviation
for gp =1:2
    for i=1:size(avedata,2)
        gpave(gp,i,1) = mean(avedata(:,i,gp));
        gpave(gp,i,2) = std(avedata(:,i,gp));
    end
end


%%
% Plotting 

figure(1)

% Overshoot
subplot(2,3,[1 2])
dir1 = [1 2];
label = {'Tangential','Normal'};
bar(dir1,gpave(:,1:2,1)',1)
set(gca, 'XTickLabel',label, 'XTick',1:numel(dir1))
ngroups = size(gpave(:,1:2,1)', 1);
nbars = size(gpave(:,1:2,1)', 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
hold on
y=gpave(:,1:2,1)';
error= gpave(:,1:2,2)';
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x,y(:,i),error(:,i) , '.k');
end
legend('Positive','Variable')
title('Overshoot')
ylabel('Overshoot (m)');
axis([0 3 -0.002 0.02]);

% Enclosing Ellipses
subplot(2,3,3)
for j=1:2
bar(j,gpave(j,3,1))
set(gca, 'XTickLabel',label, 'XTick',[])
hold on
error= gpave(j,3,2);
errorbar(j,gpave(j,3,1),error , '.k');
end


title('Enclosing Ellipses')
ylabel('Area (m^2)');
axis([0 3 0 1.1e-4]);

% Variability Time
subplot(2,3,4)
for j=1:size(gpave(:,4,1)*10e-3,1)
bar(j,gpave(j,4,1)*10e-3)
set(gca, 'XTickLabel',label, 'XTick',[])
hold on
error= gpave(j,4,2)*10e-3;
errorbar(j,gpave(j,4,1)*10e-3,error , '.k');
end


title('Variability Time')
ylabel('Time (s)');
axis([0 3 -2 12]);

% Mean speed
subplot(2,3,5)
for j=1:size(gpave(:,5,1),1)
bar(j,gpave(j,5,1))
set(gca, 'XTickLabel',label, 'XTick',[])
hold on
error= gpave(j,5,2);
errorbar(j,gpave(j,5,1),error , '.k');
end


title('Mean Speed')
ylabel('Velocity (m/s)');
axis([0 3 0 0.3]);

% Max speed
subplot(2,3,6)
for j=1:size(gpave(:,6,1),1)
bar(j,gpave(j,6,1))
set(gca, 'XTickLabel',label, 'XTick',[])
hold on
error= gpave(j,6,2);
errorbar(j,gpave(j,6,1),error , '.k');
end


title('Max Speed')
ylabel('Velocity (m/s)');
axis([0 3 0 0.6]);

%%
figure(2)
subplot(1,2,1)
for j=1:size(gpave(:,7,1),1)
bar(j,gpave(j,7,1))
set(gca, 'XTickLabel',label, 'XTick',[])
hold on
error= gpave(j,7,2);
errorbar(j,gpave(j,7,1),error , '.k');
end
ylabel('Force (N)');
xlabel('MeanRMS');
title('Interaction Force')
axis([0 3 0 3]);

subplot(1,2,2)
for j=1:size(gpave(:,8,1),1)
p(j) = bar(j,gpave(j,8,1));
set(gca, 'XTickLabel',label, 'XTick',[])
hold on
error= gpave(j,8,2);
errorbar(j,gpave(j,8,1),error , '.k');
end

xlabel('MaxRMS');
legend([p(1) p(2)],'Positive','Variable')
title('Interaction Force')
axis([0 3 0 30]);

%%

figure

% Overshoot
subplot(3,2,[1 2])
dir1 = [1 2];
label = {'Tangential','Normal'};
bar(dir1,gpave(:,1:2,1)'*100,1)
set(gca, 'XTickLabel',label, 'XTick',1:numel(dir1))
ngroups = size(gpave(:,1:2,1)', 1);
nbars = size(gpave(:,1:2,1)', 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
hold on
y=gpave(:,1:2,1)'*100;
error= gpave(:,1:2,2)'*100;
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x,y(:,i),zeros(size(error(:,i))),error(:,i) , '.k');
end
legend('Positive','Variable')
title('Overshoot')
ylabel('Overshoot (cm)');
axis([0 3 -0.002*100 0.02*100]);
box off

% % Enclosing Ellipses
% subplot(2,3,3)
% for j=1:2
% bar(j,gpave(j,3,1))
% set(gca, 'XTickLabel',label, 'XTick',[])
% hold on
% error= gpave(j,3,2);
% errorbar(j,gpave(j,3,1),error , '.k');
% end
% 
% 
% title('Enclosing Ellipses')
% ylabel('Area (m^2)');
% axis([0 3 0 1.1e-4]);
% 
% % Variability Time
% subplot(2,3,4)
% for j=1:size(gpave(:,4,1)*10e-3,1)
% bar(j,gpave(j,4,1)*10e-3)
% set(gca, 'XTickLabel',label, 'XTick',[])
% hold on
% error= gpave(j,4,2)*10e-3;
% errorbar(j,gpave(j,4,1)*10e-3,error , '.k');
% end
% 
% 
% title('Variability Time')
% ylabel('Time (s)');
% axis([0 3 -2 12]);

% Mean speed
subplot(3,2,3)
for j=1:size(gpave(:,5,1),1)
bar(j,gpave(j,5,1)*100)
set(gca, 'XTickLabel',label, 'XTick',[])
hold on
error= gpave(j,5,2)*100;
errorbar(j,gpave(j,5,1)*100,zeros(size(error)),error , '.k');
end


title('Mean Speed')
ylabel('Velocity (cm/s)');
axis([0 3 0 0.3*100]);
box off

% Max speed
subplot(3,2,4)
for j=1:size(gpave(:,6,1),1)
bar(j,gpave(j,6,1)*100)
set(gca, 'XTickLabel',label, 'XTick',[])
hold on
error= gpave(j,6,2)*100;
errorbar(j,gpave(j,6,1)*100,zeros(size(error)),error , '.k');
end


title('Max Speed')
%ylabel('Velocity (m/s)');
axis([0 3 0 0.6*100]);
box off

%force

subplot(3,2,5)
for j=1:size(gpave(:,7,1),1)
bar(j,gpave(j,7,1))
set(gca, 'XTickLabel',label, 'XTick',[])
hold on
error= gpave(j,7,2);
errorbar(j,gpave(j,7,1),zeros(size(error)),error , '.k');
end
ylabel('Force (N)');
xlabel('MeanRMS');
title('Interaction Force')
axis([0 3 0 3]);
box off

subplot(3,2,6)
for j=1:size(gpave(:,8,1),1)
p(j) = bar(j,gpave(j,8,1));
set(gca, 'XTickLabel',label, 'XTick',[])
hold on
error= gpave(j,8,2);
errorbar(j,gpave(j,8,1),zeros(size(error)),error , '.k');
end

xlabel('MaxRMS');
%legend([p(1) p(2)],'Positive','Variable')
title('Interaction Force')
axis([0 3 0 27]);
box off

%%
% For Dr. Lee
figure

% Overshoot
% subplot(3,2,[1 2])
dir1 = [1 2];
label = {'Tangential','Normal'};
bar(dir1,gpave(:,1:2,1)'*100,1)
set(gca, 'XTickLabel',label, 'XTick',1:numel(dir1))
ngroups = size(gpave(:,1:2,1)', 1);
nbars = size(gpave(:,1:2,1)', 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
hold on
y=gpave(:,1:2,1)'*100;
error= gpave(:,1:2,2)'*100;
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x,y(:,i),zeros(size(error(:,i))),error(:,i) , '.k');
end
legend('Positive','Variable')
title('Overshoot')
ylabel('Overshoot (cm)');
axis([0 3 -0.002*100 0.02*100]);
box off

% Mean speed
figure
for j=1:size(gpave(:,5,1),1)
bar(j,gpave(j,5,1)*100)
set(gca, 'XTickLabel',label, 'XTick',[])
hold on
error= gpave(j,5,2)*100;
errorbar(j,gpave(j,5,1)*100,zeros(size(error)),error , '.k');
end


title('Mean Speed')
ylabel('Velocity (cm/s)');
axis([0 3 0 0.3*100]);
box off

% Max speed
figure
for j=1:size(gpave(:,6,1),1)
bar(j,gpave(j,6,1)*100)
set(gca, 'XTickLabel',label, 'XTick',[])
hold on
error= gpave(j,6,2)*100;
errorbar(j,gpave(j,6,1)*100,zeros(size(error)),error , '.k');
end


title('Max Speed')
ylabel('Velocity (m/s)');
axis([0 3 0 0.6*100]);
box off

%force

figure
for j=1:size(gpave(:,7,1),1)
bar(j,gpave(j,7,1))
set(gca, 'XTickLabel',label, 'XTick',[])
hold on
error= gpave(j,7,2);
errorbar(j,gpave(j,7,1),zeros(size(error)),error , '.k');
end
ylabel('Force (N)');
xlabel('MeanRMS');
title('Interaction Force')
axis([0 3 0 3]);
box off

figure
for j=1:size(gpave(:,8,1),1)
p(j) = bar(j,gpave(j,8,1));
set(gca, 'XTickLabel',label, 'XTick',[])
hold on
error= gpave(j,8,2);
errorbar(j,gpave(j,8,1),zeros(size(error)),error , '.k');
end

xlabel('MaxRMS');
ylabel('Force (N)');
%legend([p(1) p(2)],'Positive','Variable')
title('Interaction Force')
axis([0 3 0 27]);
box off


%%
load('final.mat');
load('final_int.mat');

subjectnum = {'subject7'; 'subject8'; 'subject9'; 'subject10'; 'subject11'; 'subject13'; 'subject14'; 'subject17';...
    'subject18'; 'subject19'; 'subject20'; 'subject21'; 'subject22'; 'subject23'; 'subject24'; 'subject25'; 'subject26';...
    'subject27'; 'subject28'; 'subject29'; 'subject30'; 'subject32'; 'subject33'; 'subject34'; 'subject35'};
subjectsT = [7 8 9 10 11 13 14 17 18 19 20 21 22 23 24 25 26 27 28 29 30 32 33 34 35];
for gp=1:2
    for subnum=1:length(subjectsT)
        avedataT(subnum,:,gp) = groupavedata(gp,:,1,subjectsT(subnum));
    end
end
tangentialOvershoot_pos = round(avedataT(:,1,1)*100*100,0)/100;
normalOvershoot_pos = round(avedataT(:,2,1)*100*100,0)/100;
MeanSpeed_pos = round(avedataT(:,5,1)*100*100,0)/100;
MaxSpeed_pos = round(avedataT(:,6,1)*100*100,0)/100;
MeanRMSForce_pos = round(avedataT(:,7,1)*100,0)/100;
MaxRMSForce_pos = round(avedataT(:,8,1)*100,0)/100;

tangentialOvershoot_var = round(avedataT(:,1,2)*100*100,0)/100;
normalOvershoot_var = round(avedataT(:,2,2)*100*100,0)/100;
MeanSpeed_var = round(avedataT(:,5,2)*100*100,0)/100;
MaxSpeed_var = round(avedataT(:,6,2)*100*100,0)/100;
MeanRMSForce_var = round(avedataT(:,7,2)*100,0)/100;
MaxRMSForce_var = round(avedataT(:,8,2)*100,0)/100;
EMG_BL_MovementPh = round(final(:,2)*100,0)/100;
EMG_BR_MovementPh = round(final(:,3)*100,0)/100;
EMG_FL_MovementPh = round(final(:,4)*100,0)/100;
EMG_FR_MovementPh = round(final(:,5)*100,0)/100;
EMG_BL_StabilityPh = round(final(:,6)*100,0)/100;
EMG_BR_StabilityPh = round(final(:,7)*100,0)/100;
EMG_FL_StabilityPh = round(final(:,8)*100,0)/100;
EMG_FR_StabilityPh = round(final(:,9)*100,0)/100;
EMG_BL_Integral = round(final_int(:,2)*100,0)/100;
EMG_BR_Integral = round(final_int(:,3)*100,0)/100;
EMG_FL_Integral = round(final_int(:,4)*100,0)/100;
EMG_FR_Integral = round(final_int(:,5)*100,0)/100;


TableDATA = table(subjectnum,tangentialOvershoot_pos,normalOvershoot_pos,MeanSpeed_pos,MaxSpeed_pos,MeanRMSForce_pos,MaxRMSForce_pos,...
    tangentialOvershoot_var,normalOvershoot_var,MeanSpeed_var,MaxSpeed_var,MeanRMSForce_var,MaxRMSForce_var,...
    EMG_BL_MovementPh,EMG_BR_MovementPh,EMG_FL_MovementPh,EMG_FR_MovementPh,EMG_BL_StabilityPh,EMG_BR_StabilityPh,...
    EMG_FL_StabilityPh,EMG_FR_StabilityPh,EMG_BL_Integral,EMG_BR_Integral,EMG_FL_Integral,EMG_FR_Integral);

writetable(TableDATA,'Data.xlsx');

