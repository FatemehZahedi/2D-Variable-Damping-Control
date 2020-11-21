% This code is for data analysis and will find different measures

clc
clear
close all

subNUM = 23;

load(['C:\Users\Default.ASUS\OneDrive - Arizona State University\Projects\Variable Damping Control\2D basic\Data\Subject',num2str(subNUM),'\Subject',num2str(subNUM), '.mat']);
errorbound = 0.005;

for gp=1:2
    if gp==1
        CELDATA = Celofgroup1;
    else
        CELDATA = Celofgroup2;
    end
    for trial=1:length(CELDATA)
        clear DATA move hit ind_stability
        DATA = CELDATA{trial};
        
        %move = find(((DATA.x-DATA.xinit(10)).^2+(DATA.y-DATA.yinit(10)).^2)>= (errorbound)^2);
        %firstmove = move(1);
        
        move = find(((DATA.x-DATA.xinit(10)).^2+(DATA.y-DATA.yinit(10)).^2)>= (errorbound)^2);
        flag = 0;
        q = 1;
        for f=1:length(move)
                if f < length(move) && flag == 0
                    
                    if move(f+1)-move(f) > 1
                        e_length = f - q;
                        if e_length >= 500
                            firstmove = move(q);
                            flag = 1;
                        else
                            q = f+1;
                        end
                    end
                end
                if f == length(move) && flag == 0
                    firstmove = move(q);
                end
        end
        
        hit = find(((DATA.x-DATA.xtarget(10)).^2+(DATA.y-DATA.ytarget(10)).^2)<= (errorbound)^2);
        firsthit = hit(1);
        if firstmove == 1
            testtrial = trial;
        end
        ind_stability = find(((DATA.x-DATA.xtarget(10)).^2+(DATA.y-DATA.ytarget(10)).^2)<= (errorbound)^2);
        flag = 0;
        q = 1;
        for f=1:length(ind_stability)
                if f < length(ind_stability) && flag == 0
                    
                    if ind_stability(f+1)-ind_stability(f) > 1
                        e_length = f - q;
                        if e_length >= 500
                            tstability = ind_stability(q);
                            flag = 1;
                        else
                            q = f+1;
                        end
                    end
                end
                if f == length(ind_stability) && flag == 0
                    e_length = f - q;
                    if e_length >= 500
                        tstability = ind_stability(q);
                    else
                        tstability = 2000 + firsthit;
                    end
                end
        end
         
        % finding overshoot
        
        clear xrotated yrotated
        
        angleofrotation = atan2(DATA.ytarget(10)-DATA.yinit(10),DATA.xtarget(10)-DATA.xinit(10));
        angle = -angleofrotation;
        xrotated = (DATA.x-DATA.xinit(10))*cos(angle)-(DATA.y-DATA.yinit(10))*sin(angle);
        yrotated = (DATA.y-DATA.yinit(10))*cos(angle)+(DATA.x-DATA.xinit(10))*sin(angle);
        xtargetrotated = (DATA.xtarget(10)-DATA.xinit(10))*cos(angle)-(DATA.ytarget(10)-DATA.yinit(10))*sin(angle);
        ytargetrotated = (DATA.ytarget(10)-DATA.yinit(10))*cos(angle)+(DATA.xtarget(10)-DATA.xinit(10))*sin(angle);
        
        tangentialovershootpercentage = abs((max(abs(xrotated))-xtargetrotated)/xtargetrotated)*100;
        %%%tangentialovershoot = abs((max(abs(xrotated))-abs(xtargetrotated)));
        %normalovershootpercentage = abs((max(abs(yrotated))-ytargetrotated)/ytargetrotated)*100;
        %%%normalovershoot = abs((max(abs(yrotated))-abs(ytargetrotated)));
        tangentialovershoot_0 = abs((max(abs(xrotated))-abs(xtargetrotated))) - errorbound;
        if tangentialovershoot_0 <= 0
            tangentialovershoot = 0;
        else
            tangentialovershoot = tangentialovershoot_0;
        end
        normalovershoot_0 = abs((max(abs(yrotated))-abs(ytargetrotated))) - errorbound;
        if normalovershoot_0 <= 0
            normalovershoot = 0;
        else
            normalovershoot = normalovershoot_0;
        end
        
        
        Tableofmeasures(trial,1,gp) = tangentialovershoot;
        Tableofmeasures(trial,2,gp) = normalovershoot;
        
        % finding enclosing ellipses
        
        clear xrotated yrotated
        
        ellipseTollerance = 0.1;
        anglerotated = [pi/4; pi - pi/4; pi + pi/4; 2*pi - pi];
        
        xcheck = DATA.xtarget(10)-DATA.xinit(10);
        ycheck = DATA.ytarget(10)-DATA.yinit(10);
        
        if xcheck >= 0 && ycheck >= 0
            angle = -anglerotated(1);
        elseif xcheck < 0 && ycheck >=0
            angle = -anglerotated(2);
        elseif xcheck < 0 && ycheck < 0
            angle = -anglerotated(3);
        elseif xcheck >=0 && ycheck < 0
            angle = -anglerotated(4);
        end
        
        xrotated = (DATA.x(firsthit:end))*cos(angle)-(DATA.y(firsthit:end))*sin(angle);
        yrotated = (DATA.y(firsthit:end))*cos(angle)+(DATA.x(firsthit:end))*sin(angle);
        
        figure(1)
        subplot(1,2,gp)
        xlabel('x (m)');
        ylabel('y (m)');
        if gp == 1
            color = 'b';
            title('Positive')
        elseif gp == 2
            color = 'r';
            title('Variable')
        end
        plot(xrotated, yrotated, color)
        
        % Calculate and plot the enclousing ellipse
        [A, C] = MinVolEllipse([xrotated'; yrotated'], ellipseTollerance);
        hold on
        area = Ellipse_plot(A, C);
        
        grid on
        %axis equal
        axis([-0.11 0.15 -0.15 0.1])
        
        Tableofmeasures(trial,3,gp) = area;
        
        % Variability time
        
        variabilitytime = tstability - firsthit;
        Tableofmeasures(trial,4,gp) = variabilitytime;
        
        % mean speed
        velocityboth = sqrt(DATA.Vx.^2 + DATA.Vy.^2);
        meanspeed = mean(velocityboth(firstmove:firsthit));
        
        Tableofmeasures(trial,5,gp) = meanspeed;
        
        % maximum speed
        
        maxspeed = max(velocityboth);
        Tableofmeasures(trial,6,gp) = maxspeed;
        
        % User Effort
        
        ForceRMS = sqrt(0.5*(DATA.Forcex.^2 + DATA.Forcey.^2));
        
        meanRMS = mean(ForceRMS);
        maxRMS = max(ForceRMS);
        
        Tableofmeasures(trial,7,gp) = meanRMS;
        Tableofmeasures(trial,8,gp) = maxRMS;
        
        %Normalization of Data
        appending_num = 5000;
        ave_num = 200;
        if trial == 48
            xsample = DATA.x(1:end);
            ysample = DATA.y(1:end);
            uintentxsample = DATA.uintentx(1:end);
            uintentysample = DATA.uintenty(1:end);
            Bxsample = DATA.Bx(1:end);
            Bysample = DATA.By(1:end);
        end
        
        xnormalized = (DATA.x(firstmove:end)-DATA.xinit(10))/abs(xcheck);
        ynormalized = (DATA.y(firstmove:end)-DATA.yinit(10))/abs(ycheck);
        if xcheck >= 0
            normalized1x(trial,:,1,gp) = [xnormalized' mean(xnormalized(end-ave_num:end))*ones(1,appending_num-length(xnormalized))];
            normalized1x(trial,:,2,gp) = [DATA.uintentx(firstmove:end)' mean(DATA.uintentx(end-ave_num:end))*ones(1,appending_num-length(xnormalized))]/((xcheck)^2);
            normalized1x(trial,:,3,gp) = [DATA.Bx(firstmove:end)' mean(DATA.Bx(end-ave_num:end))*ones(1,appending_num-length(xnormalized))];
            FirstMOVE(trial,1,gp) = firstmove;
        elseif xcheck < 0
            normalized2x(trial,:,1,gp) = [xnormalized' mean(xnormalized(end-ave_num:end))*ones(1,appending_num-length(xnormalized))];
            normalized2x(trial,:,2,gp) = [DATA.uintentx(firstmove:end)' mean(DATA.uintentx(end-ave_num:end))*ones(1,appending_num-length(xnormalized))]/((xcheck)^2);
            normalized2x(trial,:,3,gp) = [DATA.Bx(firstmove:end)' mean(DATA.Bx(end-ave_num:end))*ones(1,appending_num-length(xnormalized))];
            FirstMOVE(trial,2,gp) = firstmove;
        end
        if ycheck >= 0
           normalized1y(trial,:,1,gp) = [ynormalized' mean(ynormalized(end-ave_num:end))*ones(1,appending_num-length(ynormalized))];
           normalized1y(trial,:,2,gp) = [DATA.uintenty(firstmove:end)' mean(DATA.uintenty(end-ave_num:end))*ones(1,appending_num-length(ynormalized))]/((ycheck)^2);
           normalized1y(trial,:,3,gp) = [DATA.By(firstmove:end)' mean(DATA.By(end-ave_num:end))*ones(1,appending_num-length(ynormalized))];
        elseif ycheck < 0
           normalized2y(trial,:,1,gp) = [ynormalized' mean(ynormalized(end-ave_num:end))*ones(1,appending_num-length(ynormalized))];
           normalized2y(trial,:,2,gp) = [DATA.uintenty(firstmove:end)' mean(DATA.uintenty(end-ave_num:end))*ones(1,appending_num-length(ynormalized))]/((ycheck)^2);
           normalized2y(trial,:,3,gp) = [DATA.By(firstmove:end)' mean(DATA.By(end-ave_num:end))*ones(1,appending_num-length(ynormalized))];
        end
           
    end
end
        
% averaging and standard deviation

for gp =1:2
    for i=1:size(Tableofmeasures,2)
        AveTable(gp,i,1) = mean(Tableofmeasures(:,i,gp));
        AveTable(gp,i,2) = std(Tableofmeasures(:,i,gp));
    end
end

%%
% Plotting 
%close
figure(2)

% Overshoot
subplot(2,3,[1 2])
dir1 = [1 2];
label = {'Tangential','Normal'};
bar(dir1,AveTable(:,1:2,1)',1)
set(gca, 'XTickLabel',label, 'XTick',1:numel(dir1))
ngroups = size(AveTable(:,1:2,1)', 1);
nbars = size(AveTable(:,1:2,1)', 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
hold on
y=AveTable(:,1:2,1)';
error= AveTable(:,1:2,2)';
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x,y(:,i),error(:,i) , '.k');
end
legend('Positive','Variable')
title('Overshoot')
ylabel('Overshoot (m)');
axis([0 3 -0.003 0.02]);

% Enclosing Ellipses
subplot(2,3,3)
for j=1:2
bar(j,AveTable(j,3,1))
set(gca, 'XTickLabel',label, 'XTick',[])
hold on
error= AveTable(j,3,2);
errorbar(j,AveTable(j,3,1),error , '.k');
end


title('Enclosing Ellipses')
ylabel('Area (m^2)');
axis([0 3 0 1.2e-4]);

% Variability Time
subplot(2,3,4)
for j=1:size(AveTable(:,4,1)*10e-3,1)
bar(j,AveTable(j,4,1)*10e-3)
set(gca, 'XTickLabel',label, 'XTick',[])
hold on
error= AveTable(j,4,2)*10e-3;
errorbar(j,AveTable(j,4,1)*10e-3,error , '.k');
end


title('Variability Time')
ylabel('Time (s)');
axis([0 3 -2 13]);

% Mean speed
subplot(2,3,5)
for j=1:size(AveTable(:,5,1),1)
bar(j,AveTable(j,5,1))
set(gca, 'XTickLabel',label, 'XTick',[])
hold on
error= AveTable(j,5,2);
errorbar(j,AveTable(j,5,1),error , '.k');
end


title('Mean Speed')
ylabel('Velocity (m/s)');
axis([0 3 0 0.3]);

% Max speed
subplot(2,3,6)
for j=1:size(AveTable(:,6,1),1)
bar(j,AveTable(j,6,1))
set(gca, 'XTickLabel',label, 'XTick',[])
hold on
error= AveTable(j,6,2);
errorbar(j,AveTable(j,6,1),error , '.k');
end


title('Max Speed')
ylabel('Velocity (m/s)');
axis([0 3 0 0.6]);

% User Effort

% figure(3)
% dir1 = [1 2];
% label = {'MeanRMS','MaxRMS'};
% bar(dir1,AveTable(:,7:8,1)',1)
% set(gca, 'XTickLabel',label, 'XTick',1:numel(dir1))
% ngroups = size(AveTable(:,7:8,1)', 1);
% nbars = size(AveTable(:,7:8,1)', 2);
% % Calculating the width for each bar group
% groupwidth = min(0.8, nbars/(nbars + 1.5));
% hold on
% y=AveTable(:,7:8,1)';
% error= AveTable(:,7:8,2)';
% for i = 1:nbars
%     x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
%     errorbar(x,y(:,i),error(:,i) , '.k');
% end
% legend('Positive','Variable')
% title('Interaction Force')
% ylabel('Force (N)');
%%
figure(3)
subplot(1,2,1)
for j=1:size(AveTable(:,7,1),1)
bar(j,AveTable(j,7,1))
set(gca, 'XTickLabel',label, 'XTick',[])
hold on
error= AveTable(j,7,2);
errorbar(j,AveTable(j,7,1),error , '.k');
end
ylabel('Force (N)');
xlabel('MeanRMS');
title('Interaction Force')

subplot(1,2,2)
for j=1:size(AveTable(:,8,1),1)
p(j) = bar(j,AveTable(j,8,1));
set(gca, 'XTickLabel',label, 'XTick',[])
hold on
error= AveTable(j,8,2);
errorbar(j,AveTable(j,8,1),error , '.k');
end

xlabel('MaxRMS');
legend([p(1) p(2)],'Positive','Variable')
title('Interaction Force')

%%
% Saving data

load('groupave.mat');
groupdata(:,:,:,subNUM)=Tableofmeasures;
groupavedata(:,:,:,subNUM) = AveTable;
save('groupave.mat','groupdata','groupavedata');
%%
% Plotting position profile

time = 1:size(normalized1x,2);

for i=1:2
    for j=1:2
        FirstMOVE2(:,i,j) = FirstMOVE(find(FirstMOVE(:,i,j)~=0),i,j);
    end
end

figure(4)
t=suptitle('Position Profile');


for direction = 1:4
    if direction == 1
        normalized = normalized1x;
        dir = 1;
    elseif direction == 2
        normalized = normalized2x;
        dir = 2;
    elseif direction == 3
        normalized = normalized1y;
        dir = 1;
    elseif direction == 4
        normalized = normalized2y;
        dir = 2;
    end
    for gp = 1:2
        
        if gp == 1
            color = 'b';
        elseif gp == 2
            color = 'r';
        end
        
        position = normalized(:,:,1,gp);
        position(~any(position,2),:,:)=[];
        if gp == 2
            userintent = normalized(:,:,2,gp);
            userintent(~any(userintent,2),:,:)=[];
            damping = normalized(:,:,3,gp);
            damping(~any(damping,2),:,:)=[];
        end
        subplot(4,4,direction+gp*(gp-1)*2)
        plot(time*0.001, mean(position,1),color);
        hold on
        plot(time*0.001,mean(position,1)+std(position),'-.k');
        hold on
        plot(time*0.001,mean(position,1)-std(position),'-.k');
        if direction+gp*(gp-1)*2 == 1 || direction+gp*(gp-1)*2 == 2
            title ('x (m)');
        elseif direction+gp*(gp-1)*2 == 3 || direction+gp*(gp-1)*2 == 4
            title ('y (m)');
        end
        if direction+gp*(gp-1)*2 == 1
            ylabel('Positive')
        elseif direction+gp*(gp-1)*2 == 5
            ylabel('Variable')
        end
        set(gca, 'XTickLabel',label, 'XTick',[])
        
        if gp == 2
            subplot(4,4,direction+gp*4)
            plot(time*0.001, mean(userintent,1),'color',[0.4940 0.1840 0.5560]);
            hold on
            plot(time*0.001,mean(userintent,1)+std(userintent),'-.k');
            hold on
            plot(time*0.001,mean(userintent,1)-std(userintent),'-.k');
            if direction ==1
                ylabel('User Intent (m^2/s^3)')
            end
            set(gca, 'XTickLabel',label, 'XTick',[])
            
            subplot(4,4,direction+gp*6)
            plot(time*0.001, mean(damping,1),'color',[0.4940 0.1840 0.5560]);
            hold on
            plot(time*0.001,mean(damping,1)+std(userintent),'-.k');
            hold on
            plot(time*0.001,mean(damping,1)-std(userintent),'-.k');
            xlabel('Time (s)')
            if direction ==1
                ylabel('damping (Ns/m)')
            end
        end
    end
end

%%
% plotting single trials

figure(5)
t=suptitle('Position Profile');

for direction = 1:4
    if direction == 1
        normalized = normalized1x;
    elseif direction == 2
        normalized = normalized2x;
    elseif direction == 3
        normalized = normalized1y;
    elseif direction == 4
        normalized = normalized2y;
    end
    for gp = 1:2
        
        if gp == 1
            color = 'b';
        elseif gp == 2
            color = 'r';
        end
        
        position = normalized(:,:,1,gp);
        position(~any(position,2),:,:)=[];
        if gp == 2
            userintent = normalized(:,:,2,gp);
            userintent(~any(userintent,2),:,:)=[];
            damping = normalized(:,:,3,gp);
            damping(~any(damping,2),:,:)=[];
        end
        subplot(4,4,direction+gp*(gp-1)*2)
        for trial =1:size(position,1)
            hold on
            plot(time*0.001, position(trial,:),'color',[0.7 0.7 0.7]);
        end
        plot(time*0.001, mean(position,1),color);
        if direction+gp*(gp-1)*2 == 1 || direction+gp*(gp-1)*2 == 2
            title ('x (m)');
        elseif direction+gp*(gp-1)*2 == 3 || direction+gp*(gp-1)*2 == 4
            title ('y (m)');
        end
        if direction+gp*(gp-1)*2 == 1
            ylabel('Positive')
        elseif direction+gp*(gp-1)*2 == 5
            ylabel('Variable')
        end
        set(gca, 'XTickLabel',label, 'XTick',[])
        
        if gp == 2
            subplot(4,4,direction+gp*4)
            for trial =1:size(position,1)
                hold on
                plot(time*0.001, userintent(trial,:),'color',[0.7 0.7 0.7]);
            end
            plot(time*0.001, mean(userintent,1),'color',[0.4940 0.1840 0.5560]);
            if direction ==1
                ylabel('User Intent (m^2/s^3)')
            end
            set(gca, 'XTickLabel',label, 'XTick',[])
            
            subplot(4,4,direction+gp*6)
            for trial =1:size(position,1)
                hold on
                plot(time*0.001, damping(trial,:),'color',[0.7 0.7 0.7]);
            end
            plot(time*0.001, mean(damping,1),'color',[0.4940 0.1840 0.5560]);
            xlabel('Time (s)')
            if direction ==1
                ylabel('damping (Ns/m)')
            end
        end
    end
end


%%
time1 = 1:size(xsample,1);
figure
subplot(3,2,1)
plot(time1*0.001,xsample)
title('ML')
ylabel('position profile (m)')
subplot(3,2,2)
plot(time1*0.001,ysample)
title('AP')
subplot(3,2,3)
plot(time1*0.001,uintentxsample)
ylabel('user intent (m^2/s^3)')
subplot(3,2,4)
plot(time1*0.001,uintentysample)
subplot(3,2,5)
plot(time1*0.001,Bxsample)
ylabel('Damping (Ns/m)')
xlabel('Time (s)');
subplot(3,2,6)
plot(time1*0.001,Bysample)
xlabel('Time (s)');
%%

figure
plot(xsample,ysample)
xlabel('x (m)')
ylabel('y (m)')


        
    
    
        


        
        
        
        
        
        
        
        