clc
clear
close all

% This code restore the data from txt files to mat format
subNUM = 30; % this is subject number

for i=1:2 % this is for group number
    count = 1;
    for j=0:6 % this is path folder should be changed to 6 later
        for k=1:10 % this is trial number
           clear data TDATA
           str = ['C:\Users\Default.ASUS\OneDrive - Arizona State University\Projects\Variable Damping Control\2D basic\Data\Subject',num2str(subNUM),'\group',num2str(i),'\path',num2str(j),'\Trial',num2str(k),'\KukaData.txt'];
           data = load(str);
           TDATA.x = data(:,27);
           TDATA.y = data(:,28);
           TDATA.xinit = data(:,25);
           TDATA.yinit = data(:,26);
           TDATA.xtarget = data(:,23);
           TDATA.ytarget = data(:,24);
           TDATA.Vx = data(:,20);
           TDATA.Vy = -data(:,19);
           TDATA.ax = data(:,22);
           TDATA.ay = -data(:,21);
           TDATA.uintentx = TDATA.Vx.*TDATA.ax;
           TDATA.uintenty = TDATA.Vy.*TDATA.ay;
           TDATA.Bx = data(:,18);
           TDATA.By = data(:,17);
           TDATA.Forcex = data(:,10);
           TDATA.Forcey = -data(:,9);
           
           if i==1
               Celofgroup1{count} = TDATA;
               count = count +1;
           else
               Celofgroup2{count} = TDATA;
               count = count +1;
           end
        end
    end
end

save(['C:\Users\Default.ASUS\OneDrive - Arizona State University\Projects\Variable Damping Control\2D basic\Data\Subject',num2str(subNUM),'\Subject',num2str(subNUM), '.mat'],'Celofgroup1','Celofgroup2');
