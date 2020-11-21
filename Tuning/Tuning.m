%clc
clear
close all

% Define the following values based on the subject number and block tuning
% number
blocktuning = 1; % 1 is when you find tuning by zero damping, 2 is variable damping, 3 is variable damping
subnum = 1;


if blocktuning == 1
    group = 0;
elseif blocktuning == 2
    group = 2;
else
    group = 2;
end

if blocktuning == 3
    countpath = 1;
else
    countpath = 0;
end

for i=0+countpath:2:2+countpath
    for j =1:10
        str = ['C:\Users\Default.ASUS\OneDrive - Arizona State University\Projects\Variable Damping Control\2D basic\Data\Subject',num2str(subnum),'\Tunning\group',num2str(group),'\path',num2str(i),'\Trial',num2str(j),'\KukaData.txt'];
        data = load(str);
        if i==0 || i==1
            n = 0;
        else
            n = 1;
        end
        
        xdot = data(:,20 - n);
        xdotdot = data(:,22 - n);
        
        [kp(n + 1,j), kn(n + 1,j)] = GetKs(xdot, xdotdot);
    end
end

kpave = mean(kp,2);
knave = mean(kn,2);

disp("kpml: ");
disp(kpave(1));
disp("knml: ");
disp(knave(1));

disp("kpap: ");
disp(kpave(2));
disp("knap: ");
disp(knave(2));