clc
clear
close all

directorypath = "C:\Users\Default.ASUS\OneDrive - Arizona State University\Projects\Variable Damping Control\2D basic\Data\Subject103\group2\path2\Trial4\KukaData.txt";

data = load(directorypath);
ycenter = 0.76;
%d1 = designfilt('lowpassiir','FilterOrder',12,'HalfPowerFrequency',20,'DesignMethod','butter','Samplerate',1000);

x = data(:,13)-ycenter;
y = -data(:,11);

vxfilt = data(:,20);
vyfilt = -data(:,19);
axfilt = data(:,22);
ayfilt = -data(:,21);

vx = diff(x)/0.001;
vy = diff(y)/0.001;

%ax = filtfilt(d1,diff(vx)/0.001);
ax = diff(vx)/0.001;
ay = diff(vy)/0.001;

userintent_x_filt = vxfilt.*axfilt;
userintent_y_filt = vyfilt.*ayfilt;

userintent_x = vx(1:end-1).*ax;
userintent_y = vy(1:end-1).*ay;

figure %related to x

subplot(4,1,1)
plot(x);
ylabel('x (m)')
subplot(4,1,2)
plot(vx);
hold on 
plot(vxfilt)
legend('unfiltered','filtered')
ylabel('velocity (m/s)')
subplot(4,1,3)
plot(ax)
hold on
plot(axfilt);
legend('unfiltered','filtered')
ylabel('acceleration (m/s^2)')
subplot(4,1,4)
plot(userintent_x)
hold on
plot(userintent_x_filt)
legend('unfiltered','filtered')
ylabel('user intent (m^2/s^3)')
xlabel('Time (ms)')

figure %related to y

subplot(4,1,1)
plot(y);
ylabel('y (m)')
subplot(4,1,2)
plot(vy);
hold on 
plot(vyfilt)
legend('unfiltered','filtered')
ylabel('velocity (m/s)')
subplot(4,1,3)
plot(ay)
hold on
plot(ayfilt);
legend('unfiltered','filtered')
ylabel('acceleration (m/s^2)')
subplot(4,1,4)
plot(userintent_y)
hold on
plot(userintent_y_filt)
legend('unfiltered','filtered')
ylabel('user intent (m^2/s^3)')
xlabel('Time (ms)')

% recalculating damping

Busedx = data(:,18);
Busedy = data(:,17);

