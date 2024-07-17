clc;
clear;
close all;


%Meassage and carrier signal generation
F2=50; %Message signal frequency
F1=4000; %Triangular carrier signal frequency
A=15; %Carrier signal amplitude
t=0:0.000001:.039999;
c=(A*sawtooth(2*pi*F1*t,0.5)); %Carrier signal
figure(1);
subplot(3,1,1);
plot(t,c,'g');
xlabel('time');
ylabel('Amplitude');
title('Carrier triangular wave');
axis([0 0.02 -20 20]);
grid on;
hold on;
m=(0.9*A*sin(2*pi*F2*t)); %Message signal
plot(t,m);
xlabel('Time');
ylabel('Amplitude');
title('Reference and triangular signal');
grid on;
n=length(c);

%Comparing message and triangular signals and calculate on-off times of pwms
for i=1:n 
if (m(i)>c(i))
pwm1(i)=1;
pwm2(i) = ~pwm1(i);
else
pwm1(i)=0;
pwm2(i) = ~pwm1(i);
end
end
k = 0;
l = 0;
for j=1:length(t)
if((pwm1(j)==1) && (k == 0))
k = 1;
tstart = t(j);
end
if((pwm1(j)==0) && (k ==1))
k = 0;
tstop= t(j);
ton = (tstop - tstart)*10^6; %calculate pwm on time
end
if((pwm2(j)==1) && (l == 0))
l = 1;
tstart1 = t(j);
end
if((pwm2(j)==0) && (l ==1))
l = 0;
tstop1= t(j);
toff = (tstop1 - tstart1)*10^6; %calculate pwm off time
end
end
subplot(3,1,2);
plot(t,pwm1,'r');
xlabel('Time');
ylabel('Amplitude');
title('plot of PWM1 for G1 and G2');
axis([0 0.02 -1 1]);
grid on;
hold on;
subplot(3,1,3);
plot(t,pwm2,'b');
xlabel('Time');
ylabel('Amplitude');
title('plot of PWM2 for G3 and G4');
axis([0 0.02 -1 1]);
grid on;

%Plot of voltage accross H-bridge(Vab)
va = 24*pwm1;
vb = -24*pwm2;
vab = va + vb;
figure(2);
plot(0:0.000001:.039999,vab);
title('Vab');
xlabel('Time');
ylabel('Amplitude');
axis([0 0.01 -25 25]);

%pwm signal(Vab) plot
Fs = 1/0.000001;
L = length(vab);
Y = fft(vab)/L;
f = linspace(-Fs/2-Fs/L,Fs/2,L);
figure(3);
stem(f,2*abs(fftshift(Y)),'color','g','Linewidth',2);
title('Single-Sided Amplitude Spectrum of Vab');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
axis([0 16000 0 25]);

%Filter design
wp = .0008*pi;
ws = .0012*pi;
tr = ws - wp;
m = ceil(8*pi/tr)+1;
hd = ideallp(0.0010*pi,m);
y = filter(hd,[1],vab);
figure(4)
plot(t-.01,y,'color','g','Linewidth',2);
xlabel('Time');
ylabel('Amplitude');
axis([0 0.03 -25 25]);
grid on;

%PSD of vab (before filter)
N=L;
fs=Fs;
xdft = Y*L;
xdft = xdft(1:N/2+1);
psdx = ((1/(fs*N))*abs(xdft).^2);
psdx(2:end-1) = 2*psdx(2:end-1);
psdx11=psdx./max(psdx);
freq = 0:fs/length(vab):fs/2;
figure(5)
plot(freq,psdx11,'color','g','Linewidth',2)
grid on
title("PSD(Normalized) Before Filtering")
xlabel("Frequency (Hz)")
ylabel("Power")
axis([0 16000 0 2]);

%THD of vab before filter
XX = abs(fft(vab)).^2;
power_bin = find(XX==max(XX), 1 )
signal_power = XX(power_bin); 
dist_power = sum([2:1:100000]*power_bin);
THD_before = (sqrt(dist_power)/sqrt(signal_power))*100

%THD after filter
xdft2=fft(y);
datavec = xdft2(51:50:251);
THD_after = norm(datavec(2:end),2)/norm(datavec(1),2)

%PSD of y(t) (after filter)
N=L;
fs=Fs;
xdft3 = xdft2;
xdft3 = xdft3(1:N/2+1);
psdx2 = ((1/(fs*N))*abs(xdft3).^2);
psdx2(2:end-1) = 2*psdx2(2:end-1);
psdx22=psdx2./max(psdx2);
freq2 = 0:fs/length(y):fs/2;
figure(6)
plot(freq2,psdx22,'color','g','Linewidth',2)
grid on
title("PSD(Normalized) After Filtering")
xlabel("Frequency (Hz)")
ylabel("Power")
axis([0 10000 0 2]);

function hd = ideallp(wc,m)
a = (m-1)/2;
n = [0:1:(m-1)];
p = n - a;
hd = wc/pi*sinc(wc*p/pi);
end