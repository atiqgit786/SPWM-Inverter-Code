clc;
clear;
close all;


%Meassage and carrier signal generation
F2=50;%message signal frequency
F1=4000;%triangular carrier signal frequency
A=15;%Carrier signal amplitude
t=0:0.000001:.039999;
c=15-abs(A*sawtooth(2*pi*F1*t,0.5));%Carrier signal
m=abs(0.9*A*sin(2*pi*F2*t)); %Message signal
n=length(c); 

%Comparing message and triangular signals,calculate on-off times of pwms
%and  plotting og gate pulses
for i=1:n 
if (m(i)>c(i))
pwm(i)=1;
else
pwm(i)=0;
end
end

xlabel('Time');
ylabel('Amplitude');
title('plot of PWM');
grid on;
k = 0;
l = 0;
for j=1:length(t)
if((pwm(j)==1) && (k == 0))
k = 1;
tstart = t(j);
end
if((pwm(j)==0) && (k ==1))
k = 0;
l = 1;
tstop= t(j);
ton = (tstop - tstart)*10^6; %calculate pwm on time
end
if((pwm(j) ==1) && (l==1))
l = 0;
tstart2 = t(j);
toff = (tstart2 - tstop) * 10^6; %calculate pwm off time
end
end
t1 = 0:0.000001:.009999;
t2 = 0.01:0.000001:.019999;
t3 = 0.02:0.000001:.029999;
t4 = 0.03:0.000001:.039999;
for i = 1:length(t1)
g2(i) = 1;
end
for i =length(t1)+1 :length(t1)+length(t2)+1
g2(i) = 0;
end
for i = 2+length(t1)+length(t2):length(t3)+length(t1)+length(t2)
g2(i) = 1;
end
for i =length(t3)+1+length(t1)+length(t2) :length(t3)+length(t4)+length(t1)+length(t2)
g2(i) = 0;
end
figure(1);
subplot(2,1,1);
plot(0:0.000001:.039999,g2,'color','r','Linewidth',3);
xlabel('Time');
ylabel('Amplitude');
title('G2');
axis([0 0.04 0 1]);
for i = 1:length(t1)+1
g4(i) = 0;
end
for i =length(t1)+2 :length(t1)+length(t2)
g4(i) = 1;
end
for i = length(t1)+length(t2):length(t3)+1+length(t1)+length(t2)
g4(i) = 0;
end
for i =length(t3)+2+length(t1)+length(t2) :length(t3)+length(t4)+length(t1)+length(t2)
g4(i) = 1;
end
subplot(2,1,2);
plot(0:0.000001:0.039999,g4,'color','r','Linewidth',3);
xlabel('Time');
ylabel('Amplitude');
title('G4');
axis([0 0.04 0 1]);
g1 = pwm & g2;
figure(2);
subplot(3,1,1);
plot(0:0.000001:.039999,g1);
xlabel('Time');
ylabel('Amplitude');
title('G1');
axis([0 0.04 0 1]);
g3 = pwm & g4;
subplot(3,1,2);
plot(0:0.000001:.039999,g3);
xlabel('Time');
ylabel('Amplitude');
title('G3');
axis([0 0.04 0 1]);

%Plot of voltage accross H-bridge(Vab)
va = 24*g1;
vb = -24*g3;
vab = va + vb;
subplot(3,1,3);
plot(0:0.000001:.039999,vab);
xlabel('Time');
ylabel('Amplitude');
title('Vab');
axis([0 0.04 -25 25]);

%pwm signal(Vab) plot
Fs = 1/0.000001;
L = length(vab);
Y = fft(vab)/L;
f = linspace(-Fs/2-Fs/L,Fs/2,L);
figure(3);
stem(f,2*abs(fftshift(Y)),'color','r','Linewidth',2);
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
plot(t-.01,y,'color','r','Linewidth',2);
xlabel('Time');
ylabel('Amplitude');
axis([0 0.03 -25 25]);
grid on;

%PSD of vab (before filtering)
N=L;
fs=Fs;
xdft = Y*L;
xdft = xdft(1:N/2+1);
psdx = ((1/(fs*N))*abs(xdft).^2);
psdx(2:end-1) = 2*psdx(2:end-1);
psdx11=psdx./max(psdx);
freq = 0:fs/length(vab):fs/2;
figure(5)
plot(freq,psdx11,'color','r','Linewidth',1)
grid on
title("PSD(Normalized) Before Filtering")
xlabel("Frequency (Hz)")
ylabel("Power")
axis([0 16000 0 2]);

%THD of Vab before filter
XX = abs(fft(vab)).^2;
power_bin = find(XX==max(XX), 1 );
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
plot(freq2,psdx22,'color','r','Linewidth',1)
grid on
title("PSD(Normalized) After Filtering")
xlabel("Frequency (Hz)")
ylabel("Power")
axis([0 16000 0 2]);

function hd = ideallp(wc,m)
a = (m-1)/2;
n = [0:1:(m-1)];
p = n - a;
hd = wc/pi*sinc(wc*p/pi);
end