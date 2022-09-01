clc; clear all; close all;

%{
Drawing the PSD of the ECG signal. Applying a moving average 
filter to this signal first, then applying low pass FIR with two 
window functions and IIR filters with Butterworth and Chebychev methods.
%}


%a
fs = 1000;
ecg_2 = importdata('ecg_2.mat');
ECG_2 = fftshift(abs(fft(ecg_2,length(ecg_2))));
f=linspace(-fs/2,fs/2,length(ecg_2)); 

t=0:1/fs:8.567;


figure
plot(f,(ECG_2)/length(ecg_2))
title('Frequency spectrum of ecg_2')
xlabel('f')
ylabel('magnitude')


%{
 Generally, parts with a magnitude below 0.1 can be called noise. 
 In addition, additional signals around frequency components that thicken 
 these lines can also be counted as noise. All signals other than those 
 around 0 Hz are also noise.
%}

%c

corr = xcorr(ecg_2);
periodogram(corr,[],length(corr),fs,'centered');

[Pxx1 f1]=periodogram(ecg_2, rectwin(length(ecg_2)), length(ecg_2), fs);
[Pxx2 f1]=pwelch(ecg_2, rectwin(300),150, length(ecg_2), fs);

figure
subplot(211)
plot(f1, 10*log10(Pxx1));
grid
title('Periodogram estimate of the PSD')
xlabel('Frequency')
ylabel('Power/frequency')

subplot(212)
plot(f1, 10*log10(Pxx2));
grid
title('Welch estimate of the PSD')
xlabel('Frequency')
ylabel('Power/frequency')


%d


%i
B = 1/10*ones(10,1);
out = filter(B,1,ecg_2);
figure
plot(t,ecg_2)
hold on
plot(t,out)
title('original signal and MA filtered signal')
grid on
legend('original signal','MA filtered signal')

%{
As stated in the previous explanation, signals other than 0 Hz are 
counted as noise. In addition, there is some noise at 0 Hz.
%}

%ii
OUT = fftshift(abs(fft(out,length(out))));
figure
plot(f,ECG_2/length(ECG_2))
hold on
plot(f,OUT/length(OUT))
title('frequency spectrums of original signal and MA filtered signal')
grid on
legend('original signal','MA filtered signal')




%e

%i
Wn1=80/500;
N = 35;
kaiser = fir1(N,Wn1,kaiser(N+1,0.5));
kaiserfiltered = filter(kaiser,1,ecg_2);
KAISERFILTERED =fftshift(fft(kaiserfiltered,length(kaiserfiltered)));     

Wn2 = 60/500;
bartlett = fir1(N,Wn2,bartlett(N+1));  
bartlettfiltered = filter(bartlett,1,ecg_2);
BARTLETTFILTERED = fftshift(fft(bartlettfiltered,length(bartlettfiltered)));    

%ii

figure
plot(t,kaiserfiltered,'g')
hold on
plot(t,ecg_2,'b')
legend('filtered signal','original signal')
xlabel('time')
title ('Kaiser Window')
ylabel ('amplitude');
figure
freqz(kaiser,1)
sum1 = sum((kaiserfiltered-ecg_2).^2)/length(ecg_2);
 
figure
plot(t,bartlettfiltered,'k')
hold on
plot(t,ecg_2,'b')
legend('filtered signal','original signal')
xlabel('time')
title ('Bartlett Window')
ylabel ('amplitude');
figure
freqz(bartlett,1)
sum2 = sum((bartlettfiltered-ecg_2).^2)/length(ecg_2);

%{
i have observed that the kaiser filter of the same rating 
is better at filtering out noise.
%}

%iii
figure
plot(f,abs(BARTLETTFILTERED)/length(ecg_2),'linewidth',1)
hold on 
plot(f,abs(ECG_2)/length(ecg_2))
xlabel('frequency')
title ('Bartlett Window')
ylabel ('amplitude');
legend('filtered signal','original signal')

figure
plot(f,abs(KAISERFILTERED)/length(ecg_2),'b')
hold on 
plot(f,abs(ECG_2)/length(ecg_2),'g')
xlabel('frequency')
title ('Kaiser Window')
ylabel ('amplitude');
legend('filtered signal','original signal')



%f
%i
Wn3 = 40/500;
N2 = 7;

[b,a]=butter(N2,Wn3,'low');
lowfilt = filter(b,a,ecg_2);
LOWFILT = fftshift(fft(lowfilt,length(lowfilt)));     


[b,a] = cheby2(N2,40,Wn3,'low'); 
lowfilt2 = filter(b,a,ecg_2);
LOWFILT2 = fftshift(fft(lowfilt2,length(lowfilt2))); 

%ii
figure
plot(t,lowfilt2,'r')
hold on
plot(t,ecg_2,'g')
xlabel('time')
title ('CHEBYCHEV 2')
ylabel ('amplitude');
legend('filtered signal','original signal')

figure
freqz(b,a)

figure
plot(t,lowfilt,'g')
hold on
plot(t,ecg_2,'b')
xlabel('time')
title ('BUTTERWORTH FILTER')
ylabel ('amplitude');
legend('filtered signal','original signal')

%{
I have observed that a butterworth filter of the same rating 
is better at filtering out noise.
%}

%iii
figure
plot(f,abs(LOWFILT2)/length(ecg_2),'linewidth',3)
hold on
plot(f,abs(ECG_2)/length(ecg_2))
xlabel('frequency')
title ('CHEBYCHEV 2 Filter')
ylabel ('amplitude');
legend('filtered signal','original signal')
figure
freqz(b,a)

figure
plot(f,abs(LOWFILT)/length(ecg_2),'r')
hold on 
plot(f,abs(ECG_2)/length(ecg_2))
xlabel('frequency')
title ('Butterworth Filter')
ylabel ('amplitude');
legend('filtered signal','original signal')
