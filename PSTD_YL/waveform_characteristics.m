% Look at Waveform in time
figure;
plot(n_use.*deltat,Ez_source);
grid on;
xlabel('time (seconds)');
ylabel('Amplitude of waveform');

% Spectral calculations
N = numel(Ez_source);
f_s = 1/deltat;         %Sampling frequency
NFFT = N; 

Y = fft(Ez_source,NFFT);
P2 = abs(Y/N);
P1 = P2(1:round(N/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
f = f_s*(0:round(N/2))/N;
phase_function = phase(Y/N);
phase_function = phase_function(1:round(N/2)+1);

% Look at Waveform in spectral domain
figure;
plot(f./1e6,10*log10(P1));
grid on;
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (MHz)')
ylabel('|P1(f)|')

figure;
plot(f./1e6,phase_function);
grid on;
title('Single-Sided Phase Spectrum of X(t)')
xlabel('f (MHz)')
ylabel('Angle P1(f)')
