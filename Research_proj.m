% Time and Signal Setup
% Reduce noise
t = linspace(0, 1, 1000);         % 1 second duration, 1000 samples
ecg_clean = 1.5*sin(2*pi*5*t) + 0.5*sin(2*pi*15*t); % Simulated clean ECG
ecg_noisy = ecg_clean + 0.4*randn(size(t));        % Add random noise

% Fourier Series Approximation
K = 50;  % Number of terms (Increase this for more details)
a0 = mean(ecg_noisy);            % DC component (average of the signal)
reconstructed = a0 * ones(size(t));  % Start with the DC component

%Fourier Series reconstruction
for k = 1:K
    cosk = cos(2*pi*k*t);           % Cosine function for kth harmonic
    sink = sin(2*pi*k*t);           % Sine function for kth harmonic
    
    % Compute Fourier coefficients ak and bk
    ak = 2 * mean(ecg_noisy .* cosk);
    bk = 2 * mean(ecg_noisy .* sink);
    
    % Add the kth harmonic to the reconstruction
    reconstructed = reconstructed + ak*cosk + bk*sink;
end


figure;
subplot(3,1,1);
plot(t, ecg_clean); 
title('Original ECG'); 
ylabel('Amplitude');
subplot(3,1,2);
plot(t, ecg_noisy); 
title('Noisy ECG'); 
ylabel('Amplitude');
subplot(3,1,3); 
plot(t, reconstructed); 
title(['Reconstructed ECG with ', num2str(K), ' Fourier terms']); 
xlabel('Time (s)'); 
ylabel('Amplitude');


%%
%Data Compression
t = linspace(0, 1, 1000); % Time vector
dt = t(2) - t(1); % Time step
ecg = 1.5*sin(2*pi*5*t) + 0.5*sin(2*pi*15*t); % same ECG clean signal

T = 1; % Period
K = 10; % Number of Fourier terms

a0 = (1/T) * sum(ecg) * dt;
reconstructed = a0 * ones(size(t));

for k = 1:K
    ak = (2/T) * sum(ecg .* cos(2*pi*k*t/T)) * dt;
    bk = (2/T) * sum(ecg .* sin(2*pi*k*t/T)) * dt;
    reconstructed = reconstructed + ak*cos(2*pi*k*t/T) + bk*sin(2*pi*k*t/T);
end


plot(t, ecg, 'k--', t, reconstructed, 'r', 'LineWidth', 1.2);
legend('Original ECG', ['Compressed with ', num2str(K), ' terms']);
title('ECG Data Compression using Fourier Series');
xlabel('Time (s)'); ylabel('Amplitude');

%%
%Cardiac disorders
t = linspace(0, 1, 1000); 
ecg_clean = 1.5*sin(2*pi*5*t) + 0.5*sin(2*pi*15*t); 

% Simulate a disorder: Ventricular Fibrillation (VF)
vf_noise = 0.5 * randn(size(t)); % Random noise to simulate VF

% Add the VF noise to the clean ECG signal
ecg_vf = ecg_clean + vf_noise;

figure;
subplot(2,1,1);
plot(t, ecg_clean, 'k-', 'LineWidth', 1.5);
title('Original Clean ECG');
xlabel('Time (s)'); ylabel('Amplitude');

subplot(2,1,2);
plot(t, ecg_vf, 'r-', 'LineWidth', 1.5);
title('ECG with Ventricular Fibrillation (VF)');
xlabel('Time (s)'); ylabel('Amplitude');



