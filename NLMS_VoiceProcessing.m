close all;
clear;

% Load audio file
[desired_audio, Fs] = audioread('晚安大小姐 (cut).wav');

% Normalize the desired audio to avoid clipping when adding noise
desired_audio = desired_audio / max(abs(desired_audio));

% Add Gaussian noise to the audio signal
noisy_audio = desired_audio + 0.1 * randn(size(desired_audio));

% Normalize the noisy audio to avoid clipping
noisy_audio = noisy_audio / max(abs(noisy_audio));

% Save the noisy audio to listen to the input
audiowrite('noisy_audio.wav', noisy_audio, Fs);

% Set parameters for NLMS
mu = 0.01; % step size
filter_size = 3; % order of the filter (M)

tic;

% Apply NLMS adaptive filter
filtered_audio = nlms_adaptive_filter(noisy_audio, desired_audio, mu, filter_size);

elapsedtime = toc;
elapsedtimeMS = elapsedtime * 1000;

disp(['Execution time of NLMS: ', num2str(elapsedtimeMS), ' ms']);

% Normalize the filtered audio to avoid clipping
filtered_audio = filtered_audio / max(abs(filtered_audio));

% Save the filtered audio
audiowrite('filtered_audio_NLMS.wav', filtered_audio, Fs);

% Calculate MSE
mse_value = immse(filtered_audio, desired_audio);
disp(['MSE: ', num2str(mse_value)]);

% Calculate PSNR
psnr_value = psnr(filtered_audio, desired_audio);
disp(['PSNR: ', num2str(psnr_value), ' dB']);

% Function definition of NLMS adaptive filter
function output = nlms_adaptive_filter(noisy_audio, desired_audio, mu, filter_size)
    N = length(noisy_audio);
    output = zeros(N, 1);
    error = zeros(N, 1);
    w = zeros(filter_size, 1); % Initialize filter weights

    for n = filter_size:N
        % Extract the current window of the noisy signal
        x = noisy_audio(n:-1:n-filter_size+1);

        % Calculate the filter output
        y = w' * x;

        % Calculate the error signal
        e = desired_audio(n) - y;
        error(n) = e;

        % Tap-weight adaption
        norm_factor = 0.0001 + x' * x; % Avoid division by zero
        w = w + (mu / norm_factor) * x * e;
        
        % Store the output signal
        output(n) = y;
    end
end