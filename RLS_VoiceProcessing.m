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

% Set parameters for RLS
lambda = 0.99; % forgetting factor
delta = 0.5; % regularization parameter
filter_size = 3; % order of the filter (M)

tic;

% Apply RLS adaptive filter
filtered_audio = rls_adaptive_filter(noisy_audio, desired_audio, lambda, delta, filter_size);

elapsedtime = toc;
elapsedtimeMS = elapsedtime * 1000;

disp(['Execution time of RLS: ', num2str(elapsedtimeMS), ' ms']);

% Normalize the filtered audio to avoid clipping
filtered_audio = filtered_audio / max(abs(filtered_audio));

% Save the filtered audio
audiowrite('filtered_audio_RLS.wav', filtered_audio, Fs);

% Calculate MSE
mse_value = immse(filtered_audio, desired_audio);
disp(['MSE: ', num2str(mse_value)]);

% Calculate PSNR
psnr_value = psnr(filtered_audio, desired_audio);
disp(['PSNR: ', num2str(psnr_value), ' dB']);

% Function definition of RLS adaptive filter
function output = rls_adaptive_filter(noisy_audio, desired_audio, lambda, delta, filter_size)
    N = length(noisy_audio);
    output = zeros(N, 1);
    error = zeros(N, 1);
    w = zeros(filter_size, 1); % Initialize filter weights
    P = eye(filter_size) / delta; % Initialize inverse correlation matrix

    for n = filter_size:N
        % Extract the current window of the noisy signal
        x = noisy_audio(n:-1:n-filter_size+1);

        % Calculate gain vector
        k = (P * x) / (lambda + x' * P * x);

        % Calculate the filter output
        y = w' * x;

        % Calculate the error signal
        e = desired_audio(n) - y;
        error(n) = e;

        % Update the filter weights
        w = w + k * e;

        % Update the inverse correlation matrix
        P = (P - k * x' * P) / lambda;
        
        % Store the output signal
        output(n) = y;
    end
end