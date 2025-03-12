close all;
clear;

% load image
img = imread('source.jpg');

% show the image
figure;
imshow(img);
title('Original Image');

% add Gaussian noise ~ N(0,0.01)
noisy_img = imnoise(img, 'gaussian', 0, 0.01);

% show the image with noise
figure;
imshow(noisy_img);
title('Image with noise');

% set parameters of RLS
lambda = 0.99; % forgetting factor, determines how fast the old data are forgotten and makes the filter suitable for a nonstationary environment 
delta = 0.5; % regularization parameter
filter_size = 3; % order of filter (M)

tic;
% apply RLS adaptive filter to each channel (red, blue, green)
filtered_img = zeros(size(noisy_img));
error_img = zeros(size(noisy_img));
for channel = 1:3
    [filtered_img(:,:,channel), error_img(:,:,channel)] = rls_adaptive_filter(double(noisy_img(:,:,channel)), double(img(:,:,channel)), lambda, delta, filter_size);
end
elapsedtime = toc;
elapsedtimeMS = elapsedtime * 1000;

% show the filtered image
figure;
imshow(uint8(filtered_img));
title('Filtered Image using RLS Adaptive Filter');

% show the error image
figure;
imshow(uint8(abs(error_img * 5)));
title('Error Image');

disp(['Execution time of RLS: ', num2str(elapsedtimeMS), ' ms']);

% calculate PSNR
psnr_R = psnr(filtered_img(:,:,1), img(:,:,1));
psnr_G = psnr(filtered_img(:,:,2), img(:,:,2));
psnr_B = psnr(filtered_img(:,:,3), img(:,:,3));
psnr_total = (psnr_R + psnr_G + psnr_B) / 3;
disp(['PSNR (R channel): ', num2str(psnr_R),' dB']);
disp(['PSNR (G channel): ', num2str(psnr_G),' dB']);
disp(['PSNR (B channel): ', num2str(psnr_B),' dB']);
disp(['Overall PSNR: ', num2str(psnr_total),' dB']);

% compute MSE
mse_R = immse(uint8(filtered_img(:,:,1)), img(:,:,1));
mse_G = immse(uint8(filtered_img(:,:,2)), img(:,:,2));
mse_B = immse(uint8(filtered_img(:,:,3)), img(:,:,3));
mse_value = (mse_R + mse_G + mse_B) / 3;
disp(['MSE of RLS: ', num2str(mse_value)]);

% function definition of RLS adaptive filter 
function [output_img, error] = rls_adaptive_filter(noisy_img, desired_img, lambda, delta, filter_size)
    [rows, cols] = size(noisy_img);
    padded_img = padarray(noisy_img, [floor(filter_size/2), floor(filter_size/2)], 'symmetric');
    output_img = zeros(rows, cols);
    error = zeros(rows, cols);
    
    % initialize RLS parameters (p.361)
    w_vec = zeros(filter_size^2, 1); % Initialize filter weights
    P_vec = (1 / delta) * eye(filter_size^2); % Initialize P matrix

    for i = 1:rows
        for j = 1:cols
            % small region
            local_region_output = padded_img(i:i+filter_size-1, j:j+filter_size-1);
            local_region_output_vector = local_region_output(:); % Convert to 1D array

            % recursive (p.361)
            % calculate output
            y = w_vec' * local_region_output_vector;
            output_img(i, j) = y;

            % calculate error
            xi = desired_img(i, j) - y;
            error(i, j) = xi;
            
            k_vec = (P_vec * local_region_output_vector) / (lambda + local_region_output_vector' * P_vec * local_region_output_vector);
            w_vec = w_vec + k_vec * xi;
            P_vec = (P_vec - k_vec * local_region_output_vector' * P_vec) / lambda;
        end
    end
    
    % Final weighted matrix (for testing)
    w_matrix = reshape(w_vec, [filter_size, filter_size]);
    disp('Final weighted matrix:');
    disp(w_matrix);
end

% function of calculate PSNR
function psnr = psnr(filtered_channel, original_channel)
    mse = mean((double(filtered_channel) - double(original_channel)).^2, 'all');
    max_pixel = 255.0;
    psnr = 20 * log10(max_pixel / sqrt(mse));
end