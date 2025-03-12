# Adaptive Filtering for Image and Audio Processing

## Overview
This repository contains MATLAB implementations of adaptive filtering algorithms, specifically the Normalized Least Mean Squares (NLMS) and Recursive Least Squares (RLS) algorithms. These filters are applied to both image and audio processing tasks to reduce noise and improve signal quality.

## Files

### Image Processing
- **NLMS_ImageProcessing.m**: Applies NLMS adaptive filtering to images with Gaussian noise.
- **RLS_ImageProcessing.m**: Applies RLS adaptive filtering to images with Gaussian noise.

### Audio Processing
- **NLMS_VoiceProcessing.m**: Applies NLMS adaptive filtering to audio signals with Gaussian noise.
- **RLS_VoiceProcessing.m**: Applies RLS adaptive filtering to audio signals with Gaussian noise.

## Features
- **Image Denoising**: Reduces Gaussian noise in images using NLMS and RLS adaptive filters.
- **Audio Denoising**: Enhances speech signals by removing added Gaussian noise.
- **Performance Evaluation**: Computes Peak Signal-to-Noise Ratio (PSNR) and Mean Squared Error (MSE) for quality assessment.
- **Execution Time Measurement**: Measures the computational time of the filtering process.

## Installation
1. Ensure you have MATLAB installed.
2. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/adaptive-filtering.git
   ```
3. Place the source image `source.jpg` and audio file `晚安大小姐 (cut).wav` in the same directory.

## Usage
### Image Processing
Run the MATLAB scripts to apply adaptive filtering to images:
```matlab
NLMS_ImageProcessing;
RLS_ImageProcessing;
```
- The output includes filtered images and error maps.
- PSNR and MSE metrics are displayed.

### Audio Processing
Run the MATLAB scripts to apply adaptive filtering to audio signals:
```matlab
NLMS_VoiceProcessing;
RLS_VoiceProcessing;
```
- The output includes filtered audio files (`filtered_audio_NLMS.wav` and `filtered_audio_RLS.wav`).
- PSNR and MSE metrics are displayed.

## Results
The scripts generate filtered images and audio files. The effectiveness of filtering is evaluated using:
- **PSNR**: Higher values indicate better noise reduction.
- **MSE**: Lower values indicate better accuracy in signal restoration.

## Future Improvements
- Implement additional adaptive filtering algorithms.
- Optimize filter parameters for better noise reduction.
- Extend support for real-time processing.

## Acknowledgments
Special thanks to researchers and developers contributing to adaptive filtering techniques.

