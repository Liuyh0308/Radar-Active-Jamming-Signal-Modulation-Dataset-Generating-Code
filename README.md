# Radar Active Jamming Signal Modulation Dataset Generating Code
## Code description
- This is a simulation **dataset generation** code for radar active jamming signals written by the author himself, please open this code with postfix `.m` in **Matlab** and then just run it directly.
## Usage
1. This code would build `.mat` file and `.png` file in your setted directory path, so please notice carefully that the *path settings* on your PC is correct.
2. This code can generate and save three kinds of eigenmode data of **11** kinds of active jamming signals in **time domain, frequency domain and time-frequency domain**, which not only have 1D sequence modality but also 2D image modality.
  > These two different modalities can be used for training an *auto radar interference detector, radar interference identifier, an end2end radar denoising system*.
3. You can also **adjust various signal parameters** in the code file, or **add multipath fading or Doppler shift** as you wish to simulate different scenarios.
  > Various signal parameters include: SNR, carrier frequency, center frequency, JSR, delay time, and parameters related to specific interference.
    > Specific interference parameters, such as the number of samples and forwardings of ISRJ, can be customized.





