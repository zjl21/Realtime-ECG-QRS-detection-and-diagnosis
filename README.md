# Realtime ECG QRS Detection and Diagnosis

## Description

Pan-Tompkins Algorithm for ECG Signal Real-Time QRS Detection and Diagnosis.

This is a project for the course "Biomedical Electronics 2 - Instruments" in the third year of Biomedical Engineering at Tsinghua University. The goal is to achieve real-time detection of ECG signals, automatically calculate heart rate, and detect premature beats.

## Usage

There are three modes of operation:

1. **Mode 1:** See the realtime detection process, run
```  python main.py --config config.ini --visualize```
2. **Mode 2:** See the result for each 
3. **Mode 3:**

## Algorithm Overview

The algorithm is based on the Pan-Tompkins method for real-time QRS detection in ECG signals. The basic principles include:

1. **Bandpass Filtering:**
2. **Differentiation:**
3. **Squaring:**
4. **Moving Window Integration:**
5. **Thresholding:**

For more detailed information, please refer to my experiment report, which has been uploaded to the remote repository.

## Results and Illustrations

You can find several illustrations demonstrating the effectiveness of the algorithm in the repository.

## Limitations

- **Edge Data Detection:** The algorithm performs poorly on the initial and final segments of the data.
- **Premature Beat Detection:** The accuracy of detecting premature beats is not high.
