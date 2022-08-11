# GNSSIR-4SAT
Evaluating GNSS-IR method's performance in measuring sea level using 4 satellite constellations

by Kayli Matsuyoshi

Introduction \
I wrote this code while working on my project titled "Improving the Global Navigation Satellite System-Interferometric Reflectometry (GNSS-IR) Method to Measure Sea Level Rise" (2022). \
It builds upon Kristine Larson's GNSS-IR method to measure sea level, and Ben Watzak's work with the NOAA Newport GNSS Station. This project aimed to improve the GNSS-IR method by using Signal-to-Noise Ratio (SNR) data from four different satellite constellations (GPS, GLONASS, Galileo, and BeiDou). The method's temporal resolution increases with data from more constellations, and when a frequency-dependent antenna phase center (APC) offset is added to correct the sea level measurements, the method's accuracy is also improved. \
With this code, one should be able to analyze raw SNR data from the four satellite constellations, and compare GNSS-IR sea level measurements with a co-located Tide Gauge sea level record.

Prior to using this code
- SNR data should be retrieved using Kristine Larson's gnssSNR code found in the Github Repository https://github.com/kristinemlarson/gnssSNR. 
- Quality Control satellite elevation and azimuth ranges should be identified using Kristine Larson's FresnelMaps code found in the Github Repository https://github.com/kristinemlarson/FresnelMaps. 
- Another tool that summarizes this process is Kristine Larson's online GNSS-IR Web App found at https://gnss-reflections.org/fancy6. \
Refer to the study Software tools for GNSS interferometric reflectometry (GNSS-IR) (2018) for more information about these resources.

Contents
- view_SNR_stats.m \
    View raw data from SNR files. \
    Requires: SNR data
- LSP_KM.m\
    Modified version of Ben Watzack's Lomb-Scargle Periodogram processing of SNR data. The script's original authors were Kristine Larson and Carolyn Roesler. The modifications made to this script aimed to speed up the LSP stage of the GNSS-IR method. \
    Requires: SNR data\
    Outputs: SNR_stats.mat files, RH S1.txt files
- snr_files_KM.m\
    Helper function for LSP_KM that extracts the SNR file and creates an output txt file with the user's input. 
- station_Metadata.m\
    Creates a .mat file using user-input parameters that summarizes distances between station benchmarks and associated files.\
    Outputs: Station_Metadata.mat file
- get_RMSD.m\
    Helper function for GNSSIR_vs_TG_SL scripts to evaluate the root-mean-square deviation. 
- get_RMSE.m\
    Helper function for GNSSIR_vs_TG_SL scripts to evaluate the root-mean-square error. 
- GNSSIR_vs_TG_SL_4SAT.m\
    Evaluates the performance of the GNSS-IR method of measuring sea level using 4 satellite constellations (GPS, GLONASS, Galileo, BeiDou). Reflector height measurements generated from the LSP_KM are converted to the same reference frame as a co-located tide gauge's sea level from MSL measurements. The script also calculates appropriate APC offsets. \
    Requires: RH S1.txt files, SNR statistics .mat files, Station_Metadata.mat file\
    Outputs: evaluations.mat file
- GNSSIR_vs_TG_SL_2SAT.m\
    Evaluates the performance of the GNSS-IR method of measuring sea level using 2 satellite constellations (GPS and GLONASS). A frequency-dependent antenna phase center offset for the GLONASS constellation must be loaded or written. \
    Requires: RH measurement .txt files, SNR statistics .mat files, Station_Metadata.mat file\

Additional Files Included
- N001_Station_Metadata.mat
- N001_evaluations.mat - contains empiracal model to find APC offsets
- N001_TG_MSL.csv
- N300_Station_Metadata.mat
- N300_evaluations.mat
- N300_TG_MSL.csv

Additional Scripts Not Included
- get_ofac_hifac.m
- get_waveL.m
- lomb.m
- median_RH.m
- open_filesX.m
- output_header.m
- peak2noise.m
- peakcheck.m
- plot_labels.m
- plot_peakcheck.m
- quickly_QC.m\
(Can be found in Kristine Larson's repository https://github.com/kristinemlarson/gnssIR_matlab_v3)

General Instructions
1. Obtain SNR files containing S1 signal strength data 
2. View SNR data using view_SNR_stats.m
3. Generate Reflector Height (RH) measurements from SNR data using LSP_KM.m
4. Input Station Metadata with station_Metadata.m
5. Evaluate the GNSS-IR method's performance at a GNSS+Tide Gauge station using GNSSIR_vs_TG_SL_4SAT.m

Acknowledgements\
For the foundation of the Matlab scripts used in this project, I would like to credit\
Kristine Larson and Ben Watzak. \
For guidance and expertise with this project, I would like to thank \
Kristine Larson (University of Colorado Boulder)\
Ben Watzak (Previous Undergraduate Research Fellow)\
Meng (Matt) Wei (URI Graduate School of Oceanography)\
Eric Breuer and Robert Heitsenrether (NOAA Tides and Currents)\
For the opportunity and funding for this project, I would like to thank \
the University of Rhode Island, Graduate School of Oceanography\
and the National Science Foundation.
