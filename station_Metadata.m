% input station metadata to generate a .mat file
% Kayli Matsuyoshi

%% setup 
close all;
clear;

%% Instructions: change variables as needed 
% GNSS station name
GNSS_station = 'N001';

% TG MSL csv data - water height from tidal datum epoch MSL
TG_MSL_data = 'N001_TG_MSL.csv';
TG_station = '8452660'; % NOAA Tide Gauge Water Level Station
TG_epoch = '1983-2001'; % Tidal Datum Analysis Period for MSL from station datum

%%%% PRIMARY STATION METRICS %%%%

% height of GNSS Antenna Phase Center (APC) above Antenna Reference Point (ARP) (m)
    % see GNSS Antenna Spec Sheet for APC to ARP measurement
        % spec sheet: 61.66 mm, measurement: 65.19 mm
APC_to_ARP = 0.0652; 

% height of GNSS ARP above Revised Local Reference (RLR) (m)
    % see GNSS Station Benchmark from SONEL
        % https://www.sonel.org/spip.php?page=nivellement&idStation=2361&campagne=594
ARP_to_RLR = 12.1550; 

% height of Tide Gauge Station Datum (STND) above RLR (m)
    % see Tide Gauge Station Datum distance from RLR from PSMSL 
        % https://www.psmsl.org/data/obtaining/rlr.diagrams/351.php
STND_to_RLR = 6.0190; 

% height of tidal datum epoch MSL above STND (m)
    % see NOAA Tides/Water Level Datums 
        % https://tidesandcurrents.noaa.gov/datums.html?datum=STND&units=1&epoch=0&id=8452660&name=Newport&state=RI
MSL_to_STND = 1.1060;

%%%% OTHER STATION METRICS %%%%
% to be used if they are more reliable than heights from RLR

% height of GNSS APC above MWWL Sensor Level Point (m)
    % see Tape Down Spreadsheet for Station
%APC_to_MWWL = 1.1492;

% height of MWWL Sensor above STND (m) 
    % NOAA-derived measurement
%MWWL_to_STND = 5.052;

% height of ARP above STND (m)
    % NOAA-derived measurement
%ARP_to_STND = 6.136; 

% custom colors to plot for each satellite
satColors = [
    10/255,191/255,0/255;   % gps
    183/255,40/255, 213/255;   % glonass
    189/255,13/255,0/255; % galileo
    0/255,135/255,237/255 % beidou
    ];  

save([GNSS_station '_Station_Metadata.mat']);