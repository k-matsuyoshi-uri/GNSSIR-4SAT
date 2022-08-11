% Authors: Kristine M. Larson & Carolyn J. Roesler
% http://kristinelarson.net
% April 16, 2018

% Contributors/Editors: Ben Watzak (2019), Kayli Matsuyoshi (2022)

% Revised notes and script from Ben's version (Kayli)

%% Introduction 

% The GPS Solutions paper that accompanies this code:
% GNSS Interferometric Reflectometry Software Tools
% Carolyn Roesler and Kristine M. Larson, GPS Solutions,TBD

% For cryosphere applications, please also cite: Larson, K.M., J. Wahr, and P. Kuipers Munneke, Constraints on Snow Accumulation and Firn Density in Greenland 
% Using GPS Receivers, J. Glaciology, Vol. 61, No. 225, doi:10.3189/2015JoG14J130, 2015

% This code should not be seen as the final word on picking proper inputs for your own GNSS-IR experiment. 
% Picking the proper azimuths and elevation angles is site dependent.

% The main outputs of this code are REFLECTOR HEIGHTS, or RH, in meters. 
% RH is the vertical distance between the GPS phase center and the ice/snow/water surface.  
% If you have field measurements of GPS antenna height, please remember that the GPS phase center is not 
% in the same place as the ground plane that most people use for field measurements of antenna height.

% Please also keep in mind that you can always pick a peak value for a periodogram, but that doesn't mean it is significant.
% It is your job to determine if your peak is meaningful. 
% I have used different ways to do this.  Possible Quality control measures: 
%  1. amplitude of the LSP peak.
%  2. amptliude of LSP peak relative to a noise level
%  3. width of LSP peak (but this depends on sample interval, so that is not implemented here).
%     do not use it here).
%  4. you don't want double peaks
%  5. you (generally) don't want a satellite track that both rises and sets. 
% I have currently set a default of peak being 3.5 times larger than the noise, but this is just to get you started.  
% This assumes you have picked a relevant interval to calculate the noise over.  

%% Editor Notes

% I modified this script by Kristine Larson to speed up the LSP model step
% of the GNSS-IR Method. Ben's edits allowed for the processing of GNSS
% data from the Newport Station N001. I aimed to apply this script to any
% specified station and data. The SNR data taken will be signal strength S1
% and GPS, GLONASS, Galileo, and BeiDou satellites are permitted. 

% The code is meant to process SNR files day by day. To most efficently use
% this code, have an SNR directory containing SNR files named with the
% format SSSS_YYYY_MMDD.snrOO
% SSSS: four character station name 
% YYYY: four digit year
% MM: two digit month
% DD: two digit day of month
% OO: SNR processing option (https://github.com/kristinemlarson/gnssSNR)

% Kayli Matsuyoshi

%% Specify Parameters

% Directory containing SNR files
SNR_directory = 'N001_SNR/'; 

% SNR processing option (for example: .snr99)
OO = '99'; 

% Satellite Elevation Angle Range (degrees)
    % My default elevation angles are usually 5 and 25. (Larson) 
    % Often station operators use an elevation mask - which is not needed, 
    % as a data analyst also sets an elevation mask. 
    % In our opinion, the mask should be applied at the analysis stage, 
    % not at the data collection stage.
    
    % !! Will need to look into the newport data to see if mask is already
    % being applied at collection. (Ben)
    
emin = 5; emax = 16;

% Satellite Azimuth Angle Range (degrees)
    % az_max must be larger than az_min
    % both values must be between 0 and 360
az_min = 15;    az_max = 145.1; % N001
% az_min = 80;    az_max = 240; % N300

% Elevation Angle Difference (degrees)
    % ediff is the elevation angle different for a given track. (Larson had
    % it set to 10) 
    % This variable is a QC metric because you don't want a bunch of tiny
    % arcs as these periodograms can be very unreliable. 
    % If you wanted to use data between elevation angles of 5 and 10
    % degrees, for example, you would have to change it. (Larson)
    
    % Looks like this is the minimum required difference between observed
    % min and max elevation angles, i.e. (maxObsE-minObsE) must be > ediff
    % (Ben)

ediff = 6 ;
 
% Maximum Possible Refector Height (m)
    % (i.e. exclude reflector heights beyond this value) 
    % (Larson set at 8, Ben set at 6.3)
maxHeight = 6.3;
      
%% Other Parameters

plt_type = 0; % (KM)

% Rule of thumb: you should not think you are correctly estimating
% RH when it is less than 2*lambda, which is 40-50 cm, depending on the
% wavelength (l1 vs l2). 
% Recall: RH is the vertical distance between the GPS phase center and the
% ice/snow surface. 
minRH   = 1.2; % meters... The water level should never, at least within 
% foreseeable future, get up to this. 
% This minRH represents the distance from GNSS phase center to levelling
% collar (MWWL), and they put the MWWL there because water level shouldn't
% ever get that high. 
maxArcTime = 1; % one hour 
% ??What exactly is the above Arctime?? (Ben)

% Mininum number of points. This value could be used as QC
% it is totally arbitrary for now.  - Larson used 25 as her value
minPoints = 15; 

% Desired Precision (m)
    % 5 mm is a reasonable level for now. (Larson)
desiredPrecision = 0.005; 

% Noise Range (m)
frange = [0 6]; % noise range is calculated between 0 and 5 meters.

%% 
%clear 
close all;

% GPS, GLONASS, Galileo, and BeiDou Satellites from SNR Data Column 1 (KM)
satlist = [1:50 100:150 200:250 300:350]; 

% plotting variables
LW = 2 ; % linewidth for quadrant plots
comL = '%'; %Not sure what this is used for?? (Ben)

pvf = 3; % polynomial order used to remove the direct signal. 
    % Ben thinks this should be changed, and I thought it should
    % be 2. But this value worked well for the available data. 
% this can be smaller, especially for short elevation angle ranges.

% the avg_maxRH variable will store a crude median reflector height for a single day/site
avg_maxRH = [];

all_peaks = [];
lsp_amps = [];

% you can turn off all the plots by changing this to false.
plot2screen = true; 

%% Process SNR to Reflector Heights

% pick up file with function snr_files_LM
[filename, year, month, dom, station, outputfile,freqtype,nofile] = snr_files_KM(SNR_directory, OO);
if nofile
  return
end

% load the SNR data into variable x, open output file for LSP results
% nofile is a boolean that knows whether data have been found
[fid,x,nofile] = open_filesX(filename,outputfile,freqtype);

if nofile
  return
end

figure
% get additional QC levels and print them to the screen
[minAmp,pknoiseCrit,frange ] = quicky_QC(freqtype, maxHeight, desiredPrecision, ediff,frange);

% get wavelength factor (lambda/2) and column where data are stored (ic)
[cf,ic] = get_waveL(freqtype); % returns wavelength factor (lambda/2) and column number (using snr format)
% make header for the output txt file
wavelength = cf*2; %this will be used to catch double peaks

output_header( fid );

% checking azimuths in 15 degree bins

az_diff = az_max - az_min;
azrange = 15; 
naz = round(az_diff/azrange); %number of bins

for a=1:naz
    
    % window by these azimuths
    azim1 = (a-1)*azrange + az_min;
    azim2 = azim1 + azrange;

    % window by satellite
    % i.e. run thru every satellite track
    for sat = satlist
        % x contains the SNR data
        % Column 1 satellite number
        % Column 2 elevation angle
        % Column 3 azimuth angle
        % Column 4 second of the day, GPS time (i.e. no leap seconds)
        % Column 7 S1
        % Column 8 S2
        % Column 9 S5
        %i=find(x(:,2) < emax & x(:,2) > emin & x(:,1) == sat & x(:,3) > azim1 & x(:,3) < azim2); 
        i=find(x(:,2) < emax & x(:,2) > emin & x(:,1) == sat & x(:,3) > azim1 & x(:,3) < azim2);% & x(:,4)>14400 & x(:,4) < (28800 )); 

        % i created a single column vector of the row indices where x meets the
        % QC criteria
        % in some cases you have both ascending and descending arcs in one
        % azimuth bin that will fulfill this find statement,  
        % but for cryosphere applications it probably won't hurt you too much
        % (Larson)

        if length(i) > minPoints
            w = x(i,:); %creates matrix of only the SNR data in x that met QC criteria for the current satellite
            elevAngles = w(:,2); % vector of elevation angles in degrees


            data = 10.^(w(:,ic)/20);   % change SNR data from dB-Hz to linear units
            time = w(:,4)/3600;    %3600 sec/hour => time is in now hours instead of sec
            % these are UTC hours  %Coordinated Universal Time (UTC)
            meanUTC = mean(time);
            dt = time(end) - time(1);    % time span of track in hours
            azm = mean(w(:,3));    %average azimuth for a track, in degrees

            % remove direct signal. polyfit value does not need to be as large as this for some arcs.
            pf=polyfit(elevAngles, data, pvf);  
            %pf is a row vector of length pvf containing the polynomial coefficients indescending powers, P(1)*X^N + P(2)*X^(N-1) +...+ P(N)*X + P(N+1).
            %polyfit(X,Y,N) finds the coefficients of a polynomial P(X) of degree N that fits the data Y best in a least-squares sense.
            pv = polyval(pf, elevAngles);
            %Y = polyval(P,X) returns the value of a polynomial P evaluated at X.
            sineE = sind(elevAngles);  % sin(elevation angles) in degrees
            saveSNR = data-pv; %remove the direct signal with a polynomial
            [sortedX,j] = sort(sineE);
            % sort the data (in ascending order) so all tracks are rising
            %[B,I] = sort(A,...) also returns a sort index I which specifies how the elements of A were rearranged to obtain the sorted output B:
            %If A is a vector, then B = A(I). 
            sortedY = saveSNR(j);
            % get the oversampling factor and hifac. see code for more information
            % ofac (oversampling factor) and hifac (high frequency factor) define the LSP frequency grid spacing and maximum frequency.
            [ofac,hifac] = get_ofac_hifac( elevAngles,cf, maxHeight, desiredPrecision);

            % call the lomb scargle code.  Input data have been scaled so that f comes out in units of reflector height (meters)
            [f,p,dd,dd2]=lomb(sortedX/cf, sortedY, ofac,hifac);
            % returned values are arrays of frequencies considered (f), the associated spectral amplitude(p), ...
            % estimated noise significance of the power values (dd), and the 95% confident level amplitude(dd2).
            [ maxRH, maxRHAmp, pknoise ] = peak2noise(f,p,frange);
            % maxRH should be more than 2*lambda, or ~40-50 cm 
            % here I am restricting arcs to be < one hour.  long dt usually means you have a track that goes over midnite -Larson

            maxObsE = max(elevAngles);
            minObsE = min(elevAngles);   

            %finally, we actually plot the LSP's
            if maxRHAmp > minAmp && maxRH > minRH  && dt < maxArcTime && pknoise > pknoiseCrit && (maxObsE-minObsE) > ediff
                fprintf(fid,'%4.0f %2.0f %2.0f %6.2f %6.2f %6.1f %6.0f %3.0f %6.2f %6.2f %6.2f %4.0f %5.2f\n', ...
                year,month,dom,maxRH,maxRHAmp,azm, sat, dt*60, minObsE, maxObsE, pknoise,freqtype,meanUTC);
                if plot2screen
                    % plot all the periodograms on top of each other in gray
                    subplot(2,1,1)
                    plot(f,p,'color',[0.5 0.5 0.5]); hold on; 
                    set(gcf, 'Position', [400, 400, 750, 500]);
                    % freq units on x-axis will be reflector heights (m) 
                    % amplitude units on y-axis will be dBHz  
                end
                avg_maxRH = [avg_maxRH; maxRH];

                %now, look for psd peaks in the satellite track on the LSP

                num_peaks = 1; %because double peaks are not expected
                [peaks, peak_RH] = findpeaks(p,f,'NPeaks',num_peaks,'SortStr','descend'); 
                %returns up to num_peaks peak values of LSP and the RH where they occur
                %record these values in the daily values vectors
                all_peaks = [all_peaks; peak_RH ]; %#ok<*AGROW>
                lsp_amps = [lsp_amps; peaks]; %#ok<*AGROW>

            else 
                fprintf(1,'%s RH %6.2f Amp %6.2f Azm %6.1f Sat %2.0f Tdiff %4.0f Emin %6.2f Emax %6.2f Peak2Noise %6.2f \n', ...
                'Fail QC',maxRH,maxRHAmp,azm, sat, dt*60, minObsE, maxObsE, pknoise);
            end % did you pass QC test loop

        end %do you have enough points loop

    end % satellite loop

end % azimuth loop

if plot2screen && plt_type == 0
    plot_labels( plt_type, station, [az_min az_max], freqtype)
end

if max(avg_maxRH) > (.2445/.1905)* min(avg_maxRH)
    fprintf('%s\n', 'There is likely a double peak and this is likely L2 data')
end

% if you are making one summary plot, figure out the median RH
% and plot it as a magenta vertical line.
median_RH( plt_type, avg_maxRH, plot2screen )
hold on;
if length(all_peaks) < 3
    %If we're dealing with a lack of data
    if freqtype ==2
        avg_RH = min(all_peaks);
        %means there's only one satellite track for L2... definitely can't expect to take an average of one thing
        %it's ok for the purposes of this script. We'll just output it to the screen
        fprintf('%s\n','There was only one satellite track. This is NOT an average then');
    else
        if length(all_peaks) < 2
            avg_RH = all_peaks(1);
            %means there's only one satellite track for L1 or L5... definitely can't expect to take an average of one thing
            %it's ok for the purposes of this script. We'll just output it to the screen
            fprintf('%s\n','There was only one satellite track. This is NOT an average then');
        end
    end     
else
    [gmm, avg_RH, std_dev] = peakcheck(all_peaks, lsp_amps, maxHeight, freqtype, plot2screen);
end

subplot(2,1,1)
hold on;
tx = ['Calculated Average H_R: ' sprintf('%4.2f ',avg_RH) '(m)'];
ylim = get(gca,'Ylim');
text(avg_RH + 0.5, ylim(2)-0.5*diff(ylim), tx,'Color','r')
plot( [avg_RH avg_RH], ylim, 'r--','linewidth',.75)

fclose(fid);

%% Save Figures and Data

date = strrep(erase(erase(filename, append(SNR_directory, station, '_')), append('.snr', OO)), '_', '/');
date = append(extractBefore(date, 8), '/', extractAfter(date, 7));

fig = figure(1);
subplot(2,1,1);
subtitle(date);
set(gca,'fontname','times','FontSize',14);
subplot(2,1,2);
xlabel('H_{R}, Reflector Height (m)');
set(gca,'fontname','times','FontSize',14);
%saveas(fig, append('figures/LSP_S1/N300/',erase(erase(filename, SNR_directory), append('.snr', OO)), '_S1_LSP'), 'png');
med_RH = median(avg_maxRH);
save_file = append('SNR_stats/', erase(erase(filename, SNR_directory), append('.snr', OO)),'_SNR_stats.mat');
RH = struct('med_RH', med_RH, 'avg_RH', avg_RH, 'std_dev_RH', std_dev);
%save(save_file, 'RH', '-append');

%close all;