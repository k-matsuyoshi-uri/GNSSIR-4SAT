function [ filename, year, month, dom, station, outputfile,freqtype,nofile] = snr_files_KM(SNR_directory, OO)
% function [ filename, year, dom, station, outputfile,freqtype,nofile ] = sample_files( )
% user for the tutorial either uses a provided sample file or inputs
% information for their own files.
% outputs are filename (SNR)
% year 
% dom
% outputfile is a txt file to save the LSP outputs
% freqtype is 1, S1
% 
% Author: Kristine M. Larson
% Editors for Newport Data and other applications: Ben Watzak (2019), Kayli Matsuyoshi (2022)

nofile = false;

% input your own site
station = input('station: ','s');   

year = input('4-digit year: ');
% 2022
if (year - 2000 < 0)
    cyr = sprintf('%04d', year + 2000); %make sure it's a 4-digit year
else
    cyr = num2str(year);
end
month = input('2-digit month: ');
% 03 to 06
cmo = sprintf('%02d', month ); %make sure it's a 2-digit month
dom = input('2-digit day of month: ');
cday = sprintf('%02d', dom ); %make sure it's a 2-digit day

% You can add a directory structure if you prefer
% This assumes snr99 files - but you can modify to allow different ones 
filename = [station '_' cyr '_' cmo cday '.snr' OO]; 
filename = [SNR_directory filename]; 

if ~exist(filename, 'file')
    disp('SNR file does not exist. Exiting.')
    % return dummy values
    outputfile = ''; freqtype = 0; year = 0; dom = 0; station = '';
    nofile = true;
    return
else
    fprintf(1,'SNR data should be in : %s \n', filename);
    freqtype = 1;
    outputfile = ['RH_S1/' station '_' cyr '_' cmo cday '_S' num2str(freqtype) '.txt'];
    fprintf(1,'LSP Output will go to : %s \n', outputfile);
end

end
