function [evals, evaluations] = compute_Evaluations(num_constellations, SL_GNSSIR_MED, SL_GNSSIR_MEAN, SL_stdev, TG_daily_mean, SL_GNSSIR, SL_TG, sat_index)

% function evals = compute_Evaluations(SL_GNSSIR_MED, SL_GNSSIR_MEAN, SL_stdev, TG_daily_mean, SL_GNSSIR, SL_TG, sat_index)
% Compute evaluation metrics (correlation, RMSD, RMSE) for GNSS-IR and Tide
% Gauge Sea Level comparison
% correlation, RMSD, and RMSE for daily measurements
% correlation and RMSD for all satellites and satellite combinations
    % Kayli Matsuyoshi

evals = struct;

% daily evaluation
evals.corr_med_daily = corr(SL_GNSSIR_MED, TG_daily_mean);
evals.corr_mean_daily = corr(SL_GNSSIR_MEAN, TG_daily_mean);
evals.RMSD_daily = get_RMSD(SL_GNSSIR_MEAN, TG_daily_mean);
evals.RMSE_daily = get_RMSE(SL_stdev);

% all constellations evalation
evals.corr = corr(SL_GNSSIR, SL_TG);
evals.RMSD = get_RMSD(SL_GNSSIR, SL_TG);

% gps evaluation
evals.corr_gps = corr(SL_GNSSIR(sat_index.gps), SL_TG(sat_index.gps));
evals.RMSD_gps = get_RMSD(SL_GNSSIR(sat_index.gps), SL_TG(sat_index.gps));
% glonass evaluation
evals.corr_glo = corr(SL_GNSSIR(sat_index.glonass), SL_TG(sat_index.glonass));
evals.RMSD_glo = get_RMSD(SL_GNSSIR(sat_index.glonass), SL_TG(sat_index.glonass));

if (num_constellations == 4)
    % galileo evaluation
    evals.corr_gal = corr(SL_GNSSIR(sat_index.galileo), SL_TG(sat_index.galileo));
    evals.RMSD_gal = get_RMSD(SL_GNSSIR(sat_index.galileo), SL_TG(sat_index.galileo));
    % beidou evaluation
    evals.corr_bei = corr(SL_GNSSIR(sat_index.beidou), SL_TG(sat_index.beidou));
    evals.RMSD_bei = get_RMSD(SL_GNSSIR(sat_index.beidou), SL_TG(sat_index.beidou));
end

% gps+glonass evaluation
evals.corr_gps_glo = corr(SL_GNSSIR(sat_index.gps_glo), SL_TG(sat_index.gps_glo));
evals.RMSD_gps_glo = get_RMSD(SL_GNSSIR(sat_index.gps_glo), SL_TG(sat_index.gps_glo));

if (num_constellations == 4)
    % gps+galileo evaluation
    evals.corr_gps_gal = corr(SL_GNSSIR(sat_index.gps_gal), SL_TG(sat_index.gps_gal));
    evals.RMSD_gps_gal = get_RMSD(SL_GNSSIR(sat_index.gps_gal), SL_TG(sat_index.gps_gal));
    % gps+beidou evaluation
    evals.corr_gps_bei = corr(SL_GNSSIR(sat_index.gps_bei), SL_TG(sat_index.gps_bei));
    evals.RMSD_gps_bei = get_RMSD(SL_GNSSIR(sat_index.gps_bei), SL_TG(sat_index.gps_bei));
    % glonass+galileo evaluation
    evals.corr_glo_gal = corr(SL_GNSSIR(sat_index.glo_gal), SL_TG(sat_index.glo_gal));
    evals.RMSD_glo_gal = get_RMSD(SL_GNSSIR(sat_index.glo_gal), SL_TG(sat_index.glo_gal));
    % glonass+beidou evaluation
    evals.corr_glo_bei = corr(SL_GNSSIR(sat_index.glo_bei), SL_TG(sat_index.glo_bei));
    evals.RMSD_glo_bei = get_RMSD(SL_GNSSIR(sat_index.glo_bei), SL_TG(sat_index.glo_bei));
    % galileo+beidou evaluation
    evals.corr_gal_bei = corr(SL_GNSSIR(sat_index.gal_bei), SL_TG(sat_index.gal_bei));
    evals.RMSD_gal_bei = get_RMSD(SL_GNSSIR(sat_index.gal_bei), SL_TG(sat_index.gal_bei));

    % gps+glonass+galileo evaluation
    evals.corr_gps_glo_gal = corr(SL_GNSSIR(sat_index.gps_glo_gal), SL_TG(sat_index.gps_glo_gal));
    evals.RMSD_gps_glo_gal = get_RMSD(SL_GNSSIR(sat_index.gps_glo_gal), SL_TG(sat_index.gps_glo_gal));
    % gps+glonass+beidou evaluation
    evals.corr_gps_glo_bei = corr(SL_GNSSIR(sat_index.gps_glo_bei), SL_TG(sat_index.gps_glo_bei));
    evals.RMSD_gps_glo_bei = get_RMSD(SL_GNSSIR(sat_index.gps_glo_bei), SL_TG(sat_index.gps_glo_bei));
    % gps+galileo+beidou evaluation
    evals.corr_gps_gal_bei = corr(SL_GNSSIR(sat_index.gps_gal_bei), SL_TG(sat_index.gps_gal_bei));
    evals.RMSD_gps_gal_bei = get_RMSD(SL_GNSSIR(sat_index.gps_gal_bei), SL_TG(sat_index.gps_gal_bei));
    % glonass+galileo+beidou evaluation
    evals.corr_glo_gal_bei = corr(SL_GNSSIR(sat_index.glo_gal_bei), SL_TG(sat_index.glo_gal_bei));
    evals.RMSD_glo_gal_bei = get_RMSD(SL_GNSSIR(sat_index.glo_gal_bei), SL_TG(sat_index.glo_gal_bei));
end

if (num_constellations == 2)
    evaluations = [evals.corr_gps, evals.RMSD_gps; % 1
        evals.corr_glo, evals.RMSD_glo; % 2
        evals.corr_gps_glo, evals.RMSD_gps_glo; % 3
        evals.corr,  evals.RMSD % 4
        evals.corr_mean_daily,  evals.RMSD_daily % 5
        ]; 
else
    evaluations = [evals.corr_gps, evals.RMSD_gps; % 1
        evals.corr_glo, evals.RMSD_glo; % 2
        evals.corr_gal, evals.RMSD_gal; % 3
        evals.corr_bei, evals.RMSD_bei; % 4
        evals.corr_gps_glo, evals.RMSD_gps_glo; % 5
        evals.corr_gps_gal, evals.RMSD_gps_gal; % 6
        evals.corr_gps_bei, evals.RMSD_gps_bei; % 7
        evals.corr_gps_glo, evals.RMSD_gps_gal; % 8
        evals.corr_gps_glo, evals.RMSD_gps_bei; % 9
        evals.corr_gps_gal, evals.RMSD_gps_bei;% 10
        evals.corr_gps_glo_gal, evals.RMSD_gps_glo_gal; % 11
        evals.corr_gps_glo_bei, evals.RMSD_gps_glo_bei; % 12
        evals.corr_gps_gal_bei, evals.RMSD_gps_gal_bei; % 13
        evals.corr_glo_gal_bei, evals.RMSD_glo_gal_bei; % 14
        evals.corr,  evals.RMSD % 15
        evals.corr_mean_daily,  evals.RMSD_daily % 16
        ]; 
end


