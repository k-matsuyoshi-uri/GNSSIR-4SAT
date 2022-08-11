function rmse = get_RMSE(stdev_values)

% function rmse = get_RMSE(stdev_values)
% Root-Mean-Square Error calculation of the average standard deviation
% value given a vector of standard deviations. 
    % Kayli Matsuyoshi

rmse = sqrt(sum(stdev_values.^2)/length(stdev_values));