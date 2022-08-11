function rmsd = get_RMSD(vector1, vector2)

% function rmsd = get_RMSD(vector1, vector2)
% Uses Root-Mean-Square Deviation method to calculate average deviation of
% values in vector 1 from vector 2. 
    % Kayli Matsuyoshi

rmsd = sqrt(sum((vector2-vector1).^2)/(length(vector1)));