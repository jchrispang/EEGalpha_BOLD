function data_filt = filter_data_1D(data, interval)
% filters data to remove NaNs and outliers (set by the percentile interval)

if nargin < 2 || isempty(interval)
    interval = [0, 100];
end

% make a copy of the original data
data_filt = data;

% Filter 1: remove the NaNs
data_filt(isnan(data)) = [];

% Filter 2: remove the outliers
prctile_x = prctile(data_filt, interval);

low_x = data_filt < prctile_x(1);
high_x = data_filt > prctile_x(2);
data_filt(((low_x + high_x) > 0)) = [];