function [data_x_filt, data_y_filt] = filter_data(data_x, data_y, interval)
% filters data_x and data_y to remove NaNs and outliers (set by the
% percentile interval)

if nargin < 3 || isempty(interval)
    interval = [0, 100];
end

% make a copy of the original data
data_x_filt = data_x;
data_y_filt = data_y;

% Filter 1: remove the NaNs
x_nan = isnan(data_x);
y_nan = isnan(data_y);
data_x_filt((x_nan+y_nan) > 0) = [];
data_y_filt((x_nan+y_nan) > 0) = [];

% Filter 2: remove the outliers
prctile_x = prctile(data_x_filt, interval);
prctile_y = prctile(data_y_filt, interval);

low_x = data_x_filt < prctile_x(1);
high_x = data_x_filt > prctile_x(2);
low_y = data_y_filt < prctile_y(1);
high_y = data_y_filt > prctile_y(2);
data_x_filt(((low_x + high_x + low_y + high_y) > 0)) = [];
data_y_filt(((low_x + high_x + low_y + high_y) > 0)) = [];
