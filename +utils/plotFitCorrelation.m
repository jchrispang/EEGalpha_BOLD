function [H, fit, r, p] = plotFitCorrelation(data_x, data_y, ...
                                        label_x, label_y, textPosition, axisPosition)

if nargin < 6 || isempty(axisPosition)
    axisPosition = get(gca, 'Position');
end

if nargin < 5 || isempty(textPosition)
    textPosition = 'NorthEast';
end

if nargin < 4 || isempty(label_y)
    label_y = 'y';
end

if nargin < 3 || isempty(label_x)
    label_x = 'x';
end

% sort the x-axis data in ascending order
data_x_sorted = sort(data_x);

% perform line fit
fit = polyfit(data_x, data_y, 1);

% get linear regression from fit results
regress = polyval(fit, data_x_sorted);

% calculate correlation coefficient and p-value
[r, p] = corrcoef(data_x, data_y);
r = r(1,2);
p = p(1,2);

% calculating figure text offset to avoid overlap
if strcmpi(textPosition, 'NorthEast')
    offset = [axisPosition(3)*0.3, axisPosition(4)*0.7, -axisPosition(3)*0.3, -axisPosition(4)*0.7];
    textpos = axisPosition + offset;
elseif strcmpi(textPosition, 'SouthEast')
    offset = [axisPosition(3)*0.3, axisPosition(4)*0.15, -axisPosition(3)*0.3, -axisPosition(4)*0.7];
    textpos = axisPosition + offset;
elseif strcmpi(textPosition, 'Best')
    if (abs(r)>=0.035 && abs(r)<=0.85) || isnan(r) 
        offset = [axisPosition(3)*0.3, axisPosition(4)*0.7, -axisPosition(3)*0.3, -axisPosition(4)*0.7];
        textpos = axisPosition + offset;
    else
        offset = [axisPosition(3)*0.3, axisPosition(4)*0.03, -axisPosition(3)*0.3, -axisPosition(4)*0.7];
        textpos = axisPosition + offset;
    end
end

% making text annotation strings
if abs(fit(1)) < 0.01
    text1 = sprintf('slope = %1.e', fit(1));
else
    text1 = sprintf('slope = %0.2f', fit(1));
end
text2 = sprintf('r = %0.2f', r);
if p~=0
    text3 = sprintf('p = %1.e', p);
else
    text3 = sprintf('p = %1.f', p);
end

% plotting the figure
plot(data_x, data_y, 'k.', 'Markersize', 5)
hold on;
if ~isnan(r)
    plot(data_x_sorted, regress, 'k-', 'linewidth', 2.5)
end
hold off;
% utils.textbp({text1; text2; text3}, 'HorizontalAlignment', 'left', 'FontSize', 13)
annotation('textbox', textpos, ...
    'Linestyle', 'none', 'FontSize', 15, ...
    'HorizontalAlignment', 'right', 'string', {text2}, 'Color', 'k',...
    'Fontweight', 'b')
% annotation('textbox', textpos, ...
%     'Linestyle', 'none', 'FontSize', 15, ...
%     'HorizontalAlignment', 'right', 'string', {text1; text2; text3}, 'Color', 'b')
set(gca, 'fontSize', 13, 'xlim', [min(data_x)-min(data_x)*0.05, max(data_x)+max(data_x)*0.05], ...
    'ylim', [min(data_y)-min(data_y)*0.1, max(data_y)+max(data_y)*0.1], ...
    'ticklength', [0.02, 0.02]);
ylabel(label_y, 'fontsize', 15, 'interpreter', 'latex')
xlabel(label_x, 'fontsize', 15, 'interpreter', 'latex')

H = get(gca);