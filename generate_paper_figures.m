function paper_figures_2ndRound_revision1
save = 0;
% figure2(save)
% figure4(save)
% figure5(save)
% figure6(save)
% figure7(save)
figure8(save)
% figure9(save)
% figure10(save)
% close all

%%
function figure2(save)
% Figure 2: Sample spectra illustrating frequency characteristics
% interpolated low frequencies

load data/fit_output.mat target_P f
bands = model.eeg_bands();

time_slice = 65;%55;

Pexp = target_P(:,time_slice);
fexp = f;
fexp = [0.1, 0.2, 0.3, 0.4, fexp];
Pexp = interp1(f, Pexp, fexp, 'linear', 'extrap');

analysis = model.get_spec_analysis(fexp, Pexp);

delta_freq = linspace(bands(1,1)+min(fexp), bands(1,2), 2);
delta_power = delta_freq.^(analysis.lf_slope)*Pexp(1)/1.3;

ylims = [2e-4, 3];

fig = figure;
pos = get(gcf, 'pos');
delete(fig)
fig = figure('pos', pos.*[1 1 1.2 1]);
loglog(fexp, Pexp, 'k-')
hold on;
plot(analysis.alpha_maxf, analysis.alpha_maxp, 'ko', 'markerfacecolor', 'k', 'markersize', 7)
plot(analysis.alpha_minf, analysis.alpha_minp, 'ko', 'markersize', 7)
% loglog(delta_freq, delta_power, 'k--', 'linewidth', 2)

text(0.55, 1.8, 'low', 'interpreter', 'latex', 'fontsize', 18, 'fontweight', 'b')
text(7.7, 1.8, 'alpha', 'interpreter', 'latex', 'fontsize', 18, 'fontweight', 'b')
% text(0.6, 1.5e-2, 'slope = $-\chi_l$', 'interpreter', 'latex', 'fontsize', 15)
text(analysis.alpha_maxf*1.53, analysis.alpha_maxp*0.4, '$[f_\alpha^{\rm max}, P(f_\alpha^{\rm max})]$', 'interpreter', 'latex', 'fontsize', 15)
text(analysis.alpha_minf*0.17, analysis.alpha_minp*3.5, '$[f_\alpha^{\rm min}, P(f_\alpha^{\rm min})]$', 'interpreter', 'latex', 'fontsize', 15)
annotation('arrow',[0.76 0.703], [0.745 0.793], 'linewidth', 1);
annotation('arrow',[0.58 0.672], [0.645 0.572], 'linewidth', 1);

fill([bands(1,1)+min(fexp) bands(1,2) bands(1,2) bands(1,1)+min(fexp)], [ylims(1) ylims(1) ylims(2) ylims(2)], ...
     [0.95,0.95,0.95],'EdgeColor','k');
fill([bands(3,1) bands(3,2) bands(3,2) bands(3,1)], [ylims(1) ylims(1) ylims(2) ylims(2)], ...
     [0.95,0.95,0.95],'EdgeColor','k');
hold off;
xlabel('$f$ {\rm (Hz)}', 'fontsize', 15, 'interpreter', 'latex')
ylabel('$P(f)$ {\rm (arbitrary units)}', 'fontsize', 15, 'interpreter', 'latex')

set(gca, 'fontsize', 15, 'XLim', [0 45], 'XTick', [0.1,0.5,1,2,3,5,10,20,40], ...
    'YLim', ylims, 'ticklength', [0.02, 0.02], 'layer', 'top')
set(gca,'children',flipud(get(gca,'children')))

if save
    set(fig, 'PaperPositionMode','auto')     %# WYSIWYG
    print(fig, '-painters', '-depsc', 'Manuscript/figures_2ndRound_revision1/Figure2.eps')
end

%%
function figure4(save)
% Figure 4: Sample experimental and fitted spectra for eyes-closd
% Electromyogram (EMG) correction added
% Extrapolated experimental data

% EC_type = '';
EC_type = 'short';

load(['data/spectrogram_data_withMuscle_', EC_type, 'EC.mat']);

data_type1 = 'fit_James_finer';
data_type2 = 'fit_James_muscle_finer';
time_slice = 20;
fexp = w_exp/(2*pi);
ffit = w_finer/(2*pi);

exp = spectrogram.exp(:, time_slice);
exp = interp1(fexp, exp, ffit, 'linear', 'extrap');

fig = figure;
loglog(ffit, exp*1.1, 'k-', 'Linewidth', 0.5)
% loglog(fexp, exp, 'k-', 'Linewidth', 0.5)
hold on;
loglog(ffit, spectrogram.(data_type1)(:, time_slice), 'b:', 'Linewidth', 3)
loglog(ffit, spectrogram.(data_type2)(:, time_slice), 'r-.', 'Linewidth', 3)
hold off;

xlabel('$f$ {\rm (Hz)}', 'fontsize', 15, 'interpreter', 'latex')
ylabel('$P(f)$ {\rm (arbitrary units)}', 'fontsize', 15, 'interpreter', 'latex')
% leg = legend('experiment', 'fit w/o EMG correction', 'fit w/ EMG correction');
% set(leg, 'FontSize', 15, 'Location', 'SouthWest');

set(gca, 'fontsize', 15, 'XLim', [0.1 45], 'XTick', [0.1, 0.5, 1,2,3,5,10,20,40], ...
    'YLim', [1e-4 2], 'ticklength', [0.02, 0.02])

if save
    set(fig, 'PaperPositionMode','auto')     %# WYSIWYG
    print(fig, '-painters', '-depsc', 'Manuscript/figures_2ndRound_revision1/Figure4.eps')
end

%%
function figure5(save)
% Figure 5: BOLD HRF spectral content

params = loadParameters;
f = linspace(0, 10, 10000);
w = 2*pi*f;
params.v_b = 1e-3;
params.Gamma = 1;
Y1 = (params.k2 - params.k3);
Y2 = (params.k1 + params.k2);
P = -params.C_z*(Y1 - params.V_0*Y2);
Q = params.C_z*(Y1*(params.eta + params.tau^(-1)) - ...
Y2*params.C_z*(params.eta - (params.tau^(-1))*(params.beta - 2)) + ...
    (params.D/params.rho_f)*(Y1 - params.V_0*Y2));
R = params.C_z*(params.D/params.rho_f)*(Y1*(params.eta + params.tau^(-1)) - ...
    Y2*params.C_z*(params.eta - (params.tau^(-1))*(params.beta - 2)));

constant = (pi/(2*params.v_b^2*params.Gamma));
num = w.^4*P^2 + w.^2*(Q^2 + 2*P*R) + R^2; 
W = constant*(1./w).*(pi/2 - atan((params.k_z^2*params.v_b^2 - w.^2)./(2*params.Gamma*w)));
L = 1./((-w.^2 + params.kappa^2/4 + params.w_f^2).^2 + w.^2*params.kappa^2);
D = 1./(w.^2 + (params.eta + params.tau^(-1)).^2);
total = num.*W.*L.*D;


fig = figure;
loglog(f, total/max(total), 'k-', 'Linewidth', 2)
xlabel('$f$ {\rm (Hz)}', 'fontsize', 15, 'interpreter', 'latex')
ylabel({'Spectral content of BOLD HRF'; '(arbitrary units)'}, 'fontsize', 15)

set(gca, 'fontsize', 15, 'XLim', [1e-2 0.5], 'XTick', [1e-3, 1e-2, 1e-1, 0.5, 1, 2, 5, 10], ...
     'YLim', [8e-4 2], 'ticklength', [0.02, 0.02])

if save
    set(fig, 'PaperPositionMode','auto')     %# WYSIWYG
    print(fig, '-painters', '-depsc', 'Manuscript/figures_2ndRound_revision1/Figure5.eps')
end

%%
function figure6(save)
% Figure 6: Relationship of delta band power with alpha properties
% 7e-235 3e-134 4e-76

band_attribs = {'int_delta', 'int_alpha', 'alpha_maxp', 'alpha_str'};

band_symbols = {'$P_{\rm low}$', '$P_\alpha$', '$P(f^{\rm{max}}_\alpha)$',...
        '$R_\alpha$'};

load data/pdb_specvalid 

bands = model.eeg_bands();

interval = [5, 95];  % low and high percentile of data to be included
                      % Abeysuriya, et al., J. Neurosci. Meth. 2015 uses
                      % [5, 95]
                      
state = 'ec';

% int_delta versus four alpha variables
y = 'int_delta';
y_ind = utils.cell_find(band_attribs, y);
data_y = analysis_output.(state).(y);
xs = {'int_alpha', 'alpha_maxp', 'alpha_str'};
    
fig = figure('Position', [200, 200, 800, 200]);
for i=1:length(xs)
    h = subplot(1,length(xs),i);
    hpos = get(h, 'Position');
%     delete(h)
%     hpos = hpos+[-0.07,0.09,0.04*length(xs)-0.12,-0.05];
    h = subplot(1,length(xs),i,'Parent',fig,'Position',hpos);
    x = xs{i};
    x_ind = utils.cell_find(band_attribs, x);
    if strcmpi(x, 'alpha_q')
        data_x = analysis_output.(state).('alpha_maxf')/(bands(3,2) - bands(3,1));
    else
        data_x = analysis_output.(state).(x);
    end
    
    [data_x_filt, data_y_filt] = utils.filter_data(data_x, data_y, interval);
    
    if strcmp(x, 'alpha_str')
        x_text = '$\log(R_\alpha)$';
        utils.plotFitCorrelation(log10(data_x_filt), data_y_filt, x_text, band_symbols{y_ind}, 'NorthEast', hpos);
    else
        utils.plotFitCorrelation(data_x_filt, data_y_filt, band_symbols{x_ind}, band_symbols{y_ind}, 'NorthEast', hpos);
    end
    
    annotation('textbox', [0.07, 0.95, 0.1, 0.1], 'string', 'a', 'edgecolor', 'none', ...
        'fontsize', 24, 'fontweight', 'b')
    annotation('textbox', [0.351, 0.95, 0.1, 0.1], 'string', 'b', 'edgecolor', 'none', ...
        'fontsize', 24, 'fontweight', 'b')
    annotation('textbox', [0.632, 0.95, 0.1, 0.1], 'string', 'c', 'edgecolor', 'none', ...
        'fontsize', 24, 'fontweight', 'b')
end

if save
    set(fig, 'PaperPositionMode','auto')     %# WYSIWYG
    print(fig, '-painters', '-depsc', ['Manuscript/figures_2ndRound_revision1/Figure6.eps'])
end

%%
function figure7(save)
% Figure 7: Experimental and fitted spectrogram

% EC_type = '';
EC_type = 'short';

load(['data/spectrogram_data_withMuscle_', EC_type, 'EC.mat']);

finer_limits = dsearchn(w_finer', [w_exp(1); w_exp(end)]);
finer_limits(2) = finer_limits(2) + 2;
fexp = w_exp/(2*pi);
ffit = w_finer(finer_limits(1):finer_limits(2))/(2*pi);

data_exp = spectrogram.exp;
data_fit = spectrogram.fit_James_finer(finer_limits(1):finer_limits(2), :);

% clims = [min([min(data_exp(:)), min(data_fit(:))]), ...
%          max([max(data_exp(:)), max(data_fit(:))])];
clims = [1e-4, ...
         max([max(data_exp(:)), max(data_fit(:))])];
     
fig = figure('Position', [200, 200, 850, 300]);
subplot(1,2,1)
imagesc(1:length(EC_index), fexp, log10(data_exp))
colormap(utils.parula)
colorbar('ytick', -4:0, 'yticklabel', {'10^{-4}', '10^{-3}', '10^{-2}', '10^{-1}', '10^{0}'})
caxis(log10(clims))
xlabel('$t$ {\rm (s)}', 'fontsize', 15, 'interpreter', 'latex')
ylabel('$f$ {\rm (Hz)}', 'fontsize', 15, 'interpreter', 'latex')
set(gca, 'ydir', 'normal', 'Fontsize', 15, 'YLim',[fexp(1) 45], ...
    'Xtick', [1,10:10:100], 'ticklength', [0.02, 0.02])

subplot(1,2,2)
imagesc(1:length(EC_index), ffit, log10(data_fit))
colormap(utils.parula)
colorbar('ytick', -4:0, 'yticklabel', {'10^{-4}', '10^{-3}', '10^{-2}', '10^{-1}', '10^{0}'})
caxis(log10(clims))
xlabel('$t$ {\rm (s)}', 'fontsize', 15, 'interpreter', 'latex')
ylabel('$f$ {\rm (Hz)}', 'fontsize', 15, 'interpreter', 'latex')

annotation('textbox', [0.075, 0.92, 0.1, 0.1], 'string', 'a', 'edgecolor', 'none', ...
        'fontsize', 24, 'fontweight', 'b')
annotation('textbox', [0.515, 0.92, 0.1, 0.1], 'string', 'b', 'edgecolor', 'none', ...
        'fontsize', 24, 'fontweight', 'b')
    
set(gca, 'ydir', 'normal', 'Fontsize', 15, 'YLim',[fexp(1) 45], ...
    'Xtick', [1,10:10:100], 'ticklength', [0.02, 0.02])

if save
    set(fig, 'PaperPositionMode','auto')     %# WYSIWYG
    print(fig, '-painters', '-depsc', 'Manuscript/figures_2ndRound_revision1/Figure7.eps')
end

%%
function figure8(save)
% Figure 8: Time series of spectral properties and feedback loops

% EC_type = '';
EC_type = 'short';

load(['data/spectralProperty_data_', EC_type, 'EC.mat']);
load(['data/phi_data_', EC_type, 'EC.mat'])
load(['data/NFTparams_data_', EC_type, 'EC.mat'])
       
data_type = 'fit_James_finer';

attribs_spectral = {'int_delta', 'int_alpha', 'alpha_maxp', 'alpha_str'};
attribs_spectral_symbols = {'$P_{\rm low}$', '$P_\alpha$',...
                   '$P(f^{\rm{max}}_\alpha)$', '$\log(R_\alpha)$'};
attribs_loop = {'X', 'Y', 'X+Y', 'Z'};
attribs_loop_symbols = {'$X$', '$Y$', '$X + Y$', '$Z$'};
              
fig = figure('Position', [200, 200, 800, 450]);
for j = 1:4
    % spectral properties
    s = subplot(4, 2, 2*j-1, 'ColorOrder',[0,0,0; 1,0,0]);
    if strcmpi(attribs_spectral{j}, 'alpha_str')
        data = log10(spectral_prop.(data_type).(attribs_spectral{j}));
        data_exp = log10(spectral_prop.exp.(attribs_spectral{j}));
    else 
        data = spectral_prop.(data_type).(attribs_spectral{j});
        data_exp = spectral_prop.exp.(attribs_spectral{j});
    end
    yyaxis(s, 'left')
    plot(data, 'k-', 'Linewidth', 2)
    xlabel('$t$ {\rm (s)}', 'fontsize', 15, 'interpreter', 'latex')
    ylabel(attribs_spectral_symbols{j}, 'fontsize', 15, 'interpreter', 'latex')
    
    set(gca, 'fontsize', 15, 'xlim', [1,length(data)], ...
        'ylim', [min(data)*0.95, max(data)*1.05], 'Xtick', [1,10:10:100],  ...
        'ticklength', [0.02, 0.02], ...
        'ytick', round(linspace(min(data)*0.98, max(data)*1.02, 3),2))
    
    [r, ~] = corrcoef(data_exp, data);
    if j==1
        text(2, (max(data)*1.05 - min(data)*0.95)*0.85 + min(data)*0.95, sprintf('r = %0.2f', r(1,2)), ...
                'fontsize', 12, 'horizontalalignment', 'left', 'color', 'k', 'fontweight', 'b')
    else
        text(length(data)-9, (max(data)*1.05 - min(data)*0.95)*0.85 + min(data)*0.95, sprintf('r = %0.2f', r(1,2)), ...
                'fontsize', 12, 'horizontalalignment', 'left', 'color', 'k', 'fontweight', 'b')
    end
    
    yyaxis(s, 'right')
    plot(data_exp, 'r--', 'Linewidth', 1)
            
    set(gca, 'fontsize', 15, 'xlim', [1,length(data)], ...
        'ylim', [min(data_exp)*0.95, max(data_exp)*1.05], ...
        'Xtick', [1,10:10:100],  'ticklength', [0.02, 0.02], ...
        'ytick', round(linspace(min(data_exp)*0.98, max(data_exp)*1.02, 3),2))
    
    
    if j~=4
        set(gca, 'XtickLabel', {})
        xlabel('')
    end
    
    % feedback loops
    subplot(4, 2, 2*j)
    if strcmpi(attribs_loop{j}, 'X+Y')
        data = NFTparams.xyz(1,:) + NFTparams.xyz(2,:);
    elseif strcmpi(attribs_loop{j}, 'Z')
        data = NFTparams.xyz(3,:);
    else
        data = NFTparams.xyz(j,:);
    end
    plot(data, 'k-', 'Linewidth', 2)

    xlabel('$t$ {\rm (s)}', 'fontsize', 15, 'interpreter', 'latex')
    ylabel(attribs_loop_symbols{j}, 'fontsize', 15, 'interpreter', 'latex')
    
    set(gca, 'fontsize', 15, 'Xtick', [1,10:10:100], ...
             'Ylim', [min(data)*0.95, max(data)*1.05],  ...
             'ticklength', [0.02, 0.02], ...
             'ytick', round(linspace(min(data)*0.98, max(data)*1.02, 3),2))

    if j~=4
        set(gca, 'XtickLabel', {})
        xlabel('')
    end
end

annotation('textbox', [0.05, 0.89, 0.1, 0.1], 'string', 'a', 'edgecolor', 'none', ...
    'fontsize', 24, 'fontweight', 'b')
annotation('textbox', [0.05, 0.67, 0.1, 0.1], 'string', 'b', 'edgecolor', 'none', ...
    'fontsize', 24, 'fontweight', 'b')
annotation('textbox', [0.05, 0.45, 0.1, 0.1], 'string', 'c', 'edgecolor', 'none', ...
    'fontsize', 24, 'fontweight', 'b')
annotation('textbox', [0.05, 0.23, 0.1, 0.1], 'string', 'd', 'edgecolor', 'none', ...
    'fontsize', 24, 'fontweight', 'b')
annotation('textbox', [0.49, 0.89, 0.1, 0.1], 'string', 'e', 'edgecolor', 'none', ...
    'fontsize', 24, 'fontweight', 'b')
annotation('textbox', [0.49, 0.67, 0.1, 0.1], 'string', 'f', 'edgecolor', 'none', ...
    'fontsize', 24, 'fontweight', 'b')
annotation('textbox', [0.49, 0.45, 0.1, 0.1], 'string', 'g', 'edgecolor', 'none', ...
    'fontsize', 24, 'fontweight', 'b')
annotation('textbox', [0.49, 0.23, 0.1, 0.1], 'string', 'h', 'edgecolor', 'none', ...
    'fontsize', 24, 'fontweight', 'b')

if save
    set(fig, 'PaperPositionMode','auto')     %# WYSIWYG
    print(fig, '-painters', '-depsc', 'Manuscript/figures_2ndRound_revision1/Figure8.eps')
end

%%
function figure9(save)
% Figure 9: Representative spectra and X+Y and Z through time

% EC_type = '';
EC_type = 'short';

load(['data/spectrogram_data_withMuscle_', EC_type, 'EC.mat']);
load(['data/NFTparams_data_', EC_type, 'EC.mat'])
bands = model.eeg_bands();

fexp = w_exp/(2*pi);
ffit = w_finer/(2*pi);
if strcmpi(EC_type, 'short')
    time_slices = [10, 24, 40];
else
    time_slices = [10, 40, 70];
end
ylims = [1e-3, 3];
xlims = [0.1, 30];
line_thick = [0.5, 1.2, 2.5];
colors = ['k', 'r', 'b'];

fig = figure;
pos = get(gcf, 'pos');
delete(fig)
fig = figure('pos', pos.*[1 1 1.2 1]);
loglog(ffit, spectrogram.fit_James_muscle_finer(:, time_slices(1)), colors(1), 'Linewidth', line_thick(1))
text(0.11, spectrogram.fit_James_muscle_finer(1, time_slices(1))*1.2, '1', 'fontsize', 16, 'fontweight', 'b')
hold on;
for j = 2:length(time_slices)
    loglog(ffit, spectrogram.fit_James_muscle_finer(:, time_slices(j)), colors(j), 'Linewidth', line_thick(j))
    text(0.11, spectrogram.fit_James_muscle_finer(1, time_slices(j))*1.2, num2str(j), 'fontsize', 16, 'fontweight', 'b')
    
    xlabel('$f$ {\rm (Hz)}', 'fontsize', 15, 'interpreter', 'latex')
    ylabel('$P(f)$ {\rm (arbitrary units)}', 'fontsize', 15, 'interpreter', 'latex')
    
    set(gca, 'fontsize', 15, 'XLim', xlims, 'XTick', [0.1, 0.5, 1,2,3,5,10,20,40], ...
        'YLim', ylims, 'ticklength', [0.02, 0.02], 'layer', 'top')    
end
text(0.55, 2.1, 'low', 'interpreter', 'latex', 'fontsize', 18, 'fontweight', 'b')
text(7.8, 2.1, 'alpha', 'interpreter', 'latex', 'fontsize', 18, 'fontweight', 'b')

fill([bands(1,1)+0.1 bands(1,2) bands(1,2) bands(1,1)+0.1], [ylims(1) ylims(1) ylims(2) ylims(2)], ...
     [0.95,0.95,0.95],'EdgeColor','k');
fill([bands(3,1) bands(3,2) bands(3,2) bands(3,1)], [ylims(1) ylims(1) ylims(2) ylims(2)], ...
     [0.95,0.95,0.95],'EdgeColor','k');
hold off;
set(gca,'children',flipud(get(gca,'children')))

% inset
HAx = axes('Position', [0.23, 0.23, 0.18, 0.18]);
data_x = NFTparams.xyz(1,:) + NFTparams.xyz(2,:);
data_y = NFTparams.xyz(3,:);
text_x = '$X + Y$';
text_y = '$Z$';

plot(data_x(time_slices), data_y(time_slices), 'ko', 'markerfacecolor', 'k', 'markersize', 7)
text(data_x(time_slices(1))*0.98, data_y(time_slices(1)), '1', 'fontsize', 10, 'fontweight', 'b')
text(data_x(time_slices(2))*1.012, data_y(time_slices(2)), '2', 'fontsize', 10, 'fontweight', 'b')
text(data_x(time_slices(3))*1.01, data_y(time_slices(3))*1.02, '3', 'fontsize', 10, 'fontweight', 'b')
utils.arrow([data_x(time_slices(1)) data_y(time_slices(1))], [data_x(time_slices(2))*0.995 data_y(time_slices(2))], ...
    'length', 7, 'baseangle', 45)
utils.arrow([data_x(time_slices(2)) data_y(time_slices(2))], [data_x(time_slices(3))*0.996 data_y(time_slices(3))*1.03], ...
    'length', 7, 'baseangle', 45)

xlabel(text_x, 'fontsize', 15, 'interpreter', 'latex')
ylabel(text_y, 'fontsize', 15, 'interpreter', 'latex')
if strcmpi(EC_type, 'short')
    set(gca, 'fontsize', 15, 'XLim', [min(data_x)*0.98 max(data_x)*1.02], ...
        'YLim', [min(data_y)*0.98 max(data_y)*1.02], 'ticklength', [0.02, 0.02])
else
    set(gca, 'fontsize', 15, 'XLim', [min(data_x)*0.85 max(data_x)*1.15], ...
        'YLim', [min(data_y)*0.85 max(data_y)*1.15], 'ticklength', [0.02, 0.02])
end

if save
    set(fig, 'PaperPositionMode','auto')     %# WYSIWYG
    print(fig, '-painters', '-depsc', 'Manuscript/figures_2ndRound_revision1/Figure9.eps')
end

%%
function figure10(save)
% Figure 10: Correlation of spectral properties and feedback loops

EC_type = 'short';

load(['data/spectralProperty_data_', EC_type, 'EC.mat']);
load(['data/phi_data_', EC_type, 'EC.mat'])
load(['data/NFTparams_data_', EC_type, 'EC.mat'])

data_type = 'fit_James_finer';

% attribs = {'int_delta', 'int_alpha', 'alpha_maxp', 'alpha_str', ...
%            'X', 'Y', 'X+Y', 'Z'};
% symbols = {'$P_{\rm low}$', '$P_\alpha$', '$P(f^{\rm{max}}_\alpha)$', ...
%            '$\log(R_\alpha)$', '$X$', '$Y$', '$X + Y$', '$Z$'};
       
attribs = {'int_delta', 'int_alpha', 'alpha_maxp', 'alpha_str', ...
           'Z', 'X+Y'};           
symbols = {'$P_{\rm low}$', '$P_\alpha$', '$P(f^{\rm{max}}_\alpha)$', ...
           '$\log(R_\alpha)$', '$Z$', '$X + Y$'};
       
correlations = zeros(length(attribs), length(attribs));
pvals = zeros(length(attribs), length(attribs));
for i = 1:length(attribs)
    if strcmpi(attribs{i}, 'alpha_str')
        data_x = log10(spectral_prop.(data_type).(attribs{i}));
    elseif strcmpi(attribs{i}, 'X')
        data_x = NFTparams.xyz(1, :);
    elseif strcmpi(attribs{i}, 'Y')
        data_x = NFTparams.xyz(2, :);
    elseif strcmpi(attribs{i}, 'X+Y')
        data_x = NFTparams.xyz(1, :) + NFTparams.xyz(2, :);    
    elseif strcmpi(attribs{i}, 'Z')
        data_x = NFTparams.xyz(3, :);
    else
        data_x = spectral_prop.(data_type).(attribs{i});
    end
    
    for j = i:length(attribs)
        if strcmpi(attribs{j}, 'alpha_str')
            data_y = log10(spectral_prop.(data_type).(attribs{j}));
        elseif strcmpi(attribs{j}, 'X')
            data_y = NFTparams.xyz(1, :);
        elseif strcmpi(attribs{j}, 'Y')
            data_y = NFTparams.xyz(2, :);
        elseif strcmpi(attribs{j}, 'X+Y')
            data_y = NFTparams.xyz(1, :) + NFTparams.xyz(2, :);    
        elseif strcmpi(attribs{j}, 'Z')
            data_y = NFTparams.xyz(3, :);
        else
            data_y = spectral_prop.(data_type).(attribs{j});
        end
        
        [r, p] = corrcoef(data_x, data_y);
        correlations(j, i) = r(1, 2);
        pvals(j, i) = p(1, 2);
    end
end

% repeat last column and row for pcolor
correlations = cat(1, correlations, correlations(end,:));
correlations = cat(2, correlations, correlations(:,end));

fig = figure;
h = pcolor(correlations);
set(h, 'Edgecolor', 'none')
caxis([-1, 1])
cmap = utils.redblue;
colormap(cmap)
c = colorbar('north');
set(c, 'position', get(c, 'position').*[0 0 0 0] + [0.53 0.8 0.35 0.05])
% title(c, 'correlation')
set(gca, 'ydir', 'reverse', 'Fontsize', 15, 'Xtick', [], 'Ytick', [], ...
    'ticklength', [0.02, 0.02], 'box', 'off');

for i = 1:length(attribs)
    for j = i:length(attribs)
        if i~=j
            text(i+0.5, j+0.3, num2str(correlations(j,i), '%0.2f'), ...
                'fontsize', 18, 'horizontalalignment', 'center', 'color', 'w')
            [number, exponent] = utils.scientific_notation(pvals(j,i));
            text(i+0.5, j+0.7, ['(',num2str(number), num2str('\times'), sprintf('10^{%d}', exponent),')'], ...
                'fontsize', 12, 'horizontalalignment', 'center', 'color', 'w')
        else
            text(i+0.5, j+0.5, num2str(correlations(j,i)), ...
                'fontsize', 18, 'horizontalalignment', 'center', 'color', 'w')
        end
    end
end

for i=1:length(attribs)
    text(0.9, i+0.5, symbols{i}, 'interpreter', 'latex', 'fontsize', 15, ...
        'HorizontalAlignment', 'right')
    text(i+0.5, length(attribs) + 1.2, symbols{i}, 'interpreter', 'latex', 'fontsize', 15, ...
        'HorizontalAlignment', 'right', 'Rotation',30)
end

for j=2:length(attribs)+1
    line([1 j], [j-1 j-1], 'color', 'k', 'linewidth', 0.1)
    line([j j], [j-1 length(attribs)+1], 'color', 'k', 'linewidth', 0.1)
end
line([1 1], [1 length(attribs)+1], 'color', 'k', 'linewidth', 0.1)
line([1 length(attribs)+1], [length(attribs)+1 length(attribs)+1], 'color', 'k', 'linewidth', 0.1)

if save
    set(fig, 'PaperPositionMode','auto')     %# WYSIWYG
    print(fig, '-painters', '-depsc', 'Manuscript/figures_2ndRound_revision1/Figure10.eps')
%     utils.fix_pcolor_eps('Manuscript/figures_2ndRound_revision1/Figure10.eps')
end

