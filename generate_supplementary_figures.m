function SupplementaryMaterial_Figures_2ndRound_revision1(subject_num, save)
% Figures for Supplementary Material

%% Load data
output = load(['data/Tracking_state=EC_subject=', num2str(subject_num), '.mat']);
experiment = load(['data/TrackingExperimentalData_state=EC_subject=', num2str(subject_num), '.mat']);

for i=1:size(output.spectrogram,2)
    Ptotal(i) = trapz(output.w_finer, output.spectrogram(:, i));
end

%% Preliminaries

if subject_num==1
    ind_start = 945;
elseif subject_num==2
    ind_start = 60;
elseif subject_num==3
    ind_start = 852;
elseif subject_num==7
    ind_start = 12;    
else
    ind_start = 1;
end
track_length = 50;
ind_end = ind_start + (track_length - 1);

new_spectrogram = output.spectrogram(:, ind_start:ind_end)./Ptotal(ind_start:ind_end);
f = output.w_finer/(2*pi);
time_slices = [10, 25, 40];
colors = ['k', 'r', 'b'];

attribs_spectral = {'int_delta', 'int_alpha', 'alpha_maxp', 'alpha_str'};
attribs_spectral_symbols = {'$P_{\rm low}$', '$P_\alpha$',...
                   '$P(f^{\rm{max}}_\alpha)$', '$\log(R_\alpha)$'};
attribs_loop = {'X', 'Y', 'X+Y', 'Z'};
attribs_loop_symbols = {'$X$', '$Y$', '$X + Y$', '$Z$'};
attribs_corr = {'int_delta', 'int_alpha', 'alpha_maxp', 'alpha_str', ...
           'Z', 'X+Y'};           
symbols_corr = {'$P_{\rm low}$', '$P_\alpha$', '$P(f^{\rm{max}}_\alpha)$', ...
           '$\log(R_\alpha)$', '$Z$', '$X + Y$'};


%% Figure 1: Spectogram
fig1 = figure;
imagesc(1:size(new_spectrogram, 2), f, log10(new_spectrogram));
colormap(gca, utils.parula)
colorbar('ytick', -5:1, 'yticklabel', {'10^{-5}', '10^{-4}', '10^{-3}', '10^{-2}', '10^{-1}', '10^{0}',  '10^{1}'})
% caxis(log10(clims))
xlabel('$t$ {\rm (s)}', 'fontsize', 15, 'interpreter', 'latex')
ylabel('$f$ {\rm (Hz)}', 'fontsize', 15, 'interpreter', 'latex')
set(gca, 'ydir', 'normal', 'Fontsize', 15, 'YLim',[f(1) 45], ...
    'Xtick', [1,10:10:100], 'ticklength', [0.02, 0.02])
annotation('textbox', [0.02, 0.9, 0.1, 0.1], 'string', 'a', 'edgecolor', 'none', ...
        'fontsize', 24, 'fontweight', 'b')
    
%% Figure 2: Representative spectra
fig2 = figure;
loglog(f, new_spectrogram(:, time_slices(1)), colors(1), 'Linewidth', 2)
hold on;
for j = 2:length(time_slices)
    loglog(f, new_spectrogram(:, time_slices(j)), colors(j), 'Linewidth', 2)
    
    xlabel('$f$ {\rm (Hz)}', 'fontsize', 15, 'interpreter', 'latex')
    ylabel('$P(f)$ {\rm (arbitrary units)}', 'fontsize', 15, 'interpreter', 'latex')
    
    set(gca, 'fontsize', 15, 'XLim', [0.1 45], 'XTick', [0.1, 0.5, 1,2,3,5,10,20,40], ...
         'YLim', [0.5e-5 1], 'ticklength', [0.02, 0.02])   
end
leg = legend(['$t =$ ', num2str(time_slices(1)), ' s'], ...
             ['$t =$ ', num2str(time_slices(2)), ' s'], ...
             ['$t =$ ', num2str(time_slices(3)), ' s']);
set(leg, 'FontSize', 20, 'Location', 'SouthWest', 'interpreter', 'latex');
hold off;
annotation('textbox', [0.02, 0.9, 0.1, 0.1], 'string', 'b', 'edgecolor', 'none', ...
        'fontsize', 24, 'fontweight', 'b')
 
%% Figure 3: Time series of spectral and NFTparams prop
% fig3 = figure('Position', [200 200, 600, 800]);
fig3 = figure;
for j = 1:4
    % spectral properties
    s = subplot(4, 2, 2*j-1, 'ColorOrder',[0,0,0; 1,0,0]);
    if strcmpi(attribs_spectral{j}, 'alpha_str')
        data_spectral = log10(output.spectral_prop.(attribs_spectral{j})./Ptotal);
        data_exp = log10(experiment.expt.spectral_prop.(attribs_spectral{j}));
    else 
        data_spectral = output.spectral_prop.(attribs_spectral{j})./Ptotal;
        data_exp = experiment.expt.spectral_prop.(attribs_spectral{j});
    end
    data_spectral = data_spectral(ind_start:ind_end);
    data_exp = data_exp(ind_start:ind_end);
    
    yyaxis(s, 'left')
    plot(data_spectral, 'k-', 'Linewidth', 2)
    xlabel('$t$ {\rm (s)}', 'fontsize', 15, 'interpreter', 'latex')
    ylabel(attribs_spectral_symbols{j}, 'fontsize', 15, 'interpreter', 'latex')
    
    if j~=4
        set(gca, 'fontsize', 15, 'xlim', [1,length(data_spectral)], ...
        'ylim', [min(data_spectral)*0.85, max(data_spectral)*1.15], ...
        'Xtick', [1,10:10:100],  'ticklength', [0.02, 0.02], ...
        'ytick', round(linspace(min(data_spectral), max(data_spectral), 3),2))
        set(gca, 'XtickLabel', {})
        xlabel('')
    else
        set(gca, 'fontsize', 15, 'xlim', [1,length(data_spectral)], ...
        'ylim', [min(data_spectral)*0.95, max(data_spectral)*1.05], ...
        'Xtick', [1,10:10:100],  'ticklength', [0.02, 0.02], ...
        'ytick', round(linspace(min(data_spectral), max(data_spectral), 3),2))
    end
    
    [r, ~] = corrcoef(data_exp, data_spectral);
    yLimits = get(gca,'YLim');
    if subject_num == 7
        if j==1
            text(3, (yLimits(2) - yLimits(1))*0.85 + yLimits(1), sprintf('r = %0.2f', r(1,2)), ...
                    'fontsize', 12, 'horizontalalignment', 'left', 'color', 'k', 'fontweight', 'b')
        else
            text(length(data_spectral)-15, (yLimits(2) - yLimits(1))*0.85 + yLimits(1)*0.95, sprintf('r = %0.2f', r(1,2)), ...
                    'fontsize', 12, 'horizontalalignment', 'left', 'color', 'k', 'fontweight', 'b')
        end
    else
        if j~=1
            text(3, (yLimits(2) - yLimits(1))*0.85 + yLimits(1), sprintf('r = %0.2f', r(1,2)), ...
                    'fontsize', 12, 'horizontalalignment', 'left', 'color', 'k', 'fontweight', 'b')
        else
            text(length(data_spectral)-15, (yLimits(2) - yLimits(1))*0.85 + yLimits(1)*0.95, sprintf('r = %0.2f', r(1,2)), ...
                    'fontsize', 12, 'horizontalalignment', 'left', 'color', 'k', 'fontweight', 'b')
        end
    end
    
    yyaxis(s, 'right')
    plot(data_exp, 'r--', 'Linewidth', 1)
    
    set(gca, 'fontsize', 15, 'xlim', [1,length(data_spectral)], ...
        'ylim', [min(data_exp)*0.95, max(data_exp)*1.05], ...
        'Xtick', [1,10:10:100],  'ticklength', [0.02, 0.02], ...
        'ytick', round(linspace(min(data_exp)*0.98, max(data_exp)*1.02, 3),2))
    
    % feedback loops
    subplot(4, 2, 2*j);
    if strcmpi(attribs_loop{j}, 'X+Y')
        data_NFTparams = output.NFTparams_prop.X + output.NFTparams_prop.Y;
    else
        data_NFTparams = output.NFTparams_prop.(attribs_loop{j});
    end
    data_NFTparams = data_NFTparams(ind_start:ind_end);
    
    plot(data_NFTparams, 'k-', 'Linewidth', 2)

    xlabel('$t$ {\rm (s)}', 'fontsize', 15, 'interpreter', 'latex')
    ylabel(attribs_loop_symbols{j}, 'fontsize', 15, 'interpreter', 'latex')
    
    set(gca, 'fontsize', 15, 'xlim', [1,length(data_NFTparams)], 'Xtick', [1,10:10:100], ...
             'Ylim', [min(data_NFTparams)*0.95, max(data_NFTparams)*1.05],  'ticklength', [0.02, 0.02], ...
             'ytick', round(linspace(min(data_NFTparams), max(data_NFTparams), 3),2))

    if j~=4
        set(gca, 'XtickLabel', {})
        xlabel('')
    end
end
annotation('textbox', [0.02, 0.9, 0.1, 0.1], 'string', 'c', 'edgecolor', 'none', ...
        'fontsize', 24, 'fontweight', 'b')
    
%% Figure 4: Correlations
correlations = zeros(length(attribs_corr), length(attribs_corr));
pvals = zeros(length(attribs_corr), length(attribs_corr));
for i = 1:length(attribs_corr)
    if strcmpi(attribs_corr{i}, 'alpha_str')
        data_x = log10(output.spectral_prop.(attribs_corr{i})./Ptotal);
    elseif strcmpi(attribs_corr{i}, 'X+Y')
        data_x = output.NFTparams_prop.X + output.NFTparams_prop.Y;
    elseif strcmpi(attribs_corr{i}, 'Z')
        data_x = output.NFTparams_prop.Z;
    else 
        data_x = output.spectral_prop.(attribs_corr{i})./Ptotal;
    end
    data_x = data_x(ind_start:ind_end);

    for j = i:length(attribs_corr)
        if strcmpi(attribs_corr{j}, 'alpha_str')
            data_y = log10(output.spectral_prop.(attribs_corr{j})./Ptotal);
        elseif strcmpi(attribs_corr{j}, 'X+Y')
            data_y = output.NFTparams_prop.X + output.NFTparams_prop.Y;
        elseif strcmpi(attribs_corr{j}, 'Z')
            data_y = output.NFTparams_prop.Z;
        else 
            data_y = output.spectral_prop.(attribs_corr{j})./Ptotal;
        end
        data_y = data_y(ind_start:ind_end);
        
        [r, p] = corrcoef(data_x, data_y);
        correlations(j, i) = r(1, 2);
        pvals(j, i) = p(1, 2);
    end
end

% repeat last column and row for pcolor
correlations = cat(1, correlations, correlations(end,:));
correlations = cat(2, correlations, correlations(:,end));

fig4 = figure;
h = pcolor(correlations);
set(h, 'Edgecolor', 'none')
caxis([-1, 1])
cmap = utils.redblue;
colormap(gca, cmap)
c = colorbar('north');
set(c, 'position', get(c, 'position').*[0 0 0 0] + [0.53 0.8 0.35 0.05])
% title(c, 'correlation')
set(gca, 'ydir', 'reverse', 'Fontsize', 15, 'Xtick', [], 'Ytick', [], ...
    'ticklength', [0.02, 0.02], 'box', 'off');

for i = 1:length(attribs_corr)
    for j = i:length(attribs_corr)
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

for i=1:length(attribs_corr)
    text(0.9, i+0.5, symbols_corr{i}, 'interpreter', 'latex', 'fontsize', 15, ...
        'HorizontalAlignment', 'right')
    text(i+0.5, length(attribs_corr) + 1.2, symbols_corr{i}, 'interpreter', 'latex', 'fontsize', 15, ...
        'HorizontalAlignment', 'right', 'Rotation',30)
end

for j=2:length(attribs_corr)+1
    line([1 j], [j-1 j-1], 'color', 'k', 'linewidth', 0.1)
    line([j j], [j-1 length(attribs_corr)+1], 'color', 'k', 'linewidth', 0.1)
end
line([1 1], [1 length(attribs_corr)+1], 'color', 'k', 'linewidth', 0.1)
line([1 length(attribs_corr)+1], [length(attribs_corr)+1 length(attribs_corr)+1], 'color', 'k', 'linewidth', 0.1)

annotation('textbox', [0.02, 0.9, 0.1, 0.1], 'string', 'd', 'edgecolor', 'none', ...
        'fontsize', 24, 'fontweight', 'b')

%% Saving the figures
if save
    set(fig1, 'PaperPositionMode','auto')     %# WYSIWYG
    set(fig2, 'PaperPositionMode','auto')     %# WYSIWYG
    set(fig3, 'PaperPositionMode','auto')     %# WYSIWYG
    set(fig4, 'PaperPositionMode','auto')     %# WYSIWYG
    print(fig1, '-painters', '-depsc', ['Manuscript/figures_2ndRound_revision1/SupplementaryFig_subject=', ...
        num2str(subject_num),'_a.eps'])
    print(fig2, '-painters', '-depsc', ['Manuscript/figures_2ndRound_revision1/SupplementaryFig_subject=', ...
        num2str(subject_num),'_b.eps'])
    print(fig3, '-painters', '-depsc', ['Manuscript/figures_2ndRound_revision1/SupplementaryFig_subject=', ...
        num2str(subject_num),'_c.eps'])
    print(fig4, '-painters', '-depsc', ['Manuscript/figures_2ndRound_revision1/SupplementaryFig_subject=', ...
        num2str(subject_num),'_d.eps'])
    
    close all
end