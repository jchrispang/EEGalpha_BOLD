function BiCo = bicoherence(f,fftx,drawbar,surrogate)
    % BiCo = bicoherence(f,fftx,drawbar,surrogate)
    % Make bicoherence plot given array of FFT values and the frequency array
    %
    % f is an array of frequencies
    % fftx is an array with FTs of the time series
    % drawbar = 1 to draw a colour bar
    % surrogate = 1 to perform surrogate data testing. In this case, the plot
    % is of bicoherence significance rather than actual bicoherence
    if nargin < 4 || isempty(surrogate)
        surrogate = 0;
    end

    if nargin < 3 || isempty(drawbar)
        drawbar = 1; % Draw colorbar by default
    end

    % Fr provides a list of frequency values to make the bicoherence plot for
    % Since we have 200Hz data and the Nyquist frequency is 100Hz, it only really
    % makes sense to go up to 50 Hz
    Fr = f(f<=50); % Calculate bicoherence for f <= 50

    if surrogate
        nsurr = 100;
        bico_base = get_bicoherence(Fr,fftx,f);
        bico_surr = zeros(size(bico_base,1),size(bico_base,2),nsurr);
        for j = 1:nsurr
            bico_surr(:,:,j) = get_bicoherence(Fr,randomize(fftx),f);
        end
        %keyboard
        bico_surr_mean = mean(bico_surr,3);
        bico_surr_std = std(bico_surr,1,3);
        BiCo = (bico_base-bico_surr_mean)./bico_surr_std;
        BiCo = abs(BiCo); % Have verified FOR SPINDLE DATA the sign of the bicoherence. But for significance, either way is significant
    else
        BiCo = get_bicoherence(Fr,fftx,f);
    end

    figure;
    imagesc([Fr(1) Fr(end)],[Fr(1) Fr(end)],BiCo)
    axis xy
    set(gca,'FontSize',28);
    xlabel('f_1 (Hz)','FontSize',28)
    ylabel('f_2 (Hz)','FontSize',28)
    axis equal
    axis tight
    if drawbar || surrogate % Surrogate data has variable colour axis so it must be displayed
        cb = colorbar;
        set(cb,'FontSize',28);
        if surrogate
            ylabel(cb,'Significance (SD)','FontSize',28);
            set(gca,'CLim',[1 15])
        else
            set(gca,'CLim',[0 1])
            ylabel(cb,'Bicoherence','FontSize',28);
        end
    end
    
    % BW plot for printing
    cm = colormap(gray);
    colormap(flipud(cm));
end

function r = randomize(fft_data)
    r = fft_data./exp(i*angle(fft_data)).*exp(i*2*pi*rand(size(fft_data)));
end

function BiCo = get_bicoherence(Fr,fft_data,f)
    % Take in the fft data, return the bicoherence matrix
    BiCo = zeros(length(Fr));
    for f1 = 1:length(Fr),
       for f2 = 1:length(Fr), 
                BiCo(f1,f2) = bico_point(fft_data,f,Fr(f1),Fr(f2));
                BiCo(f2,f1)=BiCo(f1,f2);
       end
    end
end

function hayashi = bico_point(fft_data,f,f1,f2)
    % See http://www-fusion.ciemat.es/fusionwiki/index.php/Bicoherence for notation
    [~,f1_ind] = min(abs(f-f1));
    [~,f2_ind] = min(abs(f-f2));
    [~,f3_ind] = min(abs(f-f1-f2));
    ft{1} = fft_data(f1_ind,:);
    ft{2} = fft_data(f2_ind,:);
    ft{3} = fft_data(f3_ind,:);
    
    %fusionwiki = abs(mean(fft_data(f1_ind,:).*fft_data(f2_ind,:).*conj(fft_data(f3_ind,:)))).^2;
    %fusionwiki = fusionwiki./mean(abs(fft_data(f1_ind,:).*fft_data(f2_ind,:)).^2)./mean(abs(fft_data(f3_ind,:)).^2);
    %fusionwiki = sqrt(fusionwiki);
    
    %astro = abs(sum(ft{1}.* ft{2} .* conj(ft{3}))).^2;
    %astro = astro./(sum(abs(ft{1}.*ft{2}).^2).*sum(abs(ft{3}).^2));
    %astro = sqrt(astro);
    
    %sigl = abs(sum(ft{1} .* ft{2} .* conj(ft{3})));
    %sigl = sigl./sqrt(sum(abs(ft{1}).^2.*abs(ft{2}).^2.*abs(ft{3}).^2));
    
    %hagihira = abs(sum(ft{1} .* ft{2} .* conj(ft{3})));
    %hagihira = hagihira./sum(sqrt(abs(ft{1}).^2.*abs(ft{2}).^2.*abs(ft{3}).^2));
    
    % From Hayashi 2007
    TP = ft{1} .* ft{2} .* conj(ft{3});
    hayashi = abs(sum(TP));
    hayashi = hayashi./sum(abs(TP));
    %hayashi = hayashi.^2;

    %barnett = mean(ft{1}.* ft{2} .* conj(ft{3})).^2;
    %barnett = barnett./mean(ft{1}.^2 .* ft{2}.^2.* ft{3}.^2);
    
    %fprintf('Fusionwiki: %f Astro: %f Sigl: %f Hagihira: %f Hayashi: %f Barnett: %f\n',fusionwiki,astro,sigl,hagihira,hayashi,barnett);
end
