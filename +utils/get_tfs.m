function [ts,fv,spectra] = get_tfs(x,edata,interval,decayrate)	
		% Use spectral averaging to calculate time series of spectra
		% manually set quiet = 1 for no output
		quiet = 1;
		
		if nargin < 4
			decaymode = 0;
			decayrate = 4; % Using decayrate=4 causes the last section to be about 2% weighted
		else
			decaymode = 1; % Use decay if a decay rate was specified
		end
		
		if nargin < 3
			interval = 30; % Number of seconds to average over
		end

		tic;
		
		plot_step= 1; % Time separation between successive plots (s). Default 1s

		stop = floor(x(end)); % Time to stop plots

		high = 45; % Max frequency (Hz) - Normally leave at 50, 
				   % this has been reduced to remove artifact spike
		low = 0.5; % Min frequency (Hz) % THIS HAS BEEN CHANGED
	
		% Initialise the spectrum averaging matrix
		window = 15; % Size of each window (seconds)
		stop = stop-window; % Prevent overrunning the end of the data
		overlap = window-plot_step; % Overlap between windows (seconds)
		
		n = (interval-window)/plot_step+1; % Number of slices to combine
		if quiet==0
			fprintf(1,'Processing\nSpectrum is an average of %i spectra each %i seconds long\nTotal %i s per spectrum with %.2f s separation\n',n,window,interval,plot_step);
		end
		% Offset of 1 corrects for overlap error. Note that modifying plot_step
		% really is changing the number of 
		
		rate = 1/x(2); % The sampling frequency. window*rate is the number of points needed
		if rate-round(rate) > 1e-5 || n ~= round(n)
			throw(MException('FloatError:NotEqual','ERROR with sampling rate calculation'));
		else
			rate = round(rate);
		end
		npts = length(pwelch(edata(1:rate*window),[],[],[],rate));
		data = zeros(npts,n); % Holding matrix for the spectra. This will be used later!
	
		% Prepare the time decay
		decay = exp(-[0:n-1]/(interval/decayrate)); % Using interval/4 causes the last section to be about 2% weighted
		decay = decay(end:-1:1); % Reverse the array because large-n (right side of data matrix) corresponds to the latest time
		totalweights = sum(decay); % Decay needs to be repeated to operate on each of the data points
		decay = repmat(decay,npts,1); % Now it is the same dimension as data

		
		% Populate the keyring
		for j = 1:n
			wstart = ((window-overlap) * (j-1)) * rate + 1;
			if abs(round(wstart)-wstart)>1e-7;error('Integer error (stage 1)!'); else; wstart = round(wstart); end; % Correct Matlab precision error when using non-integer plot_step
			wstop = wstart + (window*rate);
			[spec,freq] = pwelch(edata(wstart:wstop),[],[],[],rate);
			if ~all(isfinite(spec)); error('NaN error in pwelch!');end;
			data(:,j) = spec;
		end
		
		if decaymode
			s = sum(data.*decay,2)./totalweights; % Weighted average spectrum
		else
			s = sum(data,2)./n; % Average spectrum
		end
		
		s = s(1:find(freq<=high, 1, 'last' )); % truncate to <= high Hz
		freq = freq(1:find(freq<=high, 1, 'last' ));
		s = s(find(freq>=low, 1, 'first' ):end); % truncate to >= low Hz
		freq = freq(find(freq>=low, 1, 'first' ):end);
	
		% Create the fv matrices
		fv = freq;
		ts = interval:plot_step:stop;
		spectra = zeros(length(freq),length(ts));
		spectra(:,1) = s;

		% Do the loop
		count = -1;
		st = interval;

		while st < stop
			count = count+1; % Increment the counter
			st = st+plot_step;
			% What is the next slice that is required? It is the next increment which starts at
			
			wstart = st*rate;
			if abs(round(wstart)-wstart)>1e-7;error('Integer error (stage 2)!'); else; wstart = round(wstart); end; % Correct Matlab precision error when using non-integer plot_step
			
			wstop  = wstart + window*rate;
			
			[spec,freq] = pwelch(edata(wstart:wstop),[],[],[],rate);
			if ~all(isfinite(spec)); error('NaN error in pwelch!');end;
						
			data = circshift(data,[0 -1]); % Move columns left
			data(:,end) = spec;
			
			% Update the spectral data
			if decaymode
				s = sum(data.*decay,2)./totalweights; % Weighted average spectrum
			else
				s = sum(data,2)./n; % Average spectrum
			end

			s = s(1:find(freq<=high, 1, 'last' )); % truncate to <= high Hz
			s = s(find(freq>=low, 1, 'first' ):end); % truncate to >= low Hz
			
			spectra(:,count+2) = s;
			
		end
		if quiet == 0
			fprintf(1,'Finished in %.0f seconds\n',toc);
		end
