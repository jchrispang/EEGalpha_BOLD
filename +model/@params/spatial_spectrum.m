function [f,P,grid_kx,grid_ky,Pkf,output_x,output_y,Prf] = spatial_spectrum(p,n_modes,n_truncate,output_x,output_y,f,emg_component)
	% output_x and output_y should be meshgrid output
	% The return arguments are matched to nf_spatial_spectrum
	% Note that the K-spectrum is returned with appropriate fftshifts
	% Note that if truncation is specified, the output x and y arrays
	% still correspond to the original n_modes, unless manually specified!

	% Note that EMG gets added to Prw only!!

	if nargin < 6 || isempty(f)
		f = linspace(0,45,200);
	end
	
	if nargin < 2 || isempty(n_modes)
		n_modes = 6;
	end

	if nargin < 3 || isempty(n_truncate)
		n_truncate = n_modes;
	end

	% Initial setup
	[~,grid_kx,grid_ky,grid_x,grid_y,~,dk,dx]=get_frequencies(zeros(n_modes,n_modes,1),1,p.Lx,p.Ly);

	w(1,1,:) = f*2*pi;
	dw = w(2)-w(1);
	df = dw/2/pi;

	if nargin < 7 || isempty(emg_component)
		emg_f = 40;
		emg_component = (w/(2*pi*emg_f)).^2./(1+(w/(2*pi*emg_f)).^2).^2;
	end

	if nargin < 4 || isempty(output_x)
		output_x = grid_x;
		output_y = grid_y;
	end

	grid_kx = ifftshift(grid_kx);
	grid_ky = ifftshift(grid_ky); % Set up these arrays so that the forward transform doesn't need shifting

	spec_gamma_prefactor = bsxfun(@times,(1-1i*w/p.gammae).^2,ones(size(grid_x)));

	if ~isempty(p.spatial_xyz) || (isempty(p.gab) && isempty(p.gabcd))
		[X,Y,Z,Zp,L,M] = p.get_xyz(w,grid_x,grid_y);
		Afull = spec_gamma_prefactor - X - Y.*(1+Zp)./(1+Zp.*L.*L).*M.*M;
		Bfull = 1./(1+Zp.*L.*L); % Note that 1-Jei is just a number, which is the same as rolling it into phin
	else
		[Jee,Jei,Jese,Jesre,Jsrs,Jesn] = p.get_jabcd(w,grid_x,grid_y);
		Afull = spec_gamma_prefactor - Jee./(1-Jei) - (Jese + Jesre)./((1-Jei).*(1-Jsrs));
		Bfull = Jesn./((1-Jei).*(1-Jsrs));
	end
	
	Afull = fft(Afull,[],1);
	Afull = fft(Afull,[],2);
	Afull = Afull/numel(grid_kx); 

	Bfull = fft(Bfull,[],1);
	Bfull = fft(Bfull,[],2);
	Bfull = Bfull/numel(grid_kx); % verified correct normalization

	if n_truncate<n_modes
		% If truncation is required
		idx = [1 2:floor((n_truncate+1)/2) n_modes-ceil((n_truncate-1)/2)+1:n_modes];
		grid_kx = grid_kx(idx,idx);
		grid_ky = grid_ky(idx,idx);
		Afull = Afull(idx,idx,:);
		Bfull = Bfull(idx,idx,:);
	end

	grid_K2 = grid_kx.^2+grid_ky.^2;
	k0 = 10; % Volume conduction parameter
	k2_volconduct = exp(-grid_K2/k0^2);
	%k2_volconduct = ones(size(k2_volconduct)); % USE THIS TO DISABLE VOLUME CONDUCTION

	% Get the convolution matrix - see conv_compare.m in rev 1435 
	Qm = fft(eye(size(grid_K2,1)));
	Q = kron(Qm,Qm);
	Lm = fft2(reshape(1:numel(grid_K2),size(grid_K2))); % Note - need to undo the  applied to ap initially
	T = round(Q'*diag(Lm(:))*Q/numel(Lm)); 

	phin = p.phin*ones(numel(grid_K2),1);
	Pkw = zeros([size(grid_K2) length(f)]);
	Prw = zeros([size(output_x) length(f)]);
	re2 = p.re.^2;

	% global magic_P magic_P2
	% magic_P = zeros(length(f),1);
	% magic_P2 = zeros(length(f),1);

	for j = 1:size(Afull,3)
		% Go through and calculate the spectrum, one frequency at a time
		A_temp = Afull(:,:,j);
		B_temp = Bfull(:,:,j);

		A = A_temp(T);
		A = A+diag(grid_K2(:)*re2);
		B = B_temp(T);

		M = A\B; % inv(A)*B == A\B

		% Next, we need to get back the same entries from conv. This gives the same output indexes as conv2(x,y,'same')
		% Verified that this reshape gives back the equivalent result from grid_K2(:) ordering

		phie_kw = reshape(M*phin,size(grid_K2));
		phie_kw = phie_kw.*sqrt(k2_volconduct);
		Pkw(:,:,j) = abs(phie_kw).^2;

		if nargout > 5 || nargout == 0
			% mm_filter will be different on both the first and second iterations
			mm = M*M';

			% if w(j) > 8
			% 	keyboard
			% end

			mm = bsxfun(@times,mm,k2_volconduct(:)');
			
			if j == 1 || j == 2
				pt = zeros(size(output_x));
				mm_filter = abs(mm) > 0;
				kgx = bsxfun(@minus,grid_kx(:),grid_kx(:)');
				kgy = bsxfun(@minus,grid_ky(:),grid_ky(:)');
				kgx1 = kgx(mm_filter);
				kgy1 = kgy(mm_filter);
				exponent_cache = zeros(length(kgx1),numel(output_x));
				for j_1 = 1:numel(output_x)
					exponent_cache(:,j_1) = exp(1i*kgx1*output_x(j_1) + 1i*kgy1*output_y(j_1));
				end
			end

			mm_terms = mm(mm_filter);
			
			for j_1 = 1:numel(output_x) % For each position
				tmp = exponent_cache(:,j_1) .* mm_terms;			
				pt(j_1) = abs(sum(tmp(:)));
			end

			Prw(:,:,j) = pt;
			%keyboard
			% magic_P(j) = abs(p.phin).^2*trace(M*M');
			% magic_P2(j) = (abs(M)*phin)'*(abs(M)*phin);

			% keyboard

		end

	end
	
	phin_multiply = sum(abs(phin).^2)*dk*dk; % Get phin(w) from phin(k,w)
	Prw = Prw.*phin_multiply/numel(grid_kx)*dk*dk/2/pi/2/pi;
	
	if utils.isfunction(p.spatial_emg)
		emg_component = bsxfun(@times,emg_component,p.spatial_emg(output_x,output_y));
	else
		emg_component = p.emg_a*emg_component;
	end

	Prw = bsxfun(@plus,Prw,emg_component);
	Prf = Prw.*2*pi;


	% Normal power spectrum integrated over k
	Pkf = Pkw.*2*pi;
	P = sum(sum(Pkf*dk^2,1),2);
	P = P(:).'; % Convert to P(f)
	
	grid_kx = fftshift(grid_kx);
	grid_ky = fftshift(grid_ky);
	Pkf = fftshift(Pkf,1);
	Pkf = fftshift(Pkf,2);

	if nargout == 0
		utils.prw_imager(output_x,output_y,f,Prf);
	end

function [f,Kx,Ky,x,y,df,dk,dx] = get_frequencies(data,fs,Lx,Ly)
	% Take in the size of the data matrix, the sampling rate, and spatial distances
	% Return corresponding grids of frequencies and positions
	Kx = [];
	Ky = [];
	x = [];
	y = [];

	if isvector(data)
		f = (0:fs/length(data):fs/2).'; % This is the single sided frequency
		return
    end

    f = (0:fs/size(data,3):fs/2).'; % This is the single sided frequency
    df = fs/size(data,3);
    
    Kx = 2*pi*(0:1/Lx:size(data,1)/Lx/2); % This is the single sided frequency
    Ky = 2*pi*(0:1/Ly:size(data,2)/Ly/2); % This is the single sided frequency
    
    if mod(size(data,1),2) % If there is NO nyquist frequency component
    	Kx = [-Kx(end:-1:2) Kx];
    else
    	Kx = [-Kx(end:-1:2) Kx(1:end-1)];
    end

    if mod(size(data,2),2) % If there is NO nyquist frequency component
        Ky = [-Ky(end:-1:2) Ky];
    else
        Ky = [-Ky(end:-1:2) Ky(1:end-1)];
    end
    [Kx,Ky] = meshgrid(Kx,Ky);
    dk = 2*pi/Lx;

    dx = Lx/size(data,1); % Metres per pixel
    x = dx*(0:size(data,1)-1);
    dy = Ly/size(data,1); 
    y = dy*(0:size(data,1)-1);
    [x,y] = meshgrid(x,y);
