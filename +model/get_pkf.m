function output = get_pkf(f,P)
	% Find peaks in function P of f AND sanitize to pkf format
	%
	% Resulting pkf array is
	% 	[maxf, maxp, minf, minp, pkstr] for each of the 5 bands (from eeg_bands())
	%   maxf = frequency of maximum power, maxp = value of maximum power
	P = P(:);
	f = f(:);
    
	roots = find(diff(sign(diff(P))));
	roots = roots.*isfinite(roots);	
	deriv_2 = diff(diff(P));
	maxima = roots(find(deriv_2(roots) <= 0))+1; % Nature of the roots- these are peaks
	minima = roots(find(deriv_2(roots) > 0))+1; % These are minima

    % figure
    % plot(log(f),log(P));
    % hold on
    % scatter(log(f(maxima)),log(P(maxima)),30,'ro');
    % scatter(log(f(minima)),log(P(minima)),30,'go');

    if length(maxima) == 1 && isempty(minima) % The first and only maximum has no minimum
        minima = 1;
    elseif ~(isempty(minima) || isempty(maxima)) % If there are some maxima or minima
        if maxima(1) < minima(1) % If the first maximum has no minimum preceding it
            minima = [1; minima];
        end

        if minima(end) > maxima(end) % If the last minimum is not followed by a maximum
            minima = minima(1:end-1);
        end
    else % If there is a minimum with no local maximum, there are no peaks to be detected
        maxima = [];
        minima = [];
    end

    [bands,bandstr] = model.eeg_bands();

    pk_f = [f(minima),f(maxima)];
    pk_P = [P(minima),P(maxima)];
    indices = 1:length(minima);

    selected = zeros(length(bandstr),2); % first column is minimum index, second is maximum index

    for j = 1:length(bandstr)
        out = nan(1,5);
    	% PROCEDURE
    	% Pick the local maximum with the highest power
    	% And the local minimum with the lowest power that is between the candidate maximum peak and 
    	if ~isempty(pk_f)
            in_range = (pk_f(:,2) >= bands(j,1) & pk_f(:,2) <= bands(j,2));

        	if ~any(in_range)
        		selected(j,:) = [NaN NaN];
        	else
        		tmp = pk_P(:,2);
        		tmp(~in_range) = NaN;
        		[~,selected(j,2)] = max(tmp); 

        		if j >= 2 && isfinite(selected(j-1,2))
        			in_range = (indices(:)>selected(j-1,2)).*(indices(:)<=selected(j,2));
        		else
        			in_range = (indices(:)<=selected(j,2));
        		end
                        
        		tmp = pk_P(:,1);
        		tmp(~in_range) = NaN;
        		[~,selected(j,1)] = min(tmp); % The index of the local minimum with the lowest power
        	end

        	% str, maxf, maxp, minf, minp
        	if all(isfinite(selected(j,:)))
        		out = [pk_P(selected(j,2),2)/pk_P(selected(j,1),1) pk_f(selected(j,2),2) pk_P(selected(j,2),2) pk_f(selected(j,1),1) pk_P(selected(j,1),1)];
        	end
        end
    	output.(sprintf('%s_str',bandstr{j})) = out(1);
    	output.(sprintf('%s_maxf',bandstr{j})) = out(2);
    	output.(sprintf('%s_maxp',bandstr{j})) = out(3);
    	output.(sprintf('%s_minf',bandstr{j})) = out(4);
    	output.(sprintf('%s_minp',bandstr{j})) = out(5);
    end
