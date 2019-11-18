% This snippet compares subsets of the actelion data
% Examine for example the first night vs. all other nights, for each subject

% SVN rev 187: Compare actelion first night to actelion other nights


study_id = [3 29 34 48 50 63 64]; 

% Number of entries in subjects{j} must be the same
subjects{1} = {[32]}; % Actelion indices (control nights)
subjects{2} = {[1:9],[1:9],[1:9],[1:9],[1:9],[1:9],[1:9]}; % Control indices (all subjects)
subjects{3} = {[33:36]}; % Actelion other nights


% Number of entries in states{j} must be the same
states{1} = {'actelion_S3'};
states{2} = {'control_S3'};
states{3} = {'actelion_S3'};

% Number of entries here matches length(states) and length(subjects
descriptor{1} = 'Actelion controls';
descriptor{2} = 'Control';
descriptor{3} = 'Actelion treatment';

for q = 1:length(states{1}) % For each state

    for j = 1:length(subjects{1}) % For each of the Actelion subjects
        data = {};
        labels = {};

        for k = 1:length(subjects) % For each of the studies that is being examined
            % At this level, k decides which entry of subjects,states and descriptor is used
            % j decides which subject in subjects{k} is used
            % q decides which state in states{k} is used
            new_data = analyse_spectra(states{k}{q},[],[],subjects{k}{j});
            data = [data new_data];
            [~,s] = strtok(states{k}{q},'_');
            labels = [labels sprintf('%s %s (%i)',strtok(s,'_'),descriptor{k},length(new_data.scaling))];
        end
        [~,s] = strtok(states{1}{q},'_');
        make_hist_graphs(data,labels,sprintf('%s_actelion_subject_%i',strtok(s,'_'),study_id(j)));
    end
end

