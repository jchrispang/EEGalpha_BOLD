function parse_data(track_data, subject_num, sleep_state)

data = track_data((track_data.subject==subject_num & track_data.state==sleep_state), :);

save(['data/Tracking_state=EC_subject=', num2str(subject_num), '.mat'], 'data')