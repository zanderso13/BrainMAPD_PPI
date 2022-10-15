fnames = filenames('*.mat');
for sub = 1:length(fnames)
    dat = load(fnames{sub});
    onsets = {dat.onsets{3:4},dat.onsets{1:2},dat.onsets{11}};
    onsets{1} = onsets{1} - 1.9;
    onsets{2} = onsets{2} - 1.9;
    onsets{3} = onsets{3} - 1.9;
    onsets{4} = onsets{4} - 1.9;
    durations = {dat.durations{3:4},dat.durations{1:2},dat.durations{11}};
    names = {dat.names{3:4},dat.names{1:2},dat.names{11}};
    save(fullfile(strcat('anticipation/',fnames{sub})), 'names', 'durations', 'onsets')
end

clear onsets durations names

for sub = 1:length(fnames)
    dat = load(fnames{sub});
    onsets = {dat.onsets{7:8},dat.onsets{5:6},dat.onsets{11}};
    onsets{1} = onsets{1} - 1.9;
    onsets{2} = onsets{2} - 1.9;
    onsets{3} = onsets{3} - 1.9;
    onsets{4} = onsets{4} - 1.9;
    durations = {dat.durations{7:8},dat.durations{5:6},dat.durations{11}};
    names = {dat.names{7:8},dat.names{5:6},dat.names{11}};
    save(fullfile(strcat('consumption/',fnames{sub})), 'names', 'durations', 'onsets')
end