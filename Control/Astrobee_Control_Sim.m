files = dir('Matrices/*.mat');
for i=1:length(files)
    eval(['load Matrices/' files(i).name]);
end
