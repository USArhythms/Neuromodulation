function nwb_list = f_sortNWB(path_NWB)

% path_NWB = '/projectnb/devorlab/bcraus/AnalysisCode/NWB/nwb_files_HD';

files = dir(path_NWB);
files(1:2) = [];

nwb_list = struct('Mouse',[],'Date',[],'GRAB',[],'Run',[],'Path',[]);

N = numel(files);

for i = 1:N
    name = fullfile(files(i).folder,files(i).name);
    nwb_list(i).Path = name;
    nwb_list(i).GRAB = h5read(name,'/general/subject/genotype');
    name = h5read(name,'/identifier');
    name = strsplit(name,'/');

    nwb_list(i).Mouse = name{1};
    nwb_list(i).Date = name{2};
    nwb_list(i).Run = str2double(name{3}(4:5));
end

end