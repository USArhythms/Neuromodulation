function nwb_list = f_sortNWB(path_NWB)

% path_NWB = '/projectnb/devorlab/bcraus/AnalysisCode/NWB/nwb_files_HD';

files = dir(path_NWB);
files(1:2) = [];

nwb_list = struct('Mouse',[],'Date',[],'GRAB',[],'Run',[],'iRun',[],'Path',[]);

N = numel(files);

for i = 1:N
    name = files(i).name;
    name = strsplit(name,'_');

    nwb_list(i).Mouse = strrep(name{1}(5:end),'-','_');
    nwb_list(i).Date = name{2}(5:end);
    nwb_list(i).Run = str2double(name{3}(5:6));
    nwb_list(i).iRun = str2double(name{4}(6:7));
    nwb_list(i).Path = fullfile(files(i).folder,files(i).name);
    nwb_list(i).GRAB = h5read(nwb_list(i).Path,'/general/subject/genotype');
end

end