function reg = f_regImages(imgs,refParcellation,parcellationStruct,isBM)

N = numel(imgs);

reg = cell(1);

for i = 1:N
    if ~isempty(imgs{i})
        reg{i} = f_ImgReg_allen(refParcellation,parcellationStruct{i},imgs{i},isBM);
    else
        reg{i} = NaN(size(refParcellation.Masks,[1,2]));
    end
end

reg = cat(3,reg{:});

end