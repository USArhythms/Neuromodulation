function [regImg,regFactor] = f_ImgReg_allen(refParcellation,parcellationStruct,img,isBM)
%%
% parcells1 = refParcellation;
% parcells2 = comb.allen{1};
% img = comb.vis{1};

if 
L = sum(parcellationStruct.Masks(:,:,:,1),'all');
R = sum(parcellationStruct.Masks(:,:,:,2),'all');

if L && R
    hem = 2;
    side = 0;
else
    hem = 1;
    if L
        side = 0;
    else
        side = 1;
    end
end

%% find points

if hem == 1
    Vis2 = parcellationStruct.Masks(:,:,12,side+1);
    SSp_ll2 = parcellationStruct.Masks(:,:,5,side+1);
    points.p2.point1 = round(findTop(Vis2));
    points.p2.point2 = round(findTop(SSp_ll2));
    
    Vis1 = refParcellation.Masks(:,:,12,side+1);
    SSp_ll1 = refParcellation.Masks(:,:,5,side+1);
    points.p1.point1 = round(findTop(Vis1));
    points.p1.point2 = round(findTop(SSp_ll1));

else
    Vis2L = parcellationStruct.Masks(:,:,12,1);
    Vis2R = parcellationStruct.Masks(:,:,12,2);
    points.p2.point1 = round(findTop(Vis2L));
    points.p2.point2 = round(findTop(Vis2R));

    Vis1L = refParcellation.Masks(:,:,12,1);
    Vis1R = refParcellation.Masks(:,:,12,2);
    points.p1.point1 = round(findTop(Vis1L));
    points.p1.point2 = round(findTop(Vis1R));

end

%% adjust tilt

regFactor.dim1 = size(refParcellation.Masks(:,:,1,1));

tilt(1) = atan((points.p1.point1(1)-points.p1.point2(1))/(points.p1.point1(2)-points.p1.point2(2)));
tilt(2) = atan((points.p2.point1(1)-points.p2.point2(1))/(points.p2.point1(2)-points.p2.point2(2)));

d_tilt = tilt(1)-tilt(2);

dim = size(parcellationStruct.Masks(:,:,1,1));
% centroid = dim/2+0.5;
% img_padding = centroid-points.p2.point1;
% img_padding = zeros(2*abs(img_padding)+dim);
% 
% centroid = size(img_padding)/2+0.5;
% 
% c1 = centroid-points.p2.point1+1;
% c2 = c1+dim-1;
% 
% error('test')
% 
% img_padding(c1(1):c2(1),c1(2):c2(2)) = img;
% img_padding = imrotate(img_padding,d_tilt*180/pi);
% 
regFactor.tilt = d_tilt;
% img = img_padding(c1(1):c2(1),c1(2):c2(2));

%% ChatGPT rotate

cx = points.p2.point1(2);
cy = points.p2.point1(1);

angle = d_tilt*180/pi;
T1 = [1 0 -cx; 0 1 -cy; 0 0 1];
R = [cosd(angle) -sind(angle) 0;sind(angle) cosd(angle) 0; 0 0 1];
T2 = [1 0 cx;0 1 cy;0 0 1];

T = T2*R*T1;

tform = affine2d(T');
img = imwarp(img,tform,'OutputView',imref2d(size(img)));

%% adjust scale

dist(1) = sqrt(sum((points.p1.point1-points.p1.point2).^2));
dist(2) = sqrt(sum((points.p2.point1-points.p2.point2).^2));

scale = dist(1)/dist(2);

xg = 1:dim(1);
yg = 1:dim(2);
F = griddedInterpolant({xg,yg},img);

xq = (1/scale:1/scale:dim(1))';
yq = (1/scale:1/scale:dim(2))';
vq = F({xq,yq});

regFactor.scale = scale;
regFactor.dim2 = size(vq);

scaledPoint = round(points.p2.point1*scale);

%% align image sizes

regImg = NaN(regFactor.dim1);

dim = size(vq);
c1 = max([1,1;points.p1.point1-scaledPoint+1]);
c2 = min([regFactor.dim1;points.p1.point1-scaledPoint+1+dim-1]);

reg_c2 = c2-points.p1.point1+scaledPoint;
reg_c1 = scaledPoint-points.p1.point1+c1;

regImg(c1(1):c2(1),c1(2):c2(2)) = vq(reg_c1(1):reg_c2(1),reg_c1(2):reg_c2(2));

if side == 0
    regImg = fliplr(regImg);
end

if nargin == 4 && isBM
    regImg(regImg==0) = NaN;
end

%% support functions
function [point] = findTop(img)
    point = find(sum(img,2));
    point = point(1);
    point(2) = mean(find(img(point,:)));
end

end