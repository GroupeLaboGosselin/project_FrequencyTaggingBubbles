
function bulle=bubble(bubStd)
% a single bubble
maxHalfSize = 6*bubStd;
[y,x] = meshgrid(-maxHalfSize:maxHalfSize,-maxHalfSize:maxHalfSize);
bulle = exp(-(x.^2/bubStd^2)-(y.^2/bubStd^2));
clear x y tmp
end