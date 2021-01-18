function spectrum = cube2spec(cube,xx,yy,zz)
% this function change the 3d data to one dimensional spectrum
% if xx yy are or one of them is range variable, the result return the average
% spectrum in the region
% if xx yy are number, the result return the single spectrum indicated by
% xx and yy
% zz is the selected range of the spectrum
if nargin == 3
    zz = 1:size(cube,3); % 以despike3 cube.shape = 80, 400, 1340为例， zz = 1:1340(1340项的等差数列)
end

if length(xx) == 1 && length(yy) == 1 % 如果cube.shape = 1,1, zz
    spectrum = squeeze(cube(xx,yy,zz)); % 将维数为1的洗去
else 
    spectrum = mean(mean(cube(xx,yy,zz),1),2);
    spectrum = squeeze(spectrum); % 输出平均谱
end
end