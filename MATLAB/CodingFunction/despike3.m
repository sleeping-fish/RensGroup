%% desipike the data cube
function despiked = despike3(data,W,windowsize,thr)
%% parameters:
% D: spectral depth
% W: width
% data: 2D matrix
% thr: threshold 
% windowsize : medfilt window size, （3，3，3） is recommended
%% config parameters:
switch nargin
    case 1
        W = 400;
        windowsize = 3;
        thr = 1e-10;
    case 2
        windowsize = 3;
        thr = 1e-10;
    case 3
        error("4 parameters is needed, data, column pixels, window size, thr, but %d is given", nargin);
end
%% 数据整理
D = size(data,1);  %for instance, dataCube is 1340 * 32000, D = 1340
cube = reshape(data,D,W,[]); % W 这里是矩形图像的长，cube.shape = 1340, 400, 80
cube = permute(cube,[3,2,1]); % cube.shape = 80, 400, 1340

%% 中值滤波
medcube = medfilt3(cube,[windowsize, windowsize, windowsize], 'symmetric'); % 中值滤波，并存下过滤的图像
err = cube - medcube; % 计算残差（3维减3维还是3维）
%对err按平均谱进行缩放
scale = cube2spec(medcube,1:size(cube,1),1:W); % scale是过滤后cube对应的平均谱（1340 * 1？）

for i = 1:D
    err(:,:,i)=err(:,:,i)/scale(i); %取出三维残差的二维切片除以二维的平均谱
end
%% 得出阈值
ee = reshape(err,[],1); % 把残差拉成一条
m = quantile(ee,ceil(0.9999*length(ee))); % 取99.99%的分位数(向无穷大处取整)
ee(ee>=m(end)|ee<=m(1)) = []; % 把超过阈值的清空
pd = fitdist(ee,'Normal');
lw = icdf(pd,thr);
hw = icdf(pd,1-thr);

%% 筛出鬼峰并替换
for i = 1:D
    [x,y] = find(err(:,:,i)>=hw|err(:,:,i)<=lw);
    if ~isempty(x)
        fprintf('第%d个波长共有%d个鬼峰\n',i,length(x));
        cube(x,y,i) = medcube(x,y,i);
    end
end
%% 返回原始数据形式
cube = permute(cube,[3,2,1]);
despiked = reshape(cube,D,[]);
end