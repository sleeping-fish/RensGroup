%% desipike the data cube
function despiked = despike2(data,D,windowsize,thr)
%% parameters:
% D: spectral depth(wavenumber)
% data: 2D matrix (sample*feature)
% thr: threshold
% windowsize : medfilt window size, (3��3) is recommended
%% config parameters:
switch nargin
    case 2
        windowsize = 3;
        thr = 1e-10;
       case 3
        error("3 parameters is needed, data, column pixels, window size, thr, but %d is given", nargin);
end
%%
medcube = medfilt2(data,[windowsize, windowsize], 'symmetric');
err = data - medcube;
%% ����
scale = mean(medcube,1);
for i = 1:1:D
    err(:,i)=err(:,i)/scale(i);
end
%%
ee = reshape(err,[],1); %չ��һ������
m = quantile(ee,ceil(0.9999*length(ee)));
ee(ee>=m(end)|ee<=m(1)) = [];
pd = fitdist(ee,'Normal');
lw = icdf(pd,thr); %Inverse cumulative distribution function ���ۻ��ֲ�����
hw = icdf(pd,1-thr);
%%
for i = 1:1:D
    [x] = find(err(:,i)>=hw|err(:,i)<=lw);
    if ~isempty(x)
        fprintf('��%d����������%d�����\n',i,length(x));
        data(x,i) = medcube(x,i);
    end
end
despiked = data;
end
