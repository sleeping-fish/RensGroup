%% desipike the data cube
function despiked = despike3(data,W,windowsize,thr)
%% parameters:
% D: spectral depth
% W: width
% data: 2D matrix
% thr: threshold 
% windowsize : medfilt window size, ��3��3��3�� is recommended
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
%% ��������
D = size(data,1);  %for instance, dataCube is 1340 * 32000, D = 1340
cube = reshape(data,D,W,[]); % W �����Ǿ���ͼ��ĳ���cube.shape = 1340, 400, 80
cube = permute(cube,[3,2,1]); % cube.shape = 80, 400, 1340

%% ��ֵ�˲�
medcube = medfilt3(cube,[windowsize, windowsize, windowsize], 'symmetric'); % ��ֵ�˲��������¹��˵�ͼ��
err = cube - medcube; % ����в3ά��3ά����3ά��
%��err��ƽ���׽�������
scale = cube2spec(medcube,1:size(cube,1),1:W); % scale�ǹ��˺�cube��Ӧ��ƽ���ף�1340 * 1����

for i = 1:D
    err(:,:,i)=err(:,:,i)/scale(i); %ȡ����ά�в�Ķ�ά��Ƭ���Զ�ά��ƽ����
end
%% �ó���ֵ
ee = reshape(err,[],1); % �Ѳв�����һ��
m = quantile(ee,ceil(0.9999*length(ee))); % ȡ99.99%�ķ�λ��(�������ȡ��)
ee(ee>=m(end)|ee<=m(1)) = []; % �ѳ�����ֵ�����
pd = fitdist(ee,'Normal');
lw = icdf(pd,thr);
hw = icdf(pd,1-thr);

%% ɸ����岢�滻
for i = 1:D
    [x,y] = find(err(:,:,i)>=hw|err(:,:,i)<=lw);
    if ~isempty(x)
        fprintf('��%d����������%d�����\n',i,length(x));
        cube(x,y,i) = medcube(x,y,i);
    end
end
%% ����ԭʼ������ʽ
cube = permute(cube,[3,2,1]);
despiked = reshape(cube,D,[]);
end