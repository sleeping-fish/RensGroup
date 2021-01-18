%% this function denoise the hyperspectral Raman data by SVD and Statistical analysis
function [Reconstructed, Index_selected, num] = svds_auto(data,C,img_region, windowsize, count, thr1, thr2)
%% parameters 
% data: 2D matrix of the signal 
% C: the columns of the imaging area
% img_region : select the wavenumber region to determine the smooth area
% windowsize: the window size to search the image to find local smooth area, 5 is recommended
% count: how many smooth area were searching, 100 as default
% thr1: determine how many svd components are in considerate, 1e-3 as default
% thr2: determine which svd components is selected, 1e-2 as default
% Reconstructed: the denoised signal using svd denoising
% Index_selected: the selected SVD component index
% num: how many top svd components are searched in the current process
%% config the default parameters
switch nargin
    case 2
        img_region = 1280:1290;
        windowsize = 5;
        count = 100;
        thr1 = 1e-3;
        thr2 = 1e-2;
    case 1
        img_region = 1280:1290;
        C = 400;
        windowsize = 5;
        count = 100;
        thr1 = 1e-3;
        thr2 = 1e-2;
    case 7
        disp(' ');
    otherwise
        error("Need 7 parameters: data,columns,imaging region, window size, count, thr1, thr2, but only %d is given!\n", nargin);
end
%% find smooth local area on the msemap
[D,~] = size(data); % 以data.shape = 1340 * 33200 为例， D=1340
nim = mean(data(img_region,:),1); % 抓出成像波数的数据，求每一列的均值，得到成像区域的平均像素值
nim = reshape(nim,C,[])'; % C=400, nim.shape = (400 * 83).T = 83 * 400
if nargout~=1
    f1 = figure; figure(f1); subplot(211);imagesc(nim); colormap jet; title('Lipid map');set(gca,'XTick',[],'YTick',[]);
end
gap = windowsize -1; %gap = 5-1=4
[R,~] = size(nim); % R=83
msemap = zeros(R-gap,C-gap); %mesmap.shape=79 * 396 
for i = 1:R-gap
    for j = 1:C-gap % 遍历mesmap,用4 * 4的窗口扫描图像，将窗口内的每个值X 改写成 [(X-mean) ** 2] /(mean ** 2) 这个值啥意义？
        patch = nim(i:i+gap,j:j+gap);
        mv =  mean(mean(patch));
        patch = patch - mv;
        patch = patch.^2;
        patch = patch/(mv^2);
        msemap(i,j) = mean(mean(patch));
    end
end
% show the msemap
if nargout ~= 1
    figure(f1);subplot(212);imagesc(msemap);colormap jet; title('Mse Map');set(gca,'XTick',[],'YTick',[]);
end
vec = reshape(msemap,1,[]); % msemap展成一列
vec = sort(vec); %按从小到大排序
X = [];
Y = [];
for i = 1:count
    [x,y] = find(msemap == vec(i));
    X = [X,x(end)]; % 遍历vec,将vec中与msemap相等的元素坐标抓出来，将最大的坐标赋给X和Y
    Y = [Y,y(end)];
end
if nargout ~= 1
    figure(f1);subplot(212);hold on; scatter(Y(1:ceil(count/10):end),X(1:ceil(count/10):end),80,'red','s','filled');
end
X = X + gap/2;
Y = Y + gap/2;
% show the anchors position on the image
if nargout ~= 1
    figure(f1);subplot(211);hold on; scatter(Y(1:ceil(count/10):end),X(1:ceil(count/10):end),80,'red','s','filled');
end
% Selected  smooth area
for i = 1:length(X)
    sx(i,:) = X(i)-gap/2:X(i)+gap/2;
    sy(i,:) = Y(i)-gap/2:Y(i)+gap/2;
end
cube = reshape(data,D,C,[]); cube = permute(cube,[3 2 1]);
%%
if D<100
    M = 3;
else 
    M = 9;
end
for i = 1:length(X)
    rec_spec(i,:) = cube2spec(cube,sx(i,:),sy(i,:)); % matlabsx, sy区间的平均谱
    rec_spec(i,:) = wdenoise(rec_spec(i,:),M,'DenoisingMethod','BlockJS'); % 小波去噪
end
%% SVD and stastical analysis
% determine how many top SVD components are in consideration
[U,S,V] = svd(data,'econ');V = V';s = S(S~=0); ds = diff(log(s),2); ds = smooth(ds,9); % ds是什么
try
    [index, ~] = find(abs(ds)<=thr1); num = index(1)+2;  
catch
    num = length(s);
end

snrs = []; recon = 0; 
for i = 1:num
    tmp = s(i)*U(:,i)*V(i,:);
    recon = recon + tmp;
    tmp_cube = reshape(recon,D,C,[]);
    tmp_cube = permute(tmp_cube,[3 2 1]);
    for k = 1:length(X)
        newspec = cube2spec(tmp_cube,sx(k,:),sy(k,:));
        residual = rec_spec(k,:) - newspec';
        snrs(k,i) = snr(rec_spec(k,:),residual);
    end
end
% determine which SVD componets is selected to reconstruct the signal
% according to the mean SNRS matrix
pad = zeros(count,1); snrs = [pad, snrs]; dfsnr = [];
for i = 1:length(X)
    dfsnr(i,:) = diff(snrs(i,:));
end
mdfsnr = mean(dfsnr); [Index_selected,~]=find(mdfsnr'>=thr2);
if nargout ~= 1
    figure;subplot 121; plot(mdfsnr);
    subplot 122; plot(dfsnr');
end
ss = zeros(size(S)); 
for j = 1:length(Index_selected)
    ss(Index_selected(j),Index_selected(j))= s(Index_selected(j));
end
Reconstructed = U*ss*V;
end