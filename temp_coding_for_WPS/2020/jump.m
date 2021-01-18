
%1000 1200 1350 1450

clear all;
clc;
alldata = importdata('D:\pythonProject\Data\Data4WPS\SARS_antigen\point 5\SARS 1mgml gap mode-633-0.8Na-0.5%-10s-point5-time.txt');
seq = 50;%谱图数量
[m,~] = size(alldata);
wavenumber = alldata(1:m/seq,2);
intensity = alldata(:,3);
intensity = reshape(intensity,m/seq,[]);
%% 加上之前的测试谱(如果需要加上采的单个谱直接去掉所有前缀%)
%filePath = 'D:\MATLAB\raman protein'; % 这里输入数据文件的文件夹
%fileFolder=fullfile(filePath);
%for item = 1:1:5%往前累加图谱
%    path = [num2str(item), '.txt'];
%   temp = importdata(path);
%    x = temp(:,2);
%   intensity = [x,intensity];
%end
%% 扣鬼峰(有点小问题)
% intensity = despike(despike(intensity,50,1e-4),5,1e-4);
X = intensity;  % x 原始光谱
%% airPLS 扣背景
% [Xc,Z]=airPLS(intensity',10e7,3,0.1,0.5,20);
% X = Xc;
%% 去掉17 光谱(去掉异常光谱，如果需要加上采的单个谱直接去掉所有前缀%))
%figure
%plot(wavenumber,Xc(17,:)');axis tight;grid on;title('No. 17');
%X = [Xc(1:16,:);Xc(18:end,:)];
%% visualization
y = 1:1:seq;
figure
% imagesc(wavenumber',y,X);colormap jet; colorbar;title('SARS antigen corrected by airPLS');
waterfall(wavenumber',y,X');colormap jet; colorbar;title('SARS antigen corrected by airPLS');