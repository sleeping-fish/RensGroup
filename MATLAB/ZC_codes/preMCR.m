%对MCR数据进行预分析
%扣背景
%标准化
%利用PCA确定成分
clc
close all
clear
load DataWN.txt
DataWN=DataWN(3:end,:);
[r,c]=size(DataWN);
Data=DataWN(:,2:c);
WN=DataWN(:,1);
sms=zeros(r,c-1);
for i=1:c-1
[sms(:,i)] = despike(Data(:,i),50);
end
for i=1:c-1
sms(:,i)=smooth(sms(:,i),9,'sgolay',3);%进行SG平滑处理
end
%[sub,bg]= airPLS(sms',10e7,3,0.1,0.5,20); %
%s=sub';% s是扣除背景数据
%bg=bg';
save s.mat 'sms'
ss=zscore(sms);%标准化
[COEFF,SCORE,latent,tsquare]=princomp(ss');
%COEFF：载荷 SCORE：得分
explained=100*latent/sum(latent);% 每个组分的贡献率
SumContribution=cumsum(explained);