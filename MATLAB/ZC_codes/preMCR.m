%��MCR���ݽ���Ԥ����
%�۱���
%��׼��
%����PCAȷ���ɷ�
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
sms(:,i)=smooth(sms(:,i),9,'sgolay',3);%����SGƽ������
end
%[sub,bg]= airPLS(sms',10e7,3,0.1,0.5,20); %
%s=sub';% s�ǿ۳���������
%bg=bg';
save s.mat 'sms'
ss=zscore(sms);%��׼��
[COEFF,SCORE,latent,tsquare]=princomp(ss');
%COEFF���غ� SCORE���÷�
explained=100*latent/sum(latent);% ÿ����ֵĹ�����
SumContribution=cumsum(explained);