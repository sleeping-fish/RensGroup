clc
close all
clear
load sub.mat
num=10;
H=400;%水平pixel
V=101;%竖直pixel
%。。。。。。。。MCR处理 基于找到纯物种的 纯物质是给MCR的初始量 MCR会可以继续优化
%也可以用PCA获得对应的载荷作为初始量，可以把PCA的数量优化为非负 
%也可以用kmeans获得的个类的平均谱图作为初始量。
%PCA kmeans 合适与细胞未知体系。
[sp,imp]=pure(sub,num,10);%sp是每一行代表一张光谱
[copt,sopt,sdopt,ropt,areaopt,rtopt]=als(sub',sp,1,100,0.1);%s每一列代表一张光谱 
%。。。。。。。。MCR处理 基于EFA分析出浓度变化
%[e,eforward,ebackward]=efa(s,c-1);% s每一列代表波数 行代表强度, c-1必须等于s的行数 e每一列代表一种物质的浓度变化趋势
%[copt,sopt,sdopt,ropt,areaopt,rtopt]=als(s,e',1,50,0.1)
save csopt.mat
dlmwrite ('concentration.txt',copt,'\t')
dlmwrite ('MCRspec.txt',sopt,'\t')
for ii=1:num
pcopt=reshape(copt(:,ii),H,V);
pcopt=pcopt';
I=mat2gray(pcopt);%把矩阵变成灰度图
figure
imshow(I);
Imagecell{ii}=I;
name=['Image',num2str(ii)];
saveas(gcf,name,'jpg');
save(name,'I')
name2=['Image',num2str(ii),'.txt'];% 
dlmwrite(name2,pcopt,'\t')
end