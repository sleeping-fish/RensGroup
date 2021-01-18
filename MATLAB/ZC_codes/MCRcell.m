clc
close all
clear
load sub.mat
num=10;
H=400;%ˮƽpixel
V=101;%��ֱpixel
%����������������MCR���� �����ҵ������ֵ� �������Ǹ�MCR�ĳ�ʼ�� MCR����Լ����Ż�
%Ҳ������PCA��ö�Ӧ���غ���Ϊ��ʼ�������԰�PCA�������Ż�Ϊ�Ǹ� 
%Ҳ������kmeans��õĸ����ƽ����ͼ��Ϊ��ʼ����
%PCA kmeans ������ϸ��δ֪��ϵ��
[sp,imp]=pure(sub,num,10);%sp��ÿһ�д���һ�Ź���
[copt,sopt,sdopt,ropt,areaopt,rtopt]=als(sub',sp,1,100,0.1);%sÿһ�д���һ�Ź��� 
%����������������MCR���� ����EFA������Ũ�ȱ仯
%[e,eforward,ebackward]=efa(s,c-1);% sÿһ�д����� �д���ǿ��, c-1�������s������ eÿһ�д���һ�����ʵ�Ũ�ȱ仯����
%[copt,sopt,sdopt,ropt,areaopt,rtopt]=als(s,e',1,50,0.1)
save csopt.mat
dlmwrite ('concentration.txt',copt,'\t')
dlmwrite ('MCRspec.txt',sopt,'\t')
for ii=1:num
pcopt=reshape(copt(:,ii),H,V);
pcopt=pcopt';
I=mat2gray(pcopt);%�Ѿ����ɻҶ�ͼ
figure
imshow(I);
Imagecell{ii}=I;
name=['Image',num2str(ii)];
saveas(gcf,name,'jpg');
save(name,'I')
name2=['Image',num2str(ii),'.txt'];% 
dlmwrite(name2,pcopt,'\t')
end