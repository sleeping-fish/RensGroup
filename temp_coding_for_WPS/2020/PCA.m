clear all
clc;
% X_1 = load('BSA.mat','Xc');
% X_2 = load('lysozyme.mat','Xc');
% x = [X_1.Xc(:,1:30),X_2.Xc];
X_1 = load('BSA.mat', 'data_BSA');
X_2 = load('L.mat', 'data_Lysozome');
X_3 = load('S.mat', 'data_S');
x = [X_1.data_BSA(:,1:30),X_2.data_Lysozome];
x = [x, X_3.data_S];
[m,n] = size(x);
%% 数据中心化
x_mean = mean(x,2);
x_std = std(x,1,2);
std = ones(m,n);
mean = ones(m,n);
for i = 1:1:n
    mean(:,i) = x_mean;
    std (:,i) = x_std;
end
X = (x-mean)./std;
X=X';
%% PCA
[coeff,score,latent,tsquared,explained,mu] = pca(X);% X 的行代表样本数，列代表特征,波数
%% 画图
figure
scatter3(score(1:30,1),score(1:30,2),score(1:30,3),'r')
hold on
scatter3(score(31:60,1),score(31:60,2),score(31:60,3),'b')
hold on
scatter3(score(61:end,1),score(61:end,2),score(61:end,3),'g')
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
zlabel('3rd Principal Component')
legend('BSA','lysozyme','FontSize',14)
%% 
sum_explained = 0;
idx = 0;
while sum_explained < 95
    idx = idx + 1;
    sum_explained = sum_explained + explained(idx);
end
% idx
[U,S,V] = svd(X,'econ');
V=V';
X_SVD = U(:,1:45)*S(1:45,1:45)*V(1:45,:);
