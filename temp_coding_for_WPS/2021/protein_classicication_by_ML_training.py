# @Author:Xinyu_Lu
# -*- coding = 'utf-8' -*-
# @Time:2020/12/15 23:01
# @File: protein_classicication_by_ML_training.py
# @Software: PyCharm

import numpy as np
import matplotlib.pyplot as plt
import os

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn import metrics
from sklearn import svm
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import cross_val_score
import joblib

'''超参数'''
total_line_per_protein = 31
class_num = 2
feature_per_protein = 1340
proteins = ['HSA', 'BSA']  # HSA: 0  31 BSA: 1 42
marker = ['8', 'x']

'''导入数据'''


def LoadData(
        load=1, std=1, minmax=1, plot=1,
        dataset_path=r'D:\pythonProject\Data\WPS\20200106\data\data_set.txt',
        dir_path=r'D:\pythonProject\Data\CL\20201114 pathogen classification',
        wave_path=r'D:\pythonProject\Data\CL\20201114 pathogen classification\wave.txt',
        save_name='temp_dataset.txt', ):
    features_set = [0, 31, 42]
    if load:
        print('【读取数据】\n读取{}文件夹下的数据……'.format(dir_path))
        dataset = np.zeros((1, feature_per_protein + 1))  # 1 * 526+1
        for i in range(class_num):
            for _, _, files in os.walk(dir_path + '\\' + proteins[i]):
                print('读取{},数据量为：{}条'.format(proteins[i], len(files)))
                temp2 = np.zeros((len(files), feature_per_protein + 1))  # [] * 526+1
                for j in range(len(files)):
                    temp = np.loadtxt(dir_path + '\\' + proteins[i] + '\\' + files[j])
                    temp = temp[:feature_per_protein, 1].reshape(-1, 1)  # 526 * 1
                    if std:
                        std_scaler = StandardScaler()
                        temp = std_scaler.fit_transform(temp)
                    if minmax:
                        minmax_scaler = MinMaxScaler()
                        temp = minmax_scaler.fit_transform(temp)
                    temp2[j, :-1] = temp.flatten()
                    temp2[j, -1] = i
                dataset = np.concatenate((dataset, temp2), axis=0)
        dataset = dataset[1:, :]
        if std:
            print('对数据进行Z-Score标准化……')
        if minmax:
            print('对数据进行Min-Max归一化……')
        print('合并数据集完成，将其保存在{}'.format('./' + save_name))
        # 注意：dataset的最后一列是labels，除最后一列均为inputs
        np.savetxt('./' + save_name, dataset)
    else:
        dataset = np.loadtxt(dataset_path)
    labels = dataset[:, -1]
    inputs = dataset[:, :-1]

    if plot:
        wave = np.loadtxt(wave_path)
        plt.figure(figsize=(60, 10))
        print('【初始数据可视化】\n{}类的平均光谱如图所示'.format(class_num))
        for i in range(class_num):
            ax = plt.subplot(class_num, 1, i + 1)
            plt.plot(wave, inputs[i * total_line_per_protein:(i + 1) * total_line_per_protein, :].mean(0))
            plt.ylabel('a.u.')
            plt.xlabel('wavenumber($\mathregular{cm^{-1}}$)')
            plt.title('{}'.format(proteins[i]))
        plt.show()
    return inputs, labels, features_set


'划分训练集和测试集'


def Split(inputs, labels, test_size=0.5):
    train_X, test_X, train_y, test_y = train_test_split(
        inputs, labels, test_size=test_size, random_state=0, stratify=labels)

    print('【划分训练集和测试集】\n总数据量：{}条，训练集：{}条，测试集：{}条'
          .format(inputs.shape[0], train_X.shape[0], test_X.shape[0]))
    print('【建立分析模型】')
    return train_X, test_X, train_y, test_y


'PCA'


def pca_anal(train_X, test_X, features_set, plot=1, returns=0):
    print('PCA……')
    pca = PCA(n_components=0.99)
    train_X_PCA = pca.fit_transform(train_X)
    # plot PCA
    if plot:
        plt.figure()
        for i in range(class_num):
            plt.scatter(
                train_X_PCA[int(features_set[i]):int(features_set[i + 1]), 0],
                train_X_PCA[int(features_set[i]):int(features_set[i + 1]), 1],
                marker=marker[i], label=proteins[i])
        plt.legend()
        plt.title('PCA Classification')
    plt.show()
    # return results of PCA(recommend not)
    if returns:
        train_X = pca.inverse_transform(train_X_PCA)
        test_X_PCA = pca.fit_transform(test_X)
        test_X = pca.inverse_transform(test_X_PCA)

        return train_X, test_X


'''LDA'''


def lda_anal(train_X, test_X, train_y, test_y, plot=1, save_path='./lda_temp.pkl'):
    print('LDA训练中……')
    # create model
    lda = LDA()
    lda.fit(train_X, train_y)
    y_score = lda.decision_function(test_X)
    # CV
    CV(lda, train_X, train_y, test_X, test_y)
    # ROC
    PlotROC(test_y, y_score, 'LDA')
    # save model
    print('保存模型')
    joblib.dump(lda, save_path)
    # plot LDA
    if plot:
        plt.figure()
        ordered_X = np.array([])
        ordered_y = []
        label_num = []
        for i in range(class_num):
            index = []
            for j, k in enumerate(train_y):
                if k == i:
                    index.append(j)
            ordered_X = np.append(ordered_X, train_X[index])
            ordered_y = np.append(ordered_y, train_y[index])
            label_num.append(len(index))
        label_num.insert(0, 0)
        label_num = np.array(label_num)
        train_X_LDA = lda.fit_transform(ordered_X.reshape(-1, feature_per_protein), ordered_y.reshape(-1, 1))
        print(train_X_LDA.shape)
        for l in range(class_num):
            plt.scatter(
                range(label_num[l+1]),
                train_X_LDA[label_num[l]:label_num[l] + label_num[l + 1]],
                marker=marker[l], label=proteins[l])
        plt.legend()
        plt.title('LDA Classification')
        plt.show()


'''SVM'''


def svm_anal(train_X, test_X, train_y, test_y, save_path='./svm_temp.pkl'):
    print('SVM训练中……')
    # create model
    clf = svm.SVC(C=5, kernel='linear', gamma=20, decision_function_shape='ovr')
    clf.fit(train_X, train_y)
    y_score = clf.decision_function(test_X)
    # CV
    CV(clf, train_X, train_y, test_X, test_y)
    # ROC
    PlotROC(test_y, y_score, 'SVM')
    # save model
    print('保存模型')
    joblib.dump(clf, save_path)


def CV(clf, train_X, train_y, test_X, test_y):
    print('【交叉验证】')
    # 10-fold test
    score = cross_val_score(clf, train_X, train_y, cv=10)
    print('十折交叉验证的分数分别是：\n{}'.format(score))
    pred = clf.predict(test_X)
    print('交叉验证报告如下：')
    print(metrics.classification_report(test_y, pred))
    cm = metrics.confusion_matrix(test_y, pred)
    print('混淆矩阵：')
    print(cm)
    plt.figure()
    print('混淆矩阵可视化：')
    plt.imshow(cm, cmap='binary')
    plt.colorbar()
    plt.xticks(range(class_num), [proteins[i] for i in range(class_num)])
    plt.yticks(range(class_num), [proteins[i] for i in range(class_num)])
    plt.title('confusion_matrix')
    plt.show()
    print('对不同物种的分类结果：')
    for j in range(class_num):
        correct = 0
        total = 0
        for i in range(test_y.shape[0]):
            if test_y[i] == j:
                total += 1
                if pred[i] == test_y[i]:
                    correct += 1
        print('对{}的分类准确率为：{:.2f}%'.format(proteins[j], correct / total * 100))
    accuracy_svm = clf.score(test_X, test_y)
    print('分类的平均准确率为：{:.2f}%'.format(accuracy_svm * 100))


'''绘制 ROC 曲线'''


def PlotROC(test_y, y_score, title):
    # 计算相关参数
    test_y_one_hot = label_binarize(test_y, classes=[i for i in range(class_num)])
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    y_score = y_score.reshape(-1, 1)
    for i in range(class_num - 1):
        fpr[i], tpr[i], _ = roc_curve(test_y_one_hot[:, i], y_score[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])

    # 可视化
    print('【结果可视化】\nROC曲线如图所示')
    plt.figure()
    for i in range(class_num - 1):
        plt.plot(fpr[i], tpr[i],
                 label='ROC curve of %s (area = %0.4f)' % (proteins[i], roc_auc[i]))
        plt.plot([0, 1], [0, 1], color='navy', linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(title)
        plt.legend(loc="lower right")
    plt.show()


def main():
    load = 0  # 是否需要重新加载原始数据
    std = 1  # 是否需要Z-Score标准化
    minmax = 1  # 是否需要Min-Max 标准化
    plot = 0  # 是否需要绘制谱图
    test_size = 0.2  # 测试集划分率
    pca_returns = 0
    inputs, labels, features_set = LoadData(load=load, std=std, minmax=minmax, plot=plot)
    train_X, test_X, train_y, test_y = Split(inputs, labels, test_size=test_size)

    if pca_returns:
        train_X, test_X = pca_anal(train_X, test_X, features_set, plot=1, returns=pca_returns)
    else:
        pca_anal(inputs, test_X, features_set, plot=1, returns=pca_returns)

    lda_anal(train_X, test_X, train_y, test_y, plot=1)
    svm_anal(train_X, test_X, train_y, test_y)


if __name__ == '__main__':
    main()
