# @Author:Xinyu_Lu
# -*- coding = 'utf-8' -*-
# @Time:2020/12/15 23:01
# @File: germs_classicication_by_ML_predict.py
# @Software: PyCharm

import numpy as np
import joblib

'''超参数'''
class_num = 6
germs = ['Chryseo', 'E25922', 'EF29212', 'PA27853', 'SA29213', 'Salmonella']
class_dict = dict()
for i in range(class_num):
    class_dict[str(i)] = germs[i]

def load_data(label=1,
        dataset_path=r'D:\pythonProject\Python_code\temp_coding_for_CuiLi\temp_dataset.txt'):
    print('【加载数据】')
    dataset = np.loadtxt(dataset_path)
    print('共有{}条测试数据'.format(dataset.shape[0]))
    if label:
        inputs = dataset[:, :-1]
        labels = dataset[:, -1]
        return inputs, labels
    else:
        return dataset


def lda_predict(test_X, test_y=None, label=0):
    print('【LDA预测】')
    lda = joblib.load('./lda_temp.pkl')
    pred = lda.predict(test_X)
    pred_labels = []
    for item in pred:
        item = str(item.astype(int))
        pred_labels.append(class_dict[item])
    print('预测结果如下：')
    print(pred_labels)
    if label:
        print('和标签对比，输出准确率……')
        for j in range(class_num):
            correct = 0
            total = 0
            for i in range(test_y.shape[0]):
                if test_y[i] == j:
                    total += 1
                    if pred[i] == test_y[i]:
                        correct += 1
            print('对{}的分类准确率为：{:.2f}%'.format(germs[j], correct / total * 100))
        accuracy_lda = lda.score(test_X, test_y)
        print('LDA分类的平均准确率为：{:.2f}%'.format(accuracy_lda * 100))


def svm_predict(test_X, test_y=None, label=0):
    print('【SVM预测】')
    clf = joblib.load('./svm_temp.pkl')
    pred = clf.predict(test_X)
    pred_labels = []
    for item in pred:
        item = str(item.astype(int))
        pred_labels.append(class_dict[item])
    print('预测结果如下：')
    print(pred_labels)
    if label:
        print('和标签对比，输出准确率……')
        for j in range(class_num):
            correct = 0
            total = 0
            for i in range(test_y.shape[0]):
                if test_y[i] == j:
                    total += 1
                    if pred[i] == test_y[i]:
                        correct += 1
            print('对{}的分类准确率为：{:.2f}%'.format(germs[j], correct / total * 100))
        accuracy_svm = clf.score(test_X, test_y)
        print('SVM分类的平均准确率为：{:.2f}%'.format(accuracy_svm * 100))


def main():
    test_X, test_y = load_data()
    lda_predict(test_X, test_y, label=0)
    svm_predict(test_X, test_y, label=1)


if __name__ == '__main__':
    main()
