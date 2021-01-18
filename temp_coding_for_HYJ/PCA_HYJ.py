# @Author:Xinyu_Lu
# -*- coding = 'utf-8' -*-
# @Time:2020/9/19 15:02
# @File: PCA_WPS.py
# @Software: PyCharm

import numpy as np
import matplotlib.pyplot as plt

# 电位顺序：0.1 0.2 0.3 -0.1 -0.2 -0.3 0
def loadData():
    wavenumber = np.loadtxt(r'D:\pythonProject\Data\氧还原时间序列\wavenumber.txt')
    data_O2 = np.loadtxt(r'D:\pythonProject\Data\氧还原时间序列\O2\processed-O2.txt')
    data_N2 = np.loadtxt(r'D:\pythonProject\Data\氧还原时间序列\N2\processed-N2.txt')
    data = np.concatenate((data_O2, data_N2), axis=1)  # O2 和 N2的矩阵组成为一个大矩阵,一起PCA
    return data, wavenumber  # (575, 700)


def SVDAnalysis(data):

    sval_nums = 3
    U, Sigma, VT = np.linalg.svd(data)
    restruct = np.diag(Sigma[0:sval_nums]) @ VT[0:sval_nums, :]
    restruct = (restruct - restruct.mean()) / restruct.std()  # 归一化
    O2_restruct = restruct[:, 0:700]
    N2_restruct = restruct[:, 700:]
    return O2_restruct, N2_restruct, Sigma


def back2wave(data):
    U, Sigma, VT = np.linalg.svd(data)
    restruct = np.dot(U[:, 0].reshape(-1, 1) * Sigma[0], VT[0, :].reshape(1, -1))
    O2_restruct = restruct[:, 0:700]
    N2_restruct = restruct[:, 700:]
    return O2_restruct, N2_restruct


def plotFigure(O2_restruct, N2_restruct, Sigma, data='', wavenumber=''):

    tempO2_0 = O2_restruct[0, :]
    tempO2_1 = O2_restruct[1, :]
    tempO2_2 = O2_restruct[2, :]
    tempN2_0 = N2_restruct[0, :]
    tempN2_1 = N2_restruct[1, :]
    tempN2_2 = N2_restruct[2, :]
    # 画图1 3D SVD
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(tempO2_0[tempO2_0 > -40][:650], tempO2_1[tempO2_1 > -40][:650], tempO2_2[tempO2_2 > -40][:650],
               marker='x', color='blue', s=40, label='class 1')
    ax.scatter(tempN2_0[tempN2_0 > -40][:650], tempN2_1[tempN2_1 > -40][:650], tempN2_2[tempN2_2 > -40][:650],
               marker='o', color='green', s=40, label='class 2')

    ax.set_xlabel('1st Component')
    ax.set_ylabel('2nd Component')
    plt.title('3D PCA Analysis')
    plt.legend(['${O_2}$', '${N_2}$'])

    # 画图2 weights
    fig2 = plt.figure(figsize=(6, 6))
    xList = np.arange(1, 6)
    yList = Sigma[:5] / Sigma.sum() * 100
    plt.plot(xList, yList)
    plt.scatter(xList, yList)
    plt.title('Weight Analysis')
    plt.xticks(np.arange(1, 6))
    plt.xlabel('No.')
    plt.ylabel('Weights')
    for x, y in zip(xList, yList):
        plt.text(x, y + 0.3, '%.1f' % y + '%', ha='center', va='bottom', fontsize=10.5)

    # 画图3 二维SVD
    fig3 = plt.figure(figsize=(6, 6))
    plt.scatter(tempO2_0[tempO2_0 > -40][:650], tempO2_1[tempO2_1 > -40][:650], color='blue', label='class 1')
    plt.scatter(tempN2_0[tempN2_0 > -40][:650], tempN2_1[tempN2_1 > -40][:650], color='red', label='class 2')
    plt.title('2D PCA Analysis')
    plt.legend(['${O_2}$', '${N_2}$'])

    # fig4 = plt.figure()
    # O2_restruct2, N2_restruct2 = back2wave(data)
    # O2Max = np.max(O2_restruct2)
    # N2Max = np.max(N2_restruct2)
    # min_max_scaler = preprocessing.MinMaxScaler()
    # O2_temp = (O2_restruct2/O2Max).sum(1).reshape(-1, 1)
    # N2_temp = (N2_restruct2/N2Max).sum(1).reshape(-1, 1)
    # O2_minmax = min_max_scaler.fit_transform(O2_temp)
    # N2_minmax = min_max_scaler.fit_transform(N2_temp)
    # # MSE_minmax = min_max_scaler.fit_transform(((O2_restruct2.sum(1) - N2_restruct2.sum(1)) ** 4).reshape(-1, 1))
    #
    # plt.plot(wavenumber, O2_minmax, label='class 1')
    # plt.plot(wavenumber, N2_minmax, label='class 2')
    #
    # plt.title('2D PCA Analysis')
    # plt.legend(['${O_2}$', '${N_2}$', 'MSE'])
    # # plt.yticks([])
    #
    # plt.figure()
    # MSE = O2_minmax - N2_minmax
    # temp = scipy.signal.savgol_filter(MSE.flatten(), 55, 2)
    # # plt.plot(wavenumber, MSE)
    # plt.plot(wavenumber, temp)
    plt.show()


def Score(O2_restruct, N2_restruct):
    score_O2 = np.abs(O2_restruct.mean(1)[0] / O2_restruct.mean(1)[1])
    score_N2 = np.abs(N2_restruct.mean(1)[0] / N2_restruct.mean(1)[1])
    print(score_O2, score_N2)


def main():
    data, wavenumber = loadData()
    O2_restruct, N2_restruct, Sigma = SVDAnalysis(data)
    Score(O2_restruct, N2_restruct)
    plotFigure(O2_restruct, N2_restruct, Sigma, data, wavenumber)


if __name__ == '__main__':
    main()
