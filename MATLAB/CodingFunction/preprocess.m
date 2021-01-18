clear;
clc;

filePath = 'D:\pythonProject\Data\FastImaging';
fileNames = scanDir(filePath); % 读取数据
count = 0;
for i = fileNames
    temp_data = importdata([filePath '\' char(i)]); % 读取数据
%     temp_despiked = despike3(temp_data, 300); % 去鬼峰
%     temp_svded = svds_auto(temp_despiked, 300); % 去噪
%     data = (temp_svded - mean2(temp_svded))/std2(temp_svded); % 标准化
%     writematrix(data, [filePath '\processed-' char(i)], 'delimiter', '\t');
    count = count + 1;
    disp(count);
end    