% ���ģ� LIME: A Method for Low-light Image Enhancement
% ���ߣ�Xiaojie Guo
% ���ӣ�
% Author: HSW
% Date: 2018-04-27

clc;
close all;
clear;

addpath(genpath('removeHaze\')); 
addpath(genpath('BM3D\')); 

% img = imread('timg1.png');
img = imread('349293-20191013233603154-988940053.jpg');
% img = imread('timg3.png');
% img = imread('timg4.png');
figure;
imshow(img, []);
title('ԭͼ��');
Model = 'Normal'; % 'Dehaze' / 'Normal'
method = 'SA'; 
denoiseFlag = 0;  % = 0ʱ��ʾ�������������� = 1ʱ��ʾ������������
if size(img, 3) == 3
    img_in = im2double(img); 
    img_out = Lime_enhance(img_in, Model, method, 0, 0);     
    figure;
    imshow(img_out, []);
    title(['��ǿ���(', 'Model:', Model, ' | Smooth Method:', method, ')']);
else
    disp('LIMEģ�ʹ����ɫͼ��');
end