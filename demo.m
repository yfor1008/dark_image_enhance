close all
clear
clc

[file, path] = uigetfile('*.bmp; *.jpg; *.png', 'image ...');
im = imread([path file]);
im1 = ALTM_Retinex(im);
im2 = g_channel(im);
im3 = dehaze_enhance(im);

subplot(141), imshow(im)
subplot(142), imshow(im1)
subplot(143), imshow(im2)
subplot(144), imshow(im3)

im4 = nonlinear_color_enhance(im);
figure, imshow(im4)
im5 = local_exponential_correction(im);
figure, imshow(im5)
im6 = adaptive_local_tonemap(im);
figure, imshow(im6)


