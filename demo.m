close all
clear
clc

% [file, path] = uigetfile('*.bmp; *.jpg; *.png', 'image ...');
% im = imread([path file]);
% 
% % im1 = ALTM_Retinex(im);
% % im2 = g_channel(im);
% im3 = dehaze_enhance(im);
% im31 = dehaze_enhance_1(im);
% % 
% % subplot(141), imshow(im)
% % subplot(142), imshow(im1)
% % subplot(143), imshow(im2)
% % subplot(144), imshow(im3)
% 
% figure, imshow(im)
% figure, imshow(im3);
% figure, imshow(im31);
% % 
% % im4 = nonlinear_color_enhance(im);
% % figure, imshow(im4)
% % im5 = local_exponential_correction(im);
% % figure, imshow(im5)
% % im6 = adaptive_local_tonemap(im);
% % figure, imshow(im6)
% % im7 = homo_filter(im);
% % figure, imshow(im7)
% % im8 = neighbor_nonlinear_enhance(im);
% % figure, imshow(im8)
% % im9 = Ying_2017_ICCV(im);
% % figure, imshow(im9)
% % im10 = Ying_2017_CAIP(im);
% % figure, imshow(im10)
% % im11 = multi_fusion(im);
% % figure, imshow(im11)
% % [im12, R, L] = jed(im);
% % figure, imshow(im12)
% % im13 = multi_fusion_1(im);
% % figure, imshow(im13)


image_folder = 'D:\mywork\work\test_and_result\2020_02_temp\ankon_samples\';
img_type = '*.jpg';
imgs = dir([image_folder img_type]);
img_len = length(imgs);

for idx = 1 : img_len
    im = imread([image_folder imgs(idx).name]);
%     im(:,:,1) = medfilt2(im(:,:,1));
%     im(:,:,2) = medfilt2(im(:,:,2));
%     im(:,:,3) = medfilt2(im(:,:,3));
%     
    im1 = imlocalbrighten(im);
    
    imwrite(im1, ['.\ankon_samples_nne\' imgs(idx).name]);

end
