function [enhanced] = dehaze_enhance_1(im)
% dehaze_enhance_1 - matlab�Դ�ȥ�����
% 

im_inv = imcomplement(im);
im_inv1 = imreducehaze(im_inv, 'Method','approx','ContrastEnhancement','boost');
enhanced = imcomplement(im_inv1);

end

