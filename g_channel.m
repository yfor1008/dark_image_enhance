function [out] = g_channel(im)
% g_channel - ��ͼ����ǿ
%   https://blog.csdn.net/grafx/article/details/53233508
% 1. ����ɫͨ����ɫ����Ϊϵ��ֵ, �ֱ������ͨ����ˣ��õ���ͼ��
% 2. ����ͼ����ԭͼ��һ����ɫ���, f(a, b) = 1 - (1 - a)*(1 - b)
% 3. ��ν�����������, ֪��ͼ�����ȴﵽ����[90,120];
%

[h,w,c] = size(im);
out = zeros(h,w,c);

im = double(im);
r = im(:,:,1);
g = im(:,:,2);
b = im(:,:,3);

gray = 0.2989 * r + 0.5870 * g + 0.1140 * b;
V = mean(gray(:));
cnt = 0;
while (V < 80 && cnt < 5)
    
    g_alpha = 255 - g;

    r1 = r .* g_alpha / 255;
    g1 = g .* g_alpha / 255;
    b1 = b .* g_alpha / 255;

    r = 255 - (255 - r) .* (255 - r1) / 255;
    g = 255 - (255 - g) .* (255 - g1) / 255;
    b = 255 - (255 - b) .* (255 - b1) / 255;
    
    gray = 0.2989 * r + 0.5870 * g + 0.1140 * b;
    V = mean(gray(:));
    
    cnt = cnt + 1;
end

% g_alpha = 255 - g;
% 
% r1 = r .* g_alpha / 255;
% g1 = g .* g_alpha / 255;
% b1 = b .* g_alpha / 255;
% 
% r = 255 - (255 - r) .* (255 - r1) / 255;
% g = 255 - (255 - g) .* (255 - g1) / 255;
% b = 255 - (255 - b) .* (255 - b1) / 255;

out(:,:,1) = r;
out(:,:,2) = g;
out(:,:,3) = b;

out = uint8(out);
end

