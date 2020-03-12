function [enhanced] = nonlinear_color_enhance(im)
% nonlinear_color_enhance - 自适应非线性图像增强
%
% input:
%   - im: h*w*3, rgb图像
% output:
%   - enhanced: h*w*3, 增强后图像
%
% docs:
%   - 《Adaptive and integrated neighborhood-dependent approach for nonlinear enhancement of color images》
%   - 自适应亮度增强: 提高暗像素的亮度, 增强了亮度却压缩了图像的动态范围, 降低了对比度
%   - 自适应对比度增强: 根据像素值本身和邻域内平均值的比例, 比例>1时, 增大像素值, 比例<1时, 减小像素值
%   - https://www.cnblogs.com/Imageshop/p/11665100.html
%

luminance = rgb2gray(im); % 公式1
[h, w] = size(luminance);
pix_num = h * w;

im = double(im);
enhanced = zeros(h, w, 3);

% 计算直方图
histgram = zeros(256, 1);
for index = 1:pix_num
    histgram(luminance(index)+1) = histgram(luminance(index)+1) + 1;
end
im_std = std(double(luminance(:)));

% where L is the intensity level corresponding to a cumulative distribution function(CDF) of 0.1
cdf = 0;
L = 1;
ratio = 0.1 * pix_num;
for index = 1:256
    cdf = cdf + histgram(index);
    if cdf >= ratio
        L = index;
        break;
    end
end

% 计算Z
low = 50;
high = 150;
if L <= low
    % 图像中%10的像素值<low, 图像偏暗, 需增强处理  
    Z = 0;
elseif L <= high
    % 介于二者之间, 中和处理
    Z = (L - low) / (high - low);
else
    % 图像中90%像素值>high, 图像很亮, 无需处理
    Z = 1;
end

% 计算P
if im_std <= 3
    % 图像全局方差<3, 说明大部分地方颜色相同, 对比度很差, 此时P值取大
    P = 3;
elseif im_std <= 10
    % 3<图像全局方差<10, 说明原图有一定的对比度, 可以减少增强程度
    P = (27 - 2 * im_std) / 7;
else
    % 图像全局方差>10, 原图对比度较强, 可以不用增强
    P = 1;
end

% 构建映射查找表
Table = zeros(256, 256); % 行为调整后对应的亮度值, 列为原始图像对应的亮度值
for X = 1:256
    I = X / 255; % 公式2
    I = (power(I, 0.75*Z+0.25) + (1-I)*0.4*(1-Z) + power(I,2-Z)) / 2; % 公式3
    for Y = 1:256
        E = power(Y/X, P); % 公式8
        Table(Y, X) = round(255 * power(I, E)); % 公式7
    end
end

% 不同尺度高斯模糊融合: 尺度小时, 能提高局部对比度, 但可能整体看起来不是很协调; 尺度大时, 能获得整体图像更多信息, 但细节增强能力差
radius_s = 5;
radius_m = 20;
radius_l = 30;
% blur_s = imgaussfilt(luminance, 'FilterSize', radius_s*2+1);
% blur_m = imgaussfilt(luminance, 'FilterSize', radius_m*2+1);
% blur_l = imgaussfilt(luminance, 'FilterSize', radius_l*2+1);
blur_s = round(guide_filter(luminance, luminance, radius_s, 20)); % 使用guide滤波代替
blur_m = round(guide_filter(luminance, luminance, radius_m, 20));
blur_l = round(guide_filter(luminance, luminance, radius_l, 20));

luminance = double(luminance);
for Y = 1:h
    for X = 1:w
        lum = luminance(Y,X);
        if lum == 0
            enhanced(Y,X,:) = [0, 0, 0]';
        else
            val = (Table(blur_s(Y,X)+1, lum+1) + Table(blur_m(Y,X)+1, lum+1) + Table(blur_l(Y,X)+1, lum+1)) / 3; % 公式13
            cof = min(val / lum, 3); % 设置为4基本上能在增强效果和瑕疵之间达到一个平衡
            enhanced(Y,X,:) = im(Y,X,:) .* cof; % 公式14
        end
    end
end

enhanced = uint8(enhanced);

end