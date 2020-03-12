function [enhanced] = nonlinear_color_enhance(im)
% nonlinear_color_enhance - ����Ӧ������ͼ����ǿ
%
% input:
%   - im: h*w*3, rgbͼ��
% output:
%   - enhanced: h*w*3, ��ǿ��ͼ��
%
% docs:
%   - ��Adaptive and integrated neighborhood-dependent approach for nonlinear enhancement of color images��
%   - ����Ӧ������ǿ: ��߰����ص�����, ��ǿ������ȴѹ����ͼ��Ķ�̬��Χ, �����˶Աȶ�
%   - ����Ӧ�Աȶ���ǿ: ��������ֵ�����������ƽ��ֵ�ı���, ����>1ʱ, ��������ֵ, ����<1ʱ, ��С����ֵ
%   - https://www.cnblogs.com/Imageshop/p/11665100.html
%

luminance = rgb2gray(im); % ��ʽ1
[h, w] = size(luminance);
pix_num = h * w;

im = double(im);
enhanced = zeros(h, w, 3);

% ����ֱ��ͼ
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

% ����Z
low = 50;
high = 150;
if L <= low
    % ͼ����%10������ֵ<low, ͼ��ƫ��, ����ǿ����  
    Z = 0;
elseif L <= high
    % ���ڶ���֮��, �кʹ���
    Z = (L - low) / (high - low);
else
    % ͼ����90%����ֵ>high, ͼ�����, ���账��
    Z = 1;
end

% ����P
if im_std <= 3
    % ͼ��ȫ�ַ���<3, ˵���󲿷ֵط���ɫ��ͬ, �ԱȶȺܲ�, ��ʱPֵȡ��
    P = 3;
elseif im_std <= 10
    % 3<ͼ��ȫ�ַ���<10, ˵��ԭͼ��һ���ĶԱȶ�, ���Լ�����ǿ�̶�
    P = (27 - 2 * im_std) / 7;
else
    % ͼ��ȫ�ַ���>10, ԭͼ�ԱȶȽ�ǿ, ���Բ�����ǿ
    P = 1;
end

% ����ӳ����ұ�
Table = zeros(256, 256); % ��Ϊ�������Ӧ������ֵ, ��Ϊԭʼͼ���Ӧ������ֵ
for X = 1:256
    I = X / 255; % ��ʽ2
    I = (power(I, 0.75*Z+0.25) + (1-I)*0.4*(1-Z) + power(I,2-Z)) / 2; % ��ʽ3
    for Y = 1:256
        E = power(Y/X, P); % ��ʽ8
        Table(Y, X) = round(255 * power(I, E)); % ��ʽ7
    end
end

% ��ͬ�߶ȸ�˹ģ���ں�: �߶�Сʱ, ����߾ֲ��Աȶ�, ���������忴�������Ǻ�Э��; �߶ȴ�ʱ, �ܻ������ͼ�������Ϣ, ��ϸ����ǿ������
radius_s = 5;
radius_m = 20;
radius_l = 30;
% blur_s = imgaussfilt(luminance, 'FilterSize', radius_s*2+1);
% blur_m = imgaussfilt(luminance, 'FilterSize', radius_m*2+1);
% blur_l = imgaussfilt(luminance, 'FilterSize', radius_l*2+1);
blur_s = round(guide_filter(luminance, luminance, radius_s, 20)); % ʹ��guide�˲�����
blur_m = round(guide_filter(luminance, luminance, radius_m, 20));
blur_l = round(guide_filter(luminance, luminance, radius_l, 20));

luminance = double(luminance);
for Y = 1:h
    for X = 1:w
        lum = luminance(Y,X);
        if lum == 0
            enhanced(Y,X,:) = [0, 0, 0]';
        else
            val = (Table(blur_s(Y,X)+1, lum+1) + Table(blur_m(Y,X)+1, lum+1) + Table(blur_l(Y,X)+1, lum+1)) / 3; % ��ʽ13
            cof = min(val / lum, 3); % ����Ϊ4������������ǿЧ����覴�֮��ﵽһ��ƽ��
            enhanced(Y,X,:) = im(Y,X,:) .* cof; % ��ʽ14
        end
    end
end

enhanced = uint8(enhanced);

end