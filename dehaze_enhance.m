function d = dehaze_enhance(im, w)
% dehaze_enhance - 利用去雾算法增强
%
% input:
%   - im: h*w*3, rgb图像
% output:
%   - outval: h*w*3, 增强后图像
%
% docs:
%   - 《Fast efficient algorithm for enhancement of low lighting video》
%   - https://blog.csdn.net/qq_21089969/article/details/74890200
%   - https://www.jianshu.com/p/6e8ea77adf8b
%   - 取反的低亮度图像与有雾图像类似:
%   - 取反的低亮度图像中, 天空和远景部分的像素值非常大
%   - 取反的低亮度图像中, 非天空区域中, RGB通道至少有一个像素值非常低
%   - 对有天空的暗图像及有雾的暗图像效果不错!
%
if max(im(:)) > 1
    im = im2double(im);
end
isize = size(im);
im =  padarray(im,[2*5,2*5,0],'replicate');
if (~exist('w','var') || w<=0)
    w = 0.8;
end

R = 1 - im; % 取反

Ir = R(:,:,1);
Ig = R(:,:,2);
Ib = R(:,:,3);

Ir1 = Ir(:);
Ig1 = Ig(:);
Ib1 = Ib(:);
Il = (Ir1+Ig1+Ib1)./3;
n = length(Il);
N = floor(n*0.002);

Ir_d = ordfilt2(Ir,1,ones(7,7)); % 最小值滤波
Ig_d = ordfilt2(Ig,1,ones(7,7));
Ib_d = ordfilt2(Ib,1,ones(7,7));
darkc = min(min(Ir_d(:),Ig_d(:)),Ib_d(:));
[~, i] = sort(darkc,1,'descend');
temp = Il(i(1:N));
[~,j] = sort(temp,1,'descend');
Ar = Ir1(i(j(1)));
Ag = Ig1(i(j(1)));
Ab = Ib1(i(j(1)));

t = max(1-w.*min(min(Ir./Ar,Ig./Ag),Ib./Ab),10^(-7));
lc = t<0.5; % 对图像有较大影像
t(lc) = 2 * t(lc).^2;

Sr = (Ir - Ar)./t + Ar;
Sg = (Ig - Ag)./t + Ag;
Sb = (Ib - Ab)./t + Ab;

Rd = cat(3,Sr, Sg, Sb);
d = 1 - Rd(11:isize(1)+10, 11:isize(2)+10, :);

end
