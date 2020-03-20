function T_out = Illumination_filter(T_in, method, ksize, sigma1, sigma2, alpha, sharpness, maxIter)
% Inputs:
%        T_in:
%        ksize: 窗口的大小为ksize x ksize
%        method:
%        sigma1:
%        sigma2:
%        alpha:
%        sharpness:
%        maxIter:
% Outputs:
%         T_out:
%
% Author: HSW
% Date: 2018-04-27
%
if ~exist('method','var')
    method = 'Max';
end

if strcmp(method, 'Max') == 1
    if ~exist('ksize','var')
        ksize = 5;
    end
elseif strcmp(method, 'Mean') == 1
    if ~exist('ksize','var')
        ksize = 5;
    end
elseif strcmp(method, 'BF') == 1
    if ~exist('ksize','var')
        ksize = 5;
    end
    
    if ~exist('sigma1','var')
        sigma1 = 5;
    end
    
    if ~exist('sigma2','var')
        sigma2 = 7;
    end
elseif strcmp(method, 'TV') == 1
    if ~exist('alpha','var')
        alpha = 0.15;
    end
    
    if ~exist('maxIter','var')
        maxIter = 50;
    end
elseif strcmp(method, 'SA') == 1
    if ~exist('alpha','var')
        alpha = 0.15;
    end
    
    if ~exist('sigma1','var')
        sigma1 = 5;
    end
    
    if ~exist('sharpness', 'var')
        sharpness = 0.02;
    end
    
    if ~exist('maxIter','var')
        maxIter = 4;
    end
end


if strcmp(method, 'Max') == 1
    % 局部极大值
    T_out = filter_Max(T_in, ksize);
elseif strcmp(method, 'Mean') == 1
    T_out = filter_Mean(T_in, ksize);
elseif strcmp(method, 'BF') == 1
    % 双边滤波
    T_out = filter_BF(T_in, ksize, sigma1, sigma2);
elseif strcmp(method, 'TV') == 1
    T_out = filter_TV(T_in, alpha, maxIter);
elseif strcmp(method, 'SA') == 1
    T_out = filter_SA(T_in, alpha, sigma1, sharpness, maxIter);
end
end

function T_out = filter_Max(T_in, ksize)
% Inputs:
%        T_in:
%        ksize:
% Outputs:
%        T_out:
% Author: HSW
% Date: 2018-04-27
%

[m, n] = size(T_in);
hsize = floor(ksize / 2);
T_out = T_in;
for i = 1:m
    for j = 1:n
        patch = T_in(max(1, i - hsize):min(m, i + hsize), max(1, j - hsize):min(n, j + hsize));
        T_out(i,j) = max(max(patch));
    end
end
disp('run filter_max');
end

function T_out = filter_Mean(T_in, ksize)
% Inputs:
%        T_in:
%        ksize:
% Outputs:
%        T_out:
% Author: HSW
% Date: 2018-04-27
%

[m, n] = size(T_in);
hsize = floor(ksize / 2);
T_out = T_in;
for i = 1:m
    for j = 1:n
        patch = T_in(max(1, i - hsize):min(m, i + hsize), max(1, j - hsize):min(n, j + hsize));
        T_out(i,j) = mean(mean(patch));
    end
end
disp('run filter_mean');
end

function T_out = filter_BF(T_in, ksize, sigma1, sigma2)
% Inputs:
%        T_in:
%        ksize:
%        sigma1:
%        sigma2:
% Outputs:
%        T_out:
% Author: HSW
% Date: 2018-04-27
%

[m, n] = size(T_in);
hsize = floor(ksize / 2);
T_out = T_in;

xx = -hsize : hsize;
yy = -hsize : hsize;
kernel1 = zeros(ksize, ksize);
for ii = 1:ksize
    for jj = 1:ksize
        kernel1(ii, jj) = exp(-(xx(ii) * yy(jj) / sigma1)^2);
    end
end
kernel1 = kernel1 ./ sum(sum(kernel1));  % 归一化

for i = hsize + 1:m - hsize
    for j = hsize + 1:n - hsize
        kernel2 = exp(-((T_in(i - hsize : i + hsize, j - hsize:j + hsize) - repmat(T_in(i,j), ksize, ksize)) ./ sigma2).^2);
        kernel2 = kernel2 ./ sum(sum(kernel2));
        T_out(i,j) = sum(sum(kernel1 .* kernel2 .* T_in(i - hsize : i + hsize, j - hsize : j + hsize) ./ sum(sum(kernel1 .* kernel2))));
    end
end
end

function T_out = filter_TV(T_in, alpha, maxIter)
% Inputs:
%        T_in:
%        alpha:
%        maxIter:
% Outputs:
%        T_out:
% Author: HSW
% Date: 2018-04-27
%
epsilon = 0.00001;
dt = 0.1; 
J = T_in;
for iter = 1:maxIter
    
    DfJx=J([2:end end],:)-J;     %函数关于X的一阶偏导（向后差分）
    DbJx=J-J([1 1:end-1],:);     %函数关于X的一阶偏导（向前差分）
    DfJy=J(:,[2:end end])-J;     %函数关于Y的一阶偏导（向后差分）
    DbJy=J-J(:,[1 1:end-1]);     %函数关于Y的一阶偏导（向前差分）
    
    TempDJx=(epsilon+DfJx.*DfJx+((sign(DfJy)+sign(DbJy)).*min(abs(DfJy),abs(DbJy))./2).^2).^(1/2);%求梯度的模    
    TempDJy=(epsilon+DfJy.*DfJy+((sign(DfJx)+sign(DbJx)).*min(abs(DfJx),abs(DbJx))./2).^2).^(1/2);
      
    DivJx=DfJx./(TempDJx + (TempDJx == 0) * epsilon);    
    DivJy=DfJy./(TempDJy + (TempDJy == 0) * epsilon);

    %求散度    
    Div=DivJx-DivJx([1 1:end-1],:)+DivJy-DivJy(:,[1 1:end-1]);  
    J= J + dt * alpha* Div - dt * (J-T_in);                  %产生迭代
end
T_out = J;
end

function T_out = filter_SA(T_in, alpha, sigma, sharpness, maxIter)
% Inputs:
%        T_in:
%        alpha:
%        sigma:
%        maxIter:
%        sharpness:
% Outputs:
%        T_out:
% Author: HSW
% Date: 2018-04-27
%
% 论文描述的解法没有完全理解，就按照论文中的文献[17]的迭代法进行实现，只是不迭代更新Wh 和 Wv
T_out = tsmooth(T_in, alpha, sigma, sharpness, maxIter);
end


function S = tsmooth(I,lambda,sigma,sharpness,maxIter)
% 参考代码： https://blog.csdn.net/songzitea/article/details/12851723#
if (~exist('lambda','var'))
    lambda=0.01;
end
if (~exist('sigma','var'))
    sigma=3.0;
end
if (~exist('sharpness','var'))
    sharpness = 0.02;
end
if (~exist('maxIter','var'))
    maxIter=4;
end
I = im2double(I);
x = I;
sigma_iter = sigma;
lambda = lambda/2.0;
dec=2.0;
[wx, wy] = computeTextureWeights(x, sigma_iter, sharpness); % 与文献[17]不同，文献[17]每次迭代都改变wx, wy
for iter = 1:maxIter
    % [wx, wy] = computeTextureWeights(x, sigma_iter, sharpness);
    x = solveLinearEquation(I, wx, wy, lambda);
    sigma_iter = sigma_iter/dec;
    if sigma_iter < 0.5
        sigma_iter = 0.5;
    end
end
S = x;
end

function [retx, rety] = computeTextureWeights(fin, sigma,sharpness)

fx = diff(fin,1,2);
fx = padarray(fx, [0 1 0], 'post');
fy = diff(fin,1,1);
fy = padarray(fy, [1 0 0], 'post');

vareps_s = sharpness;
vareps = 0.001;

wto = max(sum(sqrt(fx.^2+fy.^2),3)/size(fin,3),vareps_s).^(-1);
fbin = lpfilter(fin, sigma);
gfx = diff(fbin,1,2);
gfx = padarray(gfx, [0 1], 'post');
gfy = diff(fbin,1,1);
gfy = padarray(gfy, [1 0], 'post');
wtbx = max(sum(abs(gfx),3)/size(fin,3),vareps).^(-1);
wtby = max(sum(abs(gfy),3)/size(fin,3),vareps).^(-1);
retx = wtbx.*wto;
rety = wtby.*wto;

retx(:,end) = 0;
rety(end,:) = 0;

end

function ret = conv2_sep(im, sigma)
ksize = bitor(round(5*sigma),1);
g = fspecial('gaussian', [1,ksize], sigma);
ret = conv2(im,g,'same');
ret = conv2(ret,g','same');
end

function FBImg = lpfilter(FImg, sigma)
FBImg = FImg;
for ic = 1:size(FBImg,3)
    FBImg(:,:,ic) = conv2_sep(FImg(:,:,ic), sigma);
end
end

function OUT = solveLinearEquation(IN, wx, wy, lambda)
[r,c,ch] = size(IN);
k = r*c;
dx = -lambda*wx(:);
dy = -lambda*wy(:);
B(:,1) = dx;
B(:,2) = dy;
d = [-r,-1];
A = spdiags(B,d,k,k);
e = dx;
w = padarray(dx, r, 'pre'); w = w(1:end-r);
s = dy;
n = padarray(dy, 1, 'pre'); n = n(1:end-1);
D = 1-(e+w+s+n);
A = A + A' + spdiags(D, 0, k, k);
if exist('ichol','builtin')
    L = ichol(A,struct('michol','on'));
    OUT = IN;
    for ii=1:ch
        tin = IN(:,:,ii);
        [tout, flag] = pcg(A, tin(:),0.1,100, L, L');
        OUT(:,:,ii) = reshape(tout, r, c);
    end
else
    OUT = IN;
    for ii=1:ch
        tin = IN(:,:,ii);
        tout = A\tin(:);
        OUT(:,:,ii) = reshape(tout, r, c);
    end
end

end
