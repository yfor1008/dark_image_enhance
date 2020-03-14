function [rgb] = multi_fusion_1(im)
% multi_fusion - 邻域非线性增强
%
% input:
%   - im: h*w*3, rgb图像
% output:
%   - enhanced: h*w*3, 增强后
%
% docs:
%   - 《A Fusion-based Enhancing Method for Weakly Illuminated Images》
%   - https://xueyangfu.github.io/projects/sp2016.html
%   - https://blog.csdn.net/weixin_44690935/article/details/104335351
%

hsv = rgb2hsv(im);
H = hsv(:,:,1);
S = hsv(:,:,2);
V = hsv(:,:,3);
[h,w] = size(V);

im1 = V .^ 0.35;

a = 50 / 255;
index = V < a;
M = sum(index(:));
N = h*w;
gamma = (N-M) / N;
im2 = 1 - (1 - V) .^ gamma;

im3 = gum(V);

sigma = 0.25;
aver = 0.3;
WE1 = exp(-(im1 - aver).^2 / (2*sigma^2));
WE2 = exp(-(im2 - aver).^2 / (2*sigma^2));
WE3 = exp(-(im3 - aver).^2 / (2*sigma^2));
sumW=WE1+WE2+WE3;
W1=WE1./sumW;
W2=WE2./sumW;
W3=WE3./sumW;

% Multi-scale fusion
level=4;
G1=gaussian_pyramid(W1,level);
G2=gaussian_pyramid(W2,level);
G3=gaussian_pyramid(W3,level);
L1=laplacian_pyramid(im1,level);
L2=laplacian_pyramid(im2,level);
L3=laplacian_pyramid(im3,level);
for i=1:level
    F{i}=G1{i}.*L1{i}+G2{i}.*L2{i}+G3{i}.*L3{i};
end
Vfinal=pyramid_reconstruct(F);
% imshow(Vfinal);figure;
hsv=cat(3,H,S,Vfinal);
rgb=hsv2rgb(hsv);

end

function [v2] = gum(x1)
x1=im2double(x1);
[ m, n, k ] = size( x1 );
x = x1;
mask = ( 1 / 25 ) * ones( 5, 5 );
y1 = conv2( x, double( mask ), 'same' );
d = x - y1;
g = 3 .* d;
v1 = y1 + g;
y2 = medfilt2( x, [ 3, 3 ] );
X = ( 1 - x ) ./ max( x, 0.01 );
Y = ( 1 - y2 ) ./ max( y2, 0.01 );
I = ones( m, n );
d1 = I ./ ( 1 + ( X ./ Y ) );
h = adapthisteq( y2, 'clipLimit', 0.005 );
c = ( 2 .* d1 ) - 1;
Gmax = 5;Gmin = 1;eta = 0.5;
beta = ( Gmax - Gmin ) / ( 1 - exp(  - 1 ) );
alpha = ( Gmax - beta );
gama = alpha + ( beta * exp(  - 1 .* abs( c ) .^ eta ) );
D = ( 1 - d1 ) ./ max( d1, 0.01 );
g = I ./ ( 1 + D .^ gama );
G = ( 1 - g ) ./ max( g, 0.01 );
H = ( 1 - h ) ./ max( h, 0.01 );
v2 = I ./ max( 0.1, ( 1 + ( H .* G ) ) );
t = v2 > 1;
v2( t ) = x( t );

end

function out = gaussian_pyramid(img, level)
h = 1/16* [1, 4, 6, 4, 1];
filt = h'*h;
out{1} = imfilter(img, filt, 'replicate', 'conv');
temp_img = img;
for i = 2 : level
    temp_img = temp_img(1 : 2 : end, 1 : 2 : end);
    out{i} = imfilter(temp_img, filt, 'replicate', 'conv');
end
end

function out = laplacian_pyramid(img, level)
h = 1/16* [1, 4, 6, 4, 1];
%filt = h'*h;
out{1} = img;
temp_img = img;
for i = 2 : level
    temp_img = temp_img(1 : 2 : end, 1 : 2 : end);
    %out{i} = imfilter(temp_img, filt, 'replicate', 'conv');
    out{i} = temp_img;
end
% calculate the DoG
for i = 1 : level - 1
    [m, n] = size(out{i});
    out{i} = out{i} - imresize(out{i+1}, [m, n]);
end
end

function out = pyramid_reconstruct(pyramid)
level = length(pyramid);
for i = level : -1 : 2
    %temp_pyramid = pyramid{i};
    [m, n] = size(pyramid{i - 1});
    %out = pyramid{i - 1} + imresize(temp_pyramid, [m, n]);
    pyramid{i - 1} = pyramid{i - 1} + imresize(pyramid{i}, [m, n]);
end
out = pyramid{1};
end