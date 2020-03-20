function img_out = Lime_enhance(img_in, Model, method, denoiseFlag, displayFlag)
% Inputs:
%        Model: 'Normal' 论文公式(1)模型, 'Dehaze'论文公式(2)模型
%        method: 'Max' / 'Mean' / 'BF' / 'TV' / 'SA'
% Outputs:
%        img_out:
%
epsilon = 0.0001;
[m, n, dims] = size(img_in);

if strcmp(Model, 'Normal') == 1
    % 估计光照亮度
    T = Illumination(img_in, 'max_c'); % 亮通道
    
    if displayFlag
        figure;
        imshow(T,[]);
        title('T');
    end
    
    % 亮通道光滑
    smoothT = Illumination_filter(T, method);
    
    if displayFlag
        figure;
        imshow(smoothT, []);
        title('smoothT');
    end
    
    % gama变换
    gama = 0.8;
    gamaT = smoothT.^gama;
    
    if displayFlag
        figure;
        imshow(gamaT, []);
        title('gamaT');
    end
    
    % I = img_in ./ repmat((gamaT + (gamaT == 0) * epsilon), [1, 1, 3]);
    I = 1 - ((ones([m, n, dims]) - img_in) - repmat(0.95 * (1 - gamaT), [1, 1, 3])) ./ repmat((gamaT + (gamaT == 0) * epsilon), [1, 1, 3]);
    I(I < 0) = 0;
    I(I > 1) = 1;
    
    if displayFlag
        figure;
        imshow(I, []);
        title('未除噪增强图像');
    end
    
    % 进行噪声处理
    if denoiseFlag == 1
        I_noise = uint8(I .* 255);
        I_YCbCr = rgb2ycbcr(I_noise);
        Y = I_YCbCr(:, :, 1);
        maxY = max(max(Y));
        minY = min(min(Y));
        % BM3D去噪
        Y_denoise = BM3D(Y, 0, 5, 1, 1);
        maxY_denoise = max(max(Y_denoise));
        minY_denoise = min(min(Y_denoise));
        Y_denoise = (double(Y_denoise) - double(minY_denoise)) ./ double(maxY_denoise - minY_denoise) .* double(maxY - minY) + double(minY);
        I_YCbCr(:, :, 1) = Y_denoise;
        I_denoise = ycbcr2rgb(I_YCbCr);
        I_denoise = double(I_denoise) ./ 255.;
        % 图像融合
        I = double(I) .* double(repmat(gamaT, [1, 1, 3])) + double(I_denoise) .* double(repmat((1 - gamaT), [1, 1, 3]));
    end
    
elseif strcmp(Model, 'Dehaze') == 1
    img_in = 1 - img_in;
    ksize = 5;
    [img_in, I, J, T_est, T, A] = removeHaze(img_in,ksize);
    I = 1 - I;
end

img_out = I;
end