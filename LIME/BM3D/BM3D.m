function img_out = BM3D(img_in, tran_mode, sigma, displayMidResult, displayProcess)
% 参考文献: 《Image denoising by sparse 3D transform-domain collaborative filtering》
%           《An Analysis and Implementation of the BM3D Image Denoising Method》
% Inputs:
%        img_in: 噪声图像，必须为矩形方阵
%        tran_mode: = 0, fft; = 1, dct; = 2, dwt, = 3, db1
%        displayMidResult: 0/1
%        displayProcess: 0/1
% Outputs: 
%        img_out: 去噪图像
% 代码下载链接： http://www.pudn.com/Download/item/id/2703664.html
%
% 
if ~exist('tran_mode', 'var') 
    tran_mode = 0; 
end 

if ~exist('sigma', 'var') 
    sigma = 10; 
end 

if ~exist('displyaMidResult', 'var') 
    displayMidResult = 0; 
end 

if ~exist('displayProcess', 'var') 
    displayProcess = 0; 
end 

img_noise = img_in;
[row,col] = size(img_noise); 
tic
% first step
kHard=8;          % 块大小
pHard=4;          % 块移动间隔
lambda_distHard=0;% 求相似的距离时，变换后，收缩的阈值
nHard=40;         % 搜索窗口大小
NHard=28;         % 最多相似块个数
tauHard=5000;     % 最大的相似距离for fft
beta=2;
% tauHard=50000;% for dct
if(tran_mode==0)        %fft
    lambda2d=400;
    lambda1d=500;
    lambda2d_wie=50;
    lambda1d_wie=500;
elseif(tran_mode == 1)  %dct
    lambda2d=50;
    lambda1d=80;
    lambda2d_wie=20;
    lambda1d_wie=60;
elseif(tran_mode == 2)  %dwt
    lambda2d=50;
    lambda1d=80;
    lambda2d_wie=20;
    lambda1d_wie=60;
end

%kaiser 窗，实际代码中可能没有用到。
kaiser_win=kaiser(kHard,1)*kaiser(kHard,1)';
%图像分块，并且做变换，为找相似块做准备
[block,tran_block,block2row_idx,block2col_idx]=im2block(img_noise,kHard,pHard,lambda_distHard,0);
% block number row:行方向上的 块的个数
bn_r=floor((row-kHard)/pHard)+1;
bn_c=floor((col-kHard)/pHard)+1;
%基础估计的图像
img_basic_sum=zeros(row,col);
img_basic_weight=zeros(row,col);

%显示处理进度
is_disp_process    = displayProcess;
process_step_total = bn_r/10;
process_step_cnt   = 0;
%基础估计
fprintf('BM3D: First Stage Start...\n');
%对每个块遍历
for i=1:1:bn_r
    if((is_disp_process) &&(i>process_step_total))
        process_step_total=process_step_total+bn_r/10;
        process_step_cnt=process_step_cnt+1;
        fprintf('  process:%d/10\n',process_step_cnt)
    end
    for j=1:1:bn_c
        [sim_blk,sim_num,sim_blk_idx]=search_similar_block(i,j,block,tran_block,kHard,floor(nHard/pHard),...
            bn_r,bn_c,tauHard,NHard);
        %toc
        %         test fine_similar block
        num_sim=size(sim_blk_idx,3);
        %tic
        %做3D变换，并且用阈值收缩
        tran3d_blk_shrink=transform_3d(sim_blk,tran_mode,lambda2d,lambda1d);
        %toc
        NHard_P=nnz(tran3d_blk_shrink);%non_zero_num
        if(NHard_P >1)
            wHard_P=1/NHard_P;
        else
            wHard_P=1;
        end
        Wwin2D= kaiser(kHard, beta) * kaiser(kHard, beta)'; % Kaiser window used in the hard-thresholding part
        wWien_P=Wwin2D*wHard_P;
        %wHard_P=wHard_P*kaiser_win;% 要不要？
        %tic
        % 3D逆变换
        blk_est=inv_transform_3d(tran3d_blk_shrink,tran_mode);
        %pause
        %可能是复数，虚部为0或接近于0，所以只取实部
        %max(abs(imag(blk_est(:))))
        blk_est=real(blk_est);
        %toc
        %tic
        for k=1:sim_num
            idx=sim_blk_idx(k);
            ir=block2row_idx(idx);
            jr=block2col_idx(idx);
            %实在算不清了。利用提前保存好 块索引和左上角坐标的 对应关系
            %ir=floor((idx-1)*pHard/col)+1;
            %jr=(idx-1)*pHard-(i-1)*col+1;
            img_basic_sum(ir:ir+kHard-1,jr:jr+kHard-1)=...
                img_basic_sum(ir:ir+kHard-1,jr:jr+kHard-1)+wHard_P*blk_est(:,:,k);
            img_basic_weight(ir:ir+kHard-1,jr:jr+kHard-1)=...
                img_basic_weight(ir:ir+kHard-1,jr:jr+kHard-1)+wHard_P;
        end
        %toc
        %pause
    end
end
fprintf('BM3D: First Stage End...\n');
img_basic=img_basic_sum./img_basic_weight;

if displayMidResult
    figure;
    imshow(img_basic,[]);
    title('BM3D:Fist Stage Result'); 
%     psnr=20*log10(255/sqrt(mean((img_basic(:)-img(:)).^2)))
end

% second step
kWien=kHard;
pWien=pHard;
lambda_distWien=lambda_distHard;
nWien=nHard;%搜索窗口大小
NWien=NHard;%最多相似块个数
tauWien=tauHard;
sigma2=sigma*sigma;

[block_basic,tran_block_basic,block2row_idx_basic,block2col_idx_basic]=im2block(img_basic,kWien,pWien,lambda_distWien,0);
bn_r=floor((row-kWien)/pWien)+1;
bn_c=floor((col-kWien)/pWien)+1;
img_wien_sum=zeros(row,col);
img_wien_weight=zeros(row,col);

process_step_total=bn_r/10;
process_step_cnt=0;
fprintf('BM3D: Second Stage Start...\n');
for i=1:1:bn_r
    if((is_disp_process) &&(i>process_step_total))
        process_step_total=process_step_total+bn_r/10;
        process_step_cnt=process_step_cnt+1;
        fprintf('  process:%d/10\n',process_step_cnt)
    end
    for j=1:1:bn_c
        [sim_blk_basic,sim_num,sim_blk_basic_idx]=search_similar_block(i,j,block_basic,tran_block_basic,kWien,floor(nWien/pWien),...
            bn_r,bn_c,tauWien,NWien);
        %对基础块进行3D变换，求得omega_P.
        tran3d_blk_basic=transform_3d(sim_blk_basic,tran_mode,lambda2d_wie,lambda1d_wie);
        omega_P=(tran3d_blk_basic.^2)./((tran3d_blk_basic.^2)+sigma2);
        %用 基础块得到的相似块的索引，来找到 噪声块的相似块，并进行3D变换
        tran3d_blk=transform_3d(block(:,:,sim_blk_basic_idx),tran_mode,lambda2d_wie,lambda1d_wie);
        blk_est=inv_transform_3d(omega_P.*tran3d_blk,tran_mode);
        %可能是复数，虚部为0或接近于0，所以只取实部
        %max(abs(imag(blk_est(:))))
        blk_est=real(blk_est);
        NWien_P=nnz(omega_P); %IPOL文中8式中矩阵范数？应该如何求？？
        if(NWien_P >1)
            wWien_P=1/(NWien_P);
        else
            wWien_P=1;
        end
        
        %         wWien_P=wWien_P/sigma2;
        for k=1:sim_num
            idx=sim_blk_basic_idx(k);
            ir=block2row_idx_basic(idx);
            jr=block2col_idx_basic(idx);
            img_wien_sum(ir:ir+kWien-1,jr:jr+kWien-1)=...
                img_wien_sum(ir:ir+kWien-1,jr:jr+kWien-1)+wWien_P*blk_est(:,:,k);
            img_wien_weight(ir:ir+kWien-1,jr:jr+kWien-1)=...
                img_wien_weight(ir:ir+kWien-1,jr:jr+kWien-1)+wWien_P;
        end
    end
end
fprintf('BM3D: Second Stage End\n');
img_wien=img_wien_sum./img_wien_weight;

if displayMidResult
    figure;
    imshow(img_wien,[]);
    title('BM3D: 去噪结果'); 
%     psnr=20*log10(255/sqrt(mean((img_wien(:)-img(:)).^2)))
end
img_out = img_wien;
toc