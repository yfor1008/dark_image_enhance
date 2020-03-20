% 3D变换，先进行 2D变换，用lambda2d阈值收缩，然后进行1D变换，
% lambda1d 阈值收缩。
function blk_tran3d=transform_3d(blk_3d,tran_mode,lambda2d,lambda1d)
    global blk_tran1d_s;
    global blk_2d_s;
    [m,n,blk_num]=size(blk_3d);
 
    %变换不同时，可能需要修改？？
    blk_2d_shrink=zeros(m,n,blk_num);
    blk_1d_shrink=zeros(m,n,blk_num);
 
    %(strcmp(tran_mode,'fft'))
    if(tran_mode==0)    %fft
        for i=1:1:blk_num
            blk_tran2d=fft2(blk_3d(:,:,i));
            blk_2d_shrink(:,:,i)=thres_shrink(blk_tran2d,lambda2d);
        end
        for i=1:1:m
            for j=1:1:n
                blk_tran1d=fft(blk_2d_shrink(i,j,:));
                blk_1d_shrink(i,j,:)=thres_shrink(blk_tran1d,lambda1d);
            end
        end
        blk_tran3d=blk_1d_shrink;
        %test 这是为了测试 还能否反变换回来。
        %blk_invtran1d=zeros(m,n,blk_num);
        %blk_est=zeros(m,n,blk_num);
        %for i=1:1:m
        %   for j=1:1:n
        %       blk_invtran1d(i,j,:)=ifft(blk_tran3d(i,j,:));
        %   end
        %end
        %for i=1:1:blk_num
        %   blk_est(:,:,i)=ifft2(blk_invtran1d(:,:,i));
        %end
         
    elseif(tran_mode==1)  %dct
        for i=1:1:blk_num
            blk_tran2d=dct2(blk_3d(:,:,i));
            blk_2d_shrink(:,:,i)=thres_shrink(blk_tran2d,lambda2d);
        end
        for i=1:1:m
            for j=1:1:n
                blk_tran1d=dct(blk_2d_shrink(i,j,:));
                blk_1d_shrink(i,j,:)=thres_shrink(blk_tran1d,lambda1d);
            end
        end
        blk_tran3d=blk_1d_shrink;
 
    elseif(tran_mode==2)    %dwt
        blk_2d_s={};
        blk_2d_shrink=[];%zeros()
        for i=1:1:blk_num
            [blk_tran2d_c,blk_tran2d_s]=wavedec2(blk_3d(:,:,i),2,'Bior1.5');
            blk_2d_shrink(:,i)=thres_shrink(blk_tran2d_c,lambda2d);
            blk_2d_s{i}=blk_tran2d_s;
        end
        %这里应该用 wavedec.因为是对1维？？
        [blk_tran1d_c,blk_tran1d_s]=wavedec2(blk_2d_shrink,1,'haar');
        blk_tran3d=thres_shrink(blk_tran1d_c,lambda1d);
%   elseif(strcmp(tran_mode,'db1')) %还未实现
%       blk_2d_s={};
%       blk_2d_shrink=[];%zeros()
%       for i=1:1:blk_num
%           [blk_tran2d_cA,blk_tran2d_cH,blk_tran2d_cV,blk_tran2d_cD]=...
%               dwt2(blk_3d(:,:,i),'db1');
%           blk_2d_shrink(:,i)=thres_shrink(blk_tran2d_c,lambda2d);
%           blk_2d_s{i}=blk_tran2d_s;
%       end
%       [blk_tran1d_c,blk_tran1d_s]=wavedec2(blk_2d_shrink,1,'haar');
%       blk_tran3d=thres_shrink(blk_tran1d_c,lambda1d);
    else
        error('tran_mode error');
    end