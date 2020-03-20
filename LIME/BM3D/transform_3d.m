% 3D�任���Ƚ��� 2D�任����lambda2d��ֵ������Ȼ�����1D�任��
% lambda1d ��ֵ������
function blk_tran3d=transform_3d(blk_3d,tran_mode,lambda2d,lambda1d)
    global blk_tran1d_s;
    global blk_2d_s;
    [m,n,blk_num]=size(blk_3d);
 
    %�任��ͬʱ��������Ҫ�޸ģ���
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
        %test ����Ϊ�˲��� ���ܷ񷴱任������
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
        %����Ӧ���� wavedec.��Ϊ�Ƕ�1ά����
        [blk_tran1d_c,blk_tran1d_s]=wavedec2(blk_2d_shrink,1,'haar');
        blk_tran3d=thres_shrink(blk_tran1d_c,lambda1d);
%   elseif(strcmp(tran_mode,'db1')) %��δʵ��
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