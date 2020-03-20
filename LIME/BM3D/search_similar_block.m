% 搜索相似块
% ik,jk: 块的左上角 坐标。
%bn_r: block_num_row,分块之后,沿着行的方向上有多少个块。
%用 tran_block 来算距离，得到的sim_blk 却是来自block.
%np: 搜索窗口的大小/块移动步距。
% tau: 最大相似距离
% max_sim_num: 最多相似块 个数
function [sim_blk,sim_num,sim_blk_idx]=search_similar_block(ik,jk,block,tran_block,k,np,bn_r,bn_c,tau,max_sim_num)
    % 搜索窗口的 左上角 和 右下角 的块的坐标。s,e:start,end.
    in_s=max(ik-floor(np/2),1);
    jn_s=max(jk-floor(np/2),1);
    in_e=min(ik+floor(np/2),bn_r);
    jn_e=min(jk+floor(np/2),bn_c);
    % 当前参考块
    ref_blk=tran_block(:,:,((ik-1)*bn_c+jk));
    k2=k*k;
    cnt=0;
     
    %dist=[];
    %blk_idx=[];
    %如果参考块在图像的边缘，那它周围的搜索区域就不会是完整的n*n
    %所以不能单纯的由cnt 反推出idx.
    %for ii=in_s:1:in_e
    %   for jj=jn_s:1:jn_e
    %       cnt=cnt+1;
    %       %idx=ii*bn_c+jj;
    %       idx=(ii-1)*bn_c+(jj-1)+1;
    %       cur_blk=tran_block(:,:,idx);
    %       blk_idx(cnt)=idx;
    %       dist(cnt)=norm(cur_blk-ref_blk);% 只是比较大小，就没必要归一化了 /k2;
    %   end
    %end
    %下面找相似块的方法要比上面的快一些
 
    ii=in_s:1:in_e;
    jj=jn_s:1:jn_e;
    [II,JJ]=meshgrid(ii,jj);
    IDX=(II-1)*bn_c+JJ;
    blk_idx=IDX(:);
    cur_blk=tran_block(:,:,blk_idx);
    cnt=size(cur_blk,3);
    ref_blk_mat=repmat(ref_blk,[1,1,cnt]);
    % shit! 要不要用norm ? 奇异值能衡量矩阵的相似程度？
    % norm(cur_blk-ref_blk_mat);
    delta_blk=cur_blk-ref_blk_mat;
    dist=sum(sum(delta_blk.*delta_blk,1),2);
 
    % 排序，找到最相似的块
    [dist_sort,dist_idx]=sort(dist);
    max_num=min(cnt,max_sim_num);%可能当前块中个数还没有max_sim_num多
    if(dist_sort(max_num)<tau)
        sim_num=max_num;
    else
        sim_num=sum(dist_sort(1:max_num)<tau);
    end
    cnt_idx=dist_idx(1:sim_num);
    sim_blk_idx=blk_idx(cnt_idx);
    sim_blk=block(:,:,sim_blk_idx);