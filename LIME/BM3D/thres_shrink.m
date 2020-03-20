 % 阈值收缩，abs小于thres的 赋值 0
function [val]=thres_shrink(data,thres)
    val=data;
    idx=find(abs(data)<thres);
    val(idx)=0;
end