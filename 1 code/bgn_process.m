function [ref] = bgn_process(A,up_limit,c)
% this function aims to calcute the background

temp = 0;
n = 0;
for k = 1:c
    if max(A(:,k)) < up_limit
        temp = temp + A(:,k);
        n = n+1;
    end
    ref = temp/n;
    ref = repmat(ref,1,c);
end
end
