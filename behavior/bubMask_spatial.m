function [posit,mask]=bubMask_spatial(qteBulles,spaceSize,masque,bubStd)
% Generate the bubbles mask
bulle=bubble(bubStd);
prob_tmp = qteBulles/sum(masque(:));
tmp=rand(spaceSize^2,1) .* masque(:);
tmp=reshape(tmp>=(1-prob_tmp),spaceSize,spaceSize); % makes the criteria probabilistic
posit=tmp;
masque2D=filter2(bulle,tmp);
masque2D = min(max(masque2D, 0), 1); % this is better

mask= repmat(masque2D,[1 1 3]);
end