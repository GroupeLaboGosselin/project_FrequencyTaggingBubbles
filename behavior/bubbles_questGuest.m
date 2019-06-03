% CHECK THIS AND CORRECT WITH bubMASK
function [BestGuess,minBulle,surfaceBulle]= bubbles_questGuest(qteBulles,spaceSize,bubStd)




% AJUSTEMENTS QUEST
% - On calcule la valeur a ajuster (variant de 0 a 1),
%  la qte de pixel revelee par les bulles.
bulle=bubble(bubStd);
tmp=rand(spaceSize ^2,1);
[y,~]=sort(tmp);
tmp=reshape(tmp<=y(1),spaceSize,spaceSize);
tmp=filter2(bulle,tmp);
tmp = (tmp-min(tmp(:)))/(max(tmp(:))-min(tmp(:)));
surfaceBulle=sum(tmp(:)); % number of pixels in a single bubble
minBulle=(surfaceBulle/spaceSize.^2);
BestGuess=((qteBulles*surfaceBulle)/spaceSize.^2)-minBulle;


end