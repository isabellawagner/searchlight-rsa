
function [pattern] = searchlight_getdata(fmri,curr_roi,indices,permute,Pseq,RSAtrue)

% ----------------------------------------
% isawag, 07-03-2013
% modified on 29-10-2017
% last change 2018-08-13
% ----------------------------------------

clear pattern

pattern.data = zeros(size(fmri.dat,4),length(curr_roi));

for i = 1:size(fmri.dat,4)
    
    clear temp
    
    if permute == 1
        
        temp = fmri.dat{1}(:,:,:,Pseq(i));
        
    else
        
        temp = fmri.dat(:,:,:,i);
        
    end
    
    pattern.data(i,:) = temp(indices)';

end

pattern.tLabel = fmri.cond';

if RSAtrue == false; pattern.cvLabel = fmri.k'; end 

%pattern.sLabel = cat.ssfold{1}';

end
