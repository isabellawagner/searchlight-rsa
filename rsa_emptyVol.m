
function [OUTPUT] = rsa_emptyVol(subj,tLabel,refvol,dim)

% ----------------------------------------
% isawag, 2015-2018
% last change 2018-08-13
% ----------------------------------------

clear wsim; wsim = unique(tLabel);

clear bsim; bsim = combnk(unique(tLabel),2);

clear OUTPUT

% ----------------------------------------
% make volumes for all within similarities
% ----------------------------------------

c = 1; 

volname = [subj,'_searchlight_WS_'];

for i = 1:numel(wsim)
    
    OUTPUT.wsim.Vout(c) = refvol;
    
    OUTPUT.wsim.Vout(c).fname = [OUTPUT.wsim.Vout(c).fname(1:find(OUTPUT.wsim.Vout(c).fname == '/',1,'last')),volname,num2str(wsim(i)),'.nii'];
    
    OUTPUT.wsim.data_out(:,:,:,c) = nan(dim(1),dim(2),dim(3));
    
    c = c + 1;
    
end

% ----------------------------------------
% make volumes for all between similarites
% ----------------------------------------

c = 1;

%empty volumes for bsim

volname = [subj,'_searchlight_BS_'];

for i = 1:size(bsim,1)
    
    OUTPUT.bsim.Vout(c) = refvol;
    
    OUTPUT.bsim.Vout(c).fname = [OUTPUT.bsim.Vout(c).fname(1:find(OUTPUT.bsim.Vout(c).fname == '/',1,'last')),volname,num2str(bsim(i,1)),'_',num2str(bsim(i,2)),'.nii'];
    
    OUTPUT.bsim.data_out(:,:,:,c) = nan(dim(1),dim(2),dim(3));
    
    c = c + 1;
    
end

