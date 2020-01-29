
function run_searchlight_rsa(jobinput)

% ---------------------------------------
% isawag, isabella.wagner@univie.ac.at
% ---------------------------------------

%%

parm = jobinput{1};

currsubj = char(parm.info.subj);

clc

% ---------------------------------------
% what we need
% ---------------------------------------

addpath(genpath(parm.dirs.spm))

addpath(genpath(parm.dirs.svm))

% ---------------------------------------
% show 
% ---------------------------------------

if parm.info.show == true
    
    disp('--------------------------------------------------')
    
    if strcmp(parm.mvpa.level,'Wsubj')
        
        disp(currsubj)
        
    end
    
    disp('--------------------------------------------------')
    
end

% ---------------------------------------
% get data
% ---------------------------------------

if strcmp(parm.mvpa.roi,'searchlight')

    get_data_searchlight

end

% ---------------------------------------
% get labels
% ---------------------------------------

%labels, spm file has 1-4 conditions in order

fmri.cond = [];

for c = 1:max(parm.mvpa.ncnd)
    
    fmri.cond = [fmri.cond, ones(1,parm.mvpa.ntrials/numel(parm.mvpa.ncnd)) * parm.mvpa.ncnd(c)];
    
end

% last trial is missing

fmri.cond(end) = [];

% ---------------------------------------
% set up cross-validation (for SVM)
% ---------------------------------------

%generally, values that were acquired within the same fold should
%also end up in the same fold. otherwise, two samples that highly correlate
%in time but are in respective training and test sets, might lead to an
%overestimation of decoding performance

%15-fold cross-validation

fmri.k = [];

if strcmp(parm.mvpa.type,'SVM')
    
    
    if strcmp(parm.mvpa.level,'Wsubj')
        
        for c = 1:max(parm.mvpa.ncnd); fmri.k = [fmri.k, 1:parm.mvpa.ntrials/max(parm.mvpa.ncnd)]; end
      
        % last trial is missing
        
        fmri.k(end) = [];
        
    else
        
        %     fmri.k = 1:numel(parm.info.subjlist{1});
        %
        %     %how many conditions do we include?
        %
        %     csize = numel(str2double(parms.train));
        %
        %     %also make a group vector ('labels')
        %
        %     fmri.group = [];
        %
        %     for gr = 1:stop
        %
        %         if strcmp(subj{1}{gr}(1),'c')
        %
        %             fmri.group = [fmri.group ones(1,numel(find(fmri.cond == str2double(parms.train))))]; %parms.ntrials/numel(parms.cond))*1];
        %
        %         elseif strcmp(subj{1}{gr}(1),'p')
        %
        %             fmri.group = [fmri.group ones(1, numel(find(fmri.cond == str2double(parms.train))))*2];%parms.ntrials/numel(parms.cond))*2];
        %
        %         end
        %
        %     end
        
    end
    
end

% ---------------------------------------
% make empty volume 
% ---------------------------------------

tmp = spm_imatrix(fmri.refvol.mat);
parm.mvpa.vx_size = abs(tmp(7:9));
parm.mvpa.dim = fmri.refvol.dim;

if strcmp(parm.mvpa.type,'SVM')
    
    parm.mvpa.volname = ['searchlight_',currsubj];
    
    [VOL,DAT] = supp_emptyvol(currsubj,fmri,parm.mvpa.dim ,parm.mvpa.volname);
    
else
    
    %parms = parm; 
    subj = currsubj;
    tLabel = fmri.cond;
    refvol = fmri.refvol;
    dim = parm.mvpa.dim;
    %P = 1; 
    
    [OUTPUT] = rsa_emptyVol(subj,tLabel,refvol,dim);
    
    DATdim = nan(parm.mvpa.dim(1),parm.mvpa.dim(2),parm.mvpa.dim(3));
    
end

% ---------------------------------------
% make a sphere
% ---------------------------------------

% first make a searchlight at coordinates 0 0 0 to start with 
% test dimensions and use this as starting point to go to coord_list(nsl,:)
% below

[roi] = searchlight_mksphere(0,0,0,str2double(parm.mvpa.radius),parm.mvpa.vx_size);

parm.mvpa.roi_size = size(roi,1);

% then, make a list of all possible searchlights that contain > x GM voxels

% [vox_proc, voxcount] = util_countspheres(fmridat,show,dim,roi,ngreyvx)

[vox_proc, voxcount, coord_list] = util_countspheres(fmri,parm.info.show,parm.mvpa.dim,roi,parm.mvpa.ngreyvx);

parm.mvpa.nrois = voxcount;

parm.mvpa.coord_list = coord_list; 

% ---------------------------------------
% run the searchlight
% ---------------------------------------

loop_start = clock;

if strcmp(parm.mvpa.type,'SVM')

dataout = NaN(size(coord_list,1),1);

elseif strcmp(parm.mvpa.type,'RSA')
    
    dataout.wsim = NaN(size(coord_list,1),numel(unique(fmri.cond)));
    
    pairs = combnk(unique(fmri.cond),2);
    
    dataout.bsim = NaN(size(coord_list,1),size(pairs,1));
    
end

test = 1:size(coord_list,1);

for nsl = 1:size(coord_list,1)
    
    % show progress
    
    vox_proc = test(nsl);
    
    if parm.info.show == true; disp([num2str(vox_proc),' / ',num2str(voxcount)]); end
    
    % we take our coord_list from above with all possible searchlights that
    % contain more than x gm voxels (check == 2), start from 0 0 0 
    
    curr_roi = roi + repmat(coord_list(nsl,:),length(roi(:,1)),1);
    
    % get the voxel coordinates of the 251 voxels contained in 8 mm sphere
    
    indices = sub2ind(size(fmri.grey),curr_roi(:,1),curr_roi(:,2),curr_roi(:,3));
    
    % check if we have to do some permutation
    
    if parm.mvpa.permute == false; Pseq = 0; end
    
    %get data in shape
    
    if strcmp(parm.mvpa.type,'SVM'); RSAtrue = false; else RSAtrue = true; end 
    
    [pattern] = searchlight_getdata(fmri,curr_roi,indices,parm.mvpa.permute,Pseq,RSAtrue);
    
    %remove NaNs
    
    pattern.data(isnan(pattern.data)) = 0;
    
    if strcmp(parm.mvpa.type,'SVM')
        
        [acc] = mvpa_classify(pattern.data,pattern.tLabel',pattern.cvLabel',parm);
        
        % assign average accuracy
        
        dataout(nsl) = mean(acc);
        
        %if parm.info.show == true; disp(mean(acc)); end
        
    elseif strcmp(parm.mvpa.type,'RSA')
        
        % zscore data within run 
        
        pattern.zdata = zscore(pattern.data,[],2);
        
        % correlate 
        
        clear corrmat; corrmat.data = corr(pattern.zdata');
        
        % Fisher's z transform Pearson's r 
        corrmat.fisher = 0.5 * log((1+corrmat.data) ./ (1-corrmat.data));
        
        % NaN lower part of matrix (incl. diagonal)
        corrmat.fisher(corrmat.fisher == tril(corrmat.fisher,0)) = NaN;
        
        % get similarity values
        [rwith, rbetw] = rsa_sim(pattern,corrmat);
        
        %assign wsim to center voxel
        for i = 1:numel(rwith)
            
            dataout.wsim(nsl,i) = rwith(i);
            
        end
        
        %assign bsim to center voxel
        for i = 1:numel(rbetw)
            
            dataout.bsim(nsl,i) = rbetw(i);
            
        end
        
    end
    
end

end_time = clock;
duration = etime(end_time,loop_start);

%%

% ---------------------------------------
% get back in shape
% ---------------------------------------

% go through accuracy list and assign to coordinate

if strcmp(parm.mvpa.type,'SVM')
    
    for nsl = 1:size(dataout,1)
        
        DAT(coord_list(nsl,1),coord_list(nsl,2),coord_list(nsl,3)) = dataout(nsl);
        
    end
    
elseif strcmp(parm.mvpa.type,'RSA')
    
    clear DAT
    
    % wsim
    
    for i = 1:numel(rwith)
        
        disp(num2str(i)); disp('----------------------');
        
        for nsl = 1:size(dataout.wsim,1)
            
            disp(num2str(nsl))
            
            if nsl == 1; DAT.wsim{i} = DATdim; end 
            
            DAT.wsim{i}(coord_list(nsl,1),coord_list(nsl,2),coord_list(nsl,3)) = dataout.wsim(nsl,i);
            
        end
        
    end
    
    % bsim
    
    for i = 1:numel(rbetw)
        
        disp(num2str(i)); disp('----------------------');
        
        for nsl = 1:size(dataout.bsim,1)
            
            disp(num2str(nsl))
            
            if nsl == 1; DAT.bsim{i} = DATdim; end
            
            DAT.bsim{i}(coord_list(nsl,1),coord_list(nsl,2),coord_list(nsl,3)) = dataout.bsim(nsl,i);
            
        end
        
    end
    
end

%% save everything 

cd(parm.dirs.out)

save([subj,'_searchlight_RSA.mat']);

%% write results to volume

if strcmp(parm.mvpa.type,'SVM')
    
    for i = 1:numel(VOL); VOL(i) = spm_write_vol(VOL(i),DAT(:,:,:,i)); end
    
elseif strcmp(parm.mvpa.type,'RSA')
    
    % wsim
    
    for i = 1:numel(OUTPUT.wsim.Vout)
        
        OUTPUT.wsim.Vout(i) = spm_write_vol(OUTPUT.wsim.Vout(i),DAT.wsim{i});
        
    end
    
    % bsim
    
    for i = 1:numel(OUTPUT.bsim.Vout)
        
        OUTPUT.bsim.Vout(i) = spm_write_vol(OUTPUT.bsim.Vout(i),DAT.bsim{i});
        
    end
    
end

%%

if parm.info.show == true; disp('done.'); end 
   

