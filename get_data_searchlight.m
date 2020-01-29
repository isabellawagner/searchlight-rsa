
% ---------------------------------------
% get_data_searchlight
% ---------------------------------------
% isawag, isabella.wagner@univie.ac.at
% ---------------------------------------

%%

fmri.dat = [];

%%

cd(parm.dirs.betas); cd(['opi13-',currsubj])

clear files; files = dir([char(parm.mvpa.datatype),'*.img']);

%%

%get only the trials that are used

if strcmp(parm.mvpa.level,'Bsubj')

    fakeseq = [...
        ones(1,parms.ntrials/numel(parms.cond)), ...
        ones(1,parms.ntrials/numel(parms.cond))*2, ...
        ones(1,parms.ntrials/numel(parms.cond))*3, ...
        ones(1,parms.ntrials/numel(parms.cond))*4];
    
    idx = find(fakeseq == str2double(parms.train));
    
    start = min(idx); stop = max(idx); 
    
    if str2double(parms.train) == 4; stop = stop -1; end 
    
else
    
    start = 1; stop = parm.mvpa.ntrials; 
    
    if ismember(currsubj,{'c001','c002','c003','c005','c006','c007','c011','c012','c013','c015'});
        
        stop = stop -1;
        
    end
    
end

%%

fprintf('get fmri data ')

counter = 1; 

for f = start:stop 
    
    clear vol; vol = spm_vol(files(f).name);
    
    if strcmp(parm.mvpa.roi,'searchlight') && ~isfield(fmri,'refvol')
        
        fmri.refvol = vol;
        
    end
    
    %if parms.show == true; fprintf('.'); end 
    
    if strcmp(parm.mvpa.roi,'searchlight')
        
        %get whole-brain data
        
        fmri.dat(:,:,:,counter) = spm_read_vols(vol);
        
    else
        
        %else only load the coordinates we need roixyz
        
        fmri.dat(:,:,:,counter) = spm_get_data(vol,roixyz);
        
    end
    
    counter = counter + 1; fprintf('.')
    
end

fprintf(' done.\n')

%% if we use a searchlight, load also the grey matter mask

if strcmp(parm.mvpa.roi,'searchlight')

    cd(fullfile(parm.dirs.project,'subjects',['opi13-',currsubj],'trio','t1_mprage_we_norm','nifti'))
    
    %check first if T1 was already coregistered (necessary first step)
    
    clear file; file = dir('c1r2betawvols.nii');
    
    %reslice it to beta dimension if necessary

    if isempty(file)
        
        jobfile = cellstr(fullfile(parm.dirs.mscripts,'job_coreg_estresl.m'));
      
        jobs = repmat(jobfile, 1, 1);
        inputs = cell(2, 1);
        
        %first coregister normlized T1 (wvols.nii) to beta > r2betawvols.nii
        
        inputs{1, 1} = cellstr(fullfile(parm.dirs.betas,['opi13-',currsubj],'beta_0001.img'));    % Coregister: Estimate & Reslice: Reference Image - cfg_files
        inputs{2, 1} = cellstr(fullfile(cd,'wvols.nii'));                                   % Coregister: Estimate & Reslice: Source Image - cfg_files
        
        spm('defaults','fMRI')
        spm_jobman('serial', jobs, '', inputs{:});
        
        %then segment coregistered, normalized T1 > c1r2betawvols.nii
        
        jobfile = cellstr(fullfile(parm.dirs.mscripts,'job_segment.m'));
      
        jobs = repmat(jobfile, 1, 1);
        inputs = cell(1, 1);
        
        inputs{1, 1} = cellstr(fullfile(cd,'r2betawvols.nii'));
        
        spm('defaults','fMRI')
        spm_jobman('serial', jobs, '', inputs{:});
        
        file = dir('c1r2betawvols.nii');
        
    end

    [roi] = supp_loadroi(file.name,0.25);
    
    fmri.grey = roi.mask;

end

