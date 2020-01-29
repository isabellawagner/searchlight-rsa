
function [vox_proc, voxcount, coord_list] = util_countspheres(fmridat,show,dim,roi,ngreyvx)

%count the number of spheres on volume and make coordinate list

if show == true; disp('count rois ...'); end 

[vox_proc, voxcount] = deal(0);

coord_list = [];

clear check

for x = 1:dim(1)
    
    for y = 1:dim(2)
    
        for z = 1:dim(3)
        
            curr_roi = roi + repmat([x,y,z],length(roi(:,1)),1);
            
            if min(curr_roi(:) > 0) && max(curr_roi(:,1)) < dim(1) && max(curr_roi(:,2)) < dim(2) && max(curr_roi(:,3)) < dim(3)
            
                check = searchlight_greyCheck(fmridat.grey,curr_roi,ngreyvx);
                
                if check == 0 || check == 2 %more than 30 GM voxels in volume (or all are GM)
                    
                    coord_list = [coord_list; x y z];
                
                    voxcount = voxcount + 1;
                    
                    if show == true; disp(num2str(voxcount)); end 
                end
                
            end
            
        end
        
    end
    
end
