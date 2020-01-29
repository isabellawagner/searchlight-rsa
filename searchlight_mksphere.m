
function [roi] = searchlight_mksphere(x,y,z,radius,vx_size)

%creates spherical ROIs with a certain radius

roi=[];

for a=(x-radius):(x+radius)
    
    for b=(y-radius):(y+radius)
        
        for c=(z-radius):(z+radius)
            
            if pdist([x,y,z;(a*vx_size(1)),(b*vx_size(2)),(c*vx_size(3))])<radius
                
                roi = [roi ; a,b,c];
              
            end
            
        end
        
    end
    
end

roi=unique(roi,'rows');

end
