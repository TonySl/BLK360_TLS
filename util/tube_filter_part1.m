function new_trunk_id=tube_filter_part1(temp_trunk_point,allpts,ind_trunk1,tube_z_slice1,tube_xysize1)

    tube_xysize=tube_xysize1;

    
    tube_x_left=temp_trunk_point(1)-tube_xysize/2;
    tube_x_right=temp_trunk_point(1)+tube_xysize/2;
    tube_y_left=temp_trunk_point(2)-tube_xysize/2;
    tube_y_right=temp_trunk_point(2)+tube_xysize/2;

    
    tube_z_slices=tube_z_slice1; % divide the tube into smaller vertical bins, to check if the points in the tube forms a line-similar 3d feature.  

    neigors_in_tube_ind=abs(allpts(:,1)-temp_trunk_point(1))<tube_xysize/2 & abs(allpts(:,2)-temp_trunk_point(2))<tube_xysize/2;

    Z_trunkpts_in_tube=allpts(ind_trunk1 & neigors_in_tube_ind,3);
    tube_zsize=range(Z_trunkpts_in_tube); %%%%%%%%%% adaptive tube z size
    tube_z_bottom=min(Z_trunkpts_in_tube);
    tube_z_floor=max(Z_trunkpts_in_tube);
    
    ind_trunk_in_tube=ind_trunk1(neigors_in_tube_ind); %%% need to check if each slice also has mode 1
    

    if sum(ind_trunk1(neigors_in_tube_ind))>floor(sum(neigors_in_tube_ind)/2)
        
        %%% wheter trunk points distributed along z 
        [z_distribution,bin_ind]=histc(Z_trunkpts_in_tube, tube_z_bottom:tube_zsize/tube_z_slices:tube_z_floor); % use this, not linespace
        

        each_bin_mode=nan(max(bin_ind),1); 
        for i_bin_ind=1:max(bin_ind)
            %if each bin has mode 1, then all the points in this tube are labled as 1, a high quality tube
            if isempty(ind_trunk_in_tube(bin_ind==i_bin_ind))
                continue
            end
%             each_bin_mode=[each_bin_mode mode(ind_trunk_in_tube(bin_ind==i_bin_ind))];
            each_bin_mode(i_bin_ind)=mode(ind_trunk_in_tube(bin_ind==i_bin_ind)); %%%%%% mode is time consuming
            
        end
        
        each_bin_mode(isnan(each_bin_mode))=[];
        

        if sum(z_distribution>0)>tube_z_slices-1 && sum(each_bin_mode==0)>0

            new_trunk_id=1;

        elseif sum(z_distribution>0)<= tube_z_slices-1 && sum(each_bin_mode==0)==0 % && range(Z_trunkpts_in_tube)>tube_zsize*0.8
            
            new_trunk_id=2;
%               new_trunk_id=find(neigors_in_tube_ind==1); %%%%% get all the points in the tube
            
        elseif sum(z_distribution>0)>tube_z_slices-1 && sum(each_bin_mode==0)==0
            
            new_trunk_id=3;
            
         else

            new_trunk_id=0;
            
        end
        
    else
        new_trunk_id=0;

    end


end
