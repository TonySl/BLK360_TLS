function new_trunk_id=tube_filter_part2(temp_trunk_point,trunkpts_afterfirstTube,tube_xysize1, height_thre) 

    %%%%% Further filter out short tubes
   
    tube_xysize=tube_xysize1;

    tube_x_left=temp_trunk_point(1)-tube_xysize/2;
    tube_x_right=temp_trunk_point(1)+tube_xysize/2;
    tube_y_left=temp_trunk_point(2)-tube_xysize/2;
    tube_y_right=temp_trunk_point(2)+tube_xysize/2;


    trunkpts_afterfirstTube_in_tube_ind=abs(trunkpts_afterfirstTube(:,1)-temp_trunk_point(1))<tube_xysize/2 & abs(trunkpts_afterfirstTube(:,2)-temp_trunk_point(2))<tube_xysize/2; % & trunkpts_afterfirstTube(:,3)>tube_z_bottom & trunkpts_afterfirstTube(:,3)<tube_z_floor;

    Z_trunkpts_in_tube=trunkpts_afterfirstTube(trunkpts_afterfirstTube_in_tube_ind,3);
    
    if range(Z_trunkpts_in_tube)< height_thre %tube_zsize/2 
        new_trunk_id=0;
        
    else
        
        new_trunk_id=1;
        
    end
   

end