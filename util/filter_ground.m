    
function ind_ground=filter_ground(pts_xyz,xyrange1,z_differ1)

    pts_xy=pts_xyz(:,1:2);
    [idx, ~] = rangesearch(pts_xy,pts_xy,xyrange1);

    min_z = cellfun(@(s) min(pts_xyz(s,3)), idx, 'UniformOutput', false);
    minZ_differ=pts_xyz(:,3)-cell2mat(min_z);
    ind_ground=minZ_differ<z_differ1;
    
end