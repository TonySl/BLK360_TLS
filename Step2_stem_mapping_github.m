%% Stem detection from TLS data. Version 0.0

%%% Shengli Tao. Juin 2020.
%%% See Tao et al. 2020. Methods in Ecology and Evolution (submitted) for details.

%%
clear
clc
warning off
addpath('./util') 


%%%% Creat a new folder to store the output files %%%%
mkdir TLidar_output


%% Input files & parameters

%%%% Several txt files are required as input files. Each txt file contains points in a specific height layers. 
%%%% In Tao et al. MEE, seven txt files were used as input, each containing points spaning a 50cm height duration(namely, 0.5-1,1-1.5,1.5-2,2-2.5,2.5-3,3-3.5,3.5-4). 
%%%% Each file has three colums (x, y, z) seperated by comma, without header line. 
%%%% Points between 0.5-1.5m height span will be used to detect small trees which are close to the ground. 
%%%% Points between 1.5-4m height span will be used for detecting other larger stems.

%%%% Use height bins each spanning 1m instead of 0.5m (such as 0.5-1.5m, 1.5-2.5m, 2.5-3.5m), we will have less noisy points
%%%% but possibly less stems detected, because we are putting a high demand on the vertical
%%%% continuity of stem points.

%%%% Please order the file from left to right according to their distance to ground (from low to high)
%%%% Please subsmaple the points to reduce computation burden. In Tao et al, the point cloud was subsampled using
%%%% CloudCompare(danielgm.net/cc), to achieve a pairwise point distance > 2 cm.

dir_input='./'; %%%%% folder containing input file

bin_filenames_05_4m={'petitplateau_bin_050_100_2cm.txt','petitplateau_bin_100_150_2cm.txt','petitplateau_bin_150_200_2cm.txt', ...
                    'petitplateau_bin_200_250_2cm.txt','petitplateau_bin_250_300_2cm.txt','petitplateau_bin_300_350_2cm.txt','petitplateau_bin_350_400_2cm.txt'};  
bin_filenames_id=1:length(bin_filenames_05_4m); %%%% Each bin will give an id from 1 to N, ordered from low to high, for visulization.

bin_id_05_4m=[1 1 2 2 2 2 2]; %%% ID 1 means this height layer will be used to detet small trees which grow close to the ground. 
                              %%% For example, [2 2 2 2 2 2 2] won't detect very small trees. [1 1 1 1 1 1 1] is not suggested (detect small trees with high bin point cloud)

                              
%%%%%%% Tube filter %%%%%%%
tube_xysize_default=0.1; %%%% 10cm tube size (width and length) as default. 15cm gave similar results.
tube_z_slices_default=5; %%%% The tube will be cut into 5 smaller bins to check whether stem points dominated in each of the small bins.
                         %%%% Higher value means a high demand on the vertical continuity of the stem points


%%%%%%% Segmentation %%%%%%%
distThreshold=0.1; % 10cm as default. Higher value will cause more under-segmentation errors, lower value more over-segmentation errors, possibly.


%% Read txt data, calculate surface normal, performe the Tube filter, get the stem points
if 1  %%%% 41525.112 s, ~11h 30min for seven height bins 
    tic
    
    disp('Read txt data, calculate surface normal, perform the Tube filter.............................')

    total_z_range_pts_higherbins=0;
    total_z_range_pts_lowerbins=0;

    number_of_pts_total_05_15m=0;
    number_of_pts_each_05_15m=[];
    
    number_of_pts_total_15_4m=0;
    number_of_pts_each_15_4m=[];

    for i_file=1:length(bin_filenames_05_4m)

        plot_bin_data=dlmread(bin_filenames_05_4m{i_file},',', 0, 0);
        % load(strcat(dir_input,bin_filenames_05_4m{i_file}))
        
        disp(bin_filenames_05_4m{i_file})

        each_bin_height=range(plot_bin_data(:,3));
        

        height_thre_secfilter=each_bin_height*0.8;


        if bin_id_05_4m(i_file)==1
            total_z_range_pts_lowerbins=total_z_range_pts_lowerbins+each_bin_height;
        else
            total_z_range_pts_higherbins=total_z_range_pts_higherbins+each_bin_height;
        end

        %%%% cut the data into ~2m*2m tiles. Didn't consider buffer here to save computation time %%%%
        original_pts_id=1:size(plot_bin_data,1);

        minx=min(plot_bin_data(:,1));
        maxx=max(plot_bin_data(:,1));
        miny=min(plot_bin_data(:,2));
        maxy=max(plot_bin_data(:,2));

        % bin the xy data of all pts into ~2m*2m tiles
        nx = round(min((maxy-miny)/2,(maxx-minx)/2)); 
        ny = nx; 

        xedge = linspace(minx,maxx,nx+1); 
        yedge = linspace(miny,maxy,ny+1);

        xbin_allpts = discretize(plot_bin_data(:,1), xedge);  % xbin index for all the pts
        ybin_allpts = discretize(plot_bin_data(:,2), yedge);  % ybin index for all the pts

%         groups = splitapply( @(x){x}, trunk_pts_withZ_05_4m, labels_10cmdist_05_4m);

%         Trunk_ind_surface_normal=nan(length(plot_bin_data),1); 
%         Trunk_ind_surface_normal_tubefilter1=nan(length(plot_bin_data),1);
        Trunk_ind_surface_normal_tubefilter2=nan(length(plot_bin_data),1);

        for i=1:length(xedge)-1

            disp(100*i/(length(xedge)-1))


            for j=1:length(yedge)-1

                pts_withintile_ind = xbin_allpts==i & ybin_allpts ==j;

                if sum(pts_withintile_ind)<10 %%%%% Needs at least 10 points inside each tile
                    continue
                end

                pts_tile=plot_bin_data(pts_withintile_ind,:); %%%%%%%%% time consuming 

                normals_MATLAB = pcnormals(pointCloud(pts_tile),10);

                v=[0 0 1];
                normals_MATLAB=num2cell(normals_MATLAB,2);
                angles_matlab = cellfun(@(s) atan2(norm(cross(s,v)),dot(s,v)).*(180/pi), normals_MATLAB, 'UniformOutput', false);
                ind_trunk=abs(cell2mat(angles_matlab)-90)<10; %%%% potential stem points, the surface normals of which are perpendicular to ground with a tolerance of 10 degree

                %%%%%%%%%%%%% make sure the 10 neighbors of every 'potential stem point' are also 'potential stem points' %%%%%%%%%%%%
                Idx = knnsearch(pts_tile,pts_tile,'K',10);
                Idx(:,1)=[];% remove itself

%                 neighbor10_trunkInd = cellfun(@(s) mode(ind_trunk(s)), num2cell(Idx,2), 'UniformOutput', false); %%% mode function is slow                 
                neighbor10_trunkInd = cellfun(@(s) sum(ind_trunk(s)), num2cell(Idx,2), 'UniformOutput', false);
%                 clear Idx


%                 neighbor10_trunkInd_tile=cell2mat(neighbor10_trunkInd);  mode function is slow  
                neighbor10_trunkInd_tile=cell2mat(neighbor10_trunkInd)>4; %%% 9 neighbors, >4 is the same with mode function

                ind_trunk2=neighbor10_trunkInd_tile==1;

                trunkpts_by_surfacenormal=pts_tile(ind_trunk2,[1 2 3]);

                %%%%% tube filter %%%%%
                trunk_points_cell=num2cell(trunkpts_by_surfacenormal,2); %%% Convert to cell array to save computation time
                tempc = cellfun(@(V) tube_filter_part1(V, pts_tile, neighbor10_trunkInd_tile,tube_z_slices_default,tube_xysize_default), trunk_points_cell,'Uniform', 0); %%% yy  


                new_trunk_id=cell2mat(tempc);
                firstTubefilter_trunkid_intile=zeros(length(ind_trunk2),1);
                firstTubefilter_trunkid_intile(ind_trunk2)=new_trunk_id;

                trunkpts_tile=pts_tile(logical(firstTubefilter_trunkid_intile),[1 2 3]);
                if isempty(trunkpts_tile)
                    continue            
                end

                trunk_points_cell=num2cell(trunkpts_tile,2);

                %%%% Filter out some potential stem points which don't have vertical neighbors 
                tempc2 = cellfun(@(V) tube_filter_part2(V, trunkpts_tile, tube_z_slices_default, height_thre_secfilter), trunk_points_cell,'Uniform', 0);


                new_trunk_id2=cell2mat(tempc2);
                secTubefilter_trunkid_intile=zeros(length(firstTubefilter_trunkid_intile),1);
                secTubefilter_trunkid_intile(logical(firstTubefilter_trunkid_intile))=new_trunk_id2;
                Trunk_ind_surface_normal_tubefilter2(original_pts_id(pts_withintile_ind))=secTubefilter_trunkid_intile; %%%%%%%%% time consuming 


            end


        end

%         sum(Trunk_ind_surface_normal_tubefilter2==1)

        %%%%%%%%% save results as temp files  %%%%%%%%%
        potential_trunk_pts=plot_bin_data(Trunk_ind_surface_normal_tubefilter2==1,1:3); % x y z
        save(strcat('./TLidar_output/tube_filter_',bin_filenames_05_4m{i_file}(1:end-4)),'potential_trunk_pts','-v7.3')

        if bin_id_05_4m(i_file)==1            
            number_of_pts_each_05_15m=[number_of_pts_each_05_15m size(potential_trunk_pts,1)];
            number_of_pts_total_05_15m=number_of_pts_total_05_15m+size(potential_trunk_pts,1);
        else    
            number_of_pts_total_15_4m=number_of_pts_total_15_4m+size(potential_trunk_pts,1);
            number_of_pts_each_15_4m=[number_of_pts_each_15_4m size(potential_trunk_pts,1)];            
        end

    end
    
    save('./TLidar_output/bin_statistics','total_z_range_pts_lowerbins','total_z_range_pts_higherbins','number_of_pts_total_15_4m','number_of_pts_total_05_15m','number_of_pts_each_05_15m','number_of_pts_each_15_4m')

    toc
end


%% DistSeg on 2D dimension, to segment the stem points into individual components (i.e. putative stem)
if 1 %%%% 

    load('./TLidar_output/bin_statistics','total_z_range_pts_lowerbins','total_z_range_pts_higherbins','number_of_pts_total_15_4m','number_of_pts_total_05_15m','number_of_pts_each_05_15m','number_of_pts_each_15_4m')
    
    height_threhold_05_15m=total_z_range_pts_lowerbins*0.8; % Two bins (totally 1m in height span) will be used to detect small trees which are close to ground. Each_bin_height*length(bin_filenames_05_15m)*0.8
    height_threhold_15_4m=total_z_range_pts_higherbins*0.8; % each_bin_height*length(bin_filenames_15_4m)*0.8
   
        
    trunk_pts_withZ_15_4m=nan(number_of_pts_total_15_4m,4); %%% Each column is: x y z binID
    trunk_pts_withZ_05_15m=nan(number_of_pts_total_05_15m,4); %%% Each column is: x y z binID
    
    number_of_pts_each_05_15m=[1 number_of_pts_each_05_15m];
    number_of_pts_each_15_4m=[1 number_of_pts_each_15_4m];
    
    number_of_pts_each_05_15m=cumsum(number_of_pts_each_05_15m);
    number_of_pts_each_15_4m=cumsum(number_of_pts_each_15_4m);
    
    i_1=0;
    i_2=0;
    
    for i_file=1:length(bin_filenames_05_4m)

        load(strcat('./TLidar_output/tube_filter_',bin_filenames_05_4m{i_file}(1:end-4)),'potential_trunk_pts')

        if bin_id_05_4m(i_file)==1
            i_1=i_1+1;
            trunk_pts_withZ_05_15m(number_of_pts_each_05_15m(i_1):number_of_pts_each_05_15m(i_1+1)-1,1:3)=potential_trunk_pts;
            trunk_pts_withZ_05_15m(number_of_pts_each_05_15m(i_1):number_of_pts_each_05_15m(i_1+1)-1,4)= ones(size(potential_trunk_pts,1),1)*i_file;
        else
            i_2=i_2+1;
            trunk_pts_withZ_15_4m(number_of_pts_each_15_4m(i_2):number_of_pts_each_15_4m(i_2+1)-1,1:3)=potential_trunk_pts;
            trunk_pts_withZ_15_4m(number_of_pts_each_15_4m(i_2):number_of_pts_each_15_4m(i_2+1)-1,4)= ones(size(potential_trunk_pts,1),1)*i_file; %%% x y z bin_id
        end

    end
    
%     if sum(sum(isnan(trunk_pts_withZ_05_15m),1))>0 | sum(sum(isnan(trunk_pts_withZ_15_4m),1))>0
%         disp('wrong rows')
%     end
    
    %%
    
    if isempty(trunk_pts_withZ_15_4m)
        disp('Need at least one bin above 1.5m height')
        return
    end
    
    trunk_pts_noZ_15_4m=trunk_pts_withZ_15_4m(:,1:3);
    trunk_pts_noZ_15_4m(:,3)=0;
    
    disp('distance segmentation for higher bins ..............................')
    tic
    [labels_10cmdist_15_4m,numClusters_10cmdist_15_4m] = pcsegdist(pointCloud(trunk_pts_noZ_15_4m),distThreshold); %%% labels_10cmdist_15_4m: stemID
    clear trunk_pts_noZ_15_4m
    
    %%% height of each segmented component
    groups = splitapply( @(x){x}, trunk_pts_withZ_15_4m(:,3), labels_10cmdist_15_4m); %%%%% group the z value of all points according to the cluster ID
    stem_bin_height_15_4m = cellfun(@range,groups);
    stem_bin_height_15_4m=[stem_bin_height_15_4m (1:max(labels_10cmdist_15_4m))']; %%% height, stemID
    
    
    in2=stem_bin_height_15_4m(:,1)>height_threhold_15_4m; %%%%%%%%%%%%%%%%% threholding the height here. Some short components are not likely to be real stems.
    in2_id=stem_bin_height_15_4m(in2,2);

    trunk_pts_withZ_15_4m = trunk_pts_withZ_15_4m(ismember(labels_10cmdist_15_4m,in2_id),:);
    labels_10cmdist_15_4m = labels_10cmdist_15_4m(ismember(labels_10cmdist_15_4m,in2_id));
    
    %%% relabel the stem ID from 1 to N consectively.    
    [~,~,labels_10cmdist_15_4m_new]=unique(labels_10cmdist_15_4m); clear labels_10cmdist_15_4m
    labels_10cmdist_15_4m=labels_10cmdist_15_4m_new; clear labels_10cmdist_15_4m_new
    
%     save('./TLidar_output/Petit_plateau_2cmsubsample_dist10cmseg_15_4m','labels_10cmdist_15_4m','trunk_pts_withZ_15_4m','-v7.3')
    toc
    
%     fileID = fopen('./TLidar_output/Petit_plateau_2cmsubsample_dist10cmseg_15_4m_stems.txt','w');
%     fprintf(fileID,'%.3f,%.3f,%.3f,%i\n',[trunk_pts_withZ_15_4m(:,1:3) labels_10cmdist_15_4m]');
%     fclose(fileID);
    
    %%
    if isempty(trunk_pts_withZ_05_15m) %%%% if no lower bin, no small tree will be detected
        disp('no small tree will be detected')
    else
        
        disp('distance segmentation for lower bins ..............................')
        tic
        trunk_pts_noZ_05_15m=trunk_pts_withZ_05_15m(:,1:3);
        trunk_pts_noZ_05_15m(:,3)=0;

        [labels_10cmdist_05_15m,numClusters_10cmdist_05_15m] = pcsegdist(pointCloud(trunk_pts_noZ_05_15m),distThreshold);
        clear trunk_pts_noZ_05_15m

        %%% height duration of each segmented component
        groups = splitapply( @(x){x}, trunk_pts_withZ_05_15m(:,3), labels_10cmdist_05_15m); %%%%% group the z value of all points according to the cluster ID
        stem_bin_height_05_15m = cellfun(@range,groups);
        stem_bin_height_05_15m=[stem_bin_height_05_15m (1:max(labels_10cmdist_05_15m))']; %%% height, stemID


        in2=stem_bin_height_05_15m(:,1)>height_threhold_05_15m; %%%%%%%%%%%%%%%%% threholding the height here. Some short components are not likely to be real stems.
        in2_id=stem_bin_height_05_15m(in2,2);

        trunk_pts_withZ_05_15m = trunk_pts_withZ_05_15m(ismember(labels_10cmdist_05_15m,in2_id),:);
        labels_10cmdist_05_15m = labels_10cmdist_05_15m(ismember(labels_10cmdist_05_15m,in2_id));

        %%% relabel the stem ID from 1 to N consectively.    
        [~,~,labels_10cmdist_05_15m_new]=unique(labels_10cmdist_05_15m); clear labels_10cmdist_05_15m
        labels_10cmdist_05_15m=labels_10cmdist_05_15m_new; clear labels_10cmdist_05_15m_new

%         save('./TLidar_output/Petit_plateau_2cmsubsample_dist10cmseg_05_15m','labels_10cmdist_05_15m','trunk_pts_withZ_05_15m','-v7.3')
        toc
        
%         fileID = fopen('./TLidar_output/Petit_plateau_2cmsubsample_dist10cmseg_05_15m_stems.txt','w');
%         fprintf(fileID,'%.3f,%.3f,%.3f,%i\n',[trunk_pts_withZ_05_15m(:,1:3) labels_10cmdist_05_15m]');
%         fclose(fileID);
        
    end
    
    %%
%     disp('Combining same stems..............')
    
    if isempty(trunk_pts_withZ_05_15m) %%%% if no lower bin, no need to combine
        
        trunk_pts_withZ_05_4m=trunk_pts_withZ_15_4m; clear trunk_pts_withZ_15_4m
        labels_10cmdist_05_4m=labels_10cmdist_15_4m; clear labels_10cmdist_15_4m
    
    else        
    
        tic
        trunk_pts_noZ_15_4m=trunk_pts_withZ_15_4m(:,1:3);
        trunk_pts_noZ_15_4m(:,3)=0;
        trunk_pts_noZ_05_15m=trunk_pts_withZ_05_15m(:,1:3);
        trunk_pts_noZ_05_15m(:,3)=0;

        trunk_pts_noZ_05_4m=[trunk_pts_noZ_05_15m;trunk_pts_noZ_15_4m]; %%% x,y   
        clear trunk_pts_noZ_05_15m   trunk_pts_noZ_15_4m
        trunk_pts_withZ_05_4m=[trunk_pts_withZ_05_15m;trunk_pts_withZ_15_4m]; %%% x,y,z,bin_id
        clear trunk_pts_withZ_05_15m   trunk_pts_withZ_15_4m

        disp('distance segmentation..............................')
        [labels_10cmdist_05_4m,numClusters_10cmdist_05_4m] = pcsegdist(pointCloud(trunk_pts_noZ_05_4m),distThreshold);
    
    end

    %%% height of each segmented component
    groups = splitapply( @(x){x}, trunk_pts_withZ_05_4m(:,3), labels_10cmdist_05_4m); %%%%% group the z value of all points according to the cluster ID
    stem_bin_height_05_4m = cellfun(@range,groups);
    stem_bin_height_05_4m=[stem_bin_height_05_4m (1:max(labels_10cmdist_05_4m))']; %%% height, stemID

    %% check again the stem height span   
    in2=stem_bin_height_05_4m(:,1)>height_threhold_05_15m; %%%%%%%%%%%%%%%%% threholding the height here. Some short components are not likely to be real trees.
    in2_id=stem_bin_height_05_4m(in2,2);

    trunk_pts_withZ_05_4m = trunk_pts_withZ_05_4m(ismember(labels_10cmdist_05_4m,in2_id),:);
    labels_10cmdist_05_4m = labels_10cmdist_05_4m(ismember(labels_10cmdist_05_4m,in2_id));
    
    %%% relabel the stem ID from 1 to N consectively.    
    [~,~,labels_10cmdist_05_4m_new]=unique(labels_10cmdist_05_4m); clear labels_10cmdist_05_4m
    labels_10cmdist_05_4m=labels_10cmdist_05_4m_new; clear labels_10cmdist_05_4m_new
    
    
    disp('DBH estimation.....')
    groups = splitapply( @(x){x}, trunk_pts_withZ_05_4m(:,1:4), labels_10cmdist_05_4m); %%%%% 'groups' is a cell_array, each cell contains the points of a stem (x,y,z,bin_id)
    groups_withDBH=DBH_TLS(groups);   
    trunk_pts_withZ_05_4m_binID_stemID_DBH=cell2mat(groups_withDBH); clear groups_withDBH groups

    
    
    save('./TLidar_output/Petit_plateau_2cmsubsample_dist10cmseg_05_4m','labels_10cmdist_05_4m','trunk_pts_withZ_05_4m','trunk_pts_withZ_05_4m_binID_stemID_DBH','height_threhold_05_15m','height_threhold_15_4m','-v7.3')

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fileID = fopen(strcat('./TLidar_output/all_stems.txt'),'a');   
    str{1}='x';
    str{2}='y';
    str{3}='z';
    str{4}='bin_id';
    str{5}='stem_id';
    str{6}='DBH';

    fprintf(fileID,'%s,',str{1:5});
    fprintf(fileID,'%s\n',str{6});


    fprintf(fileID,'%.3f,%.3f,%.3f,%i,%i,%.3f\n',trunk_pts_withZ_05_4m_binID_stemID_DBH(:,1:6)'); %%% x, y, z, bin_id, stem_id, DBH
    fclose(fileID);

    toc


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% visually check the stems with a DBH > 50cm to find and correct potential over-segmentation errors %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end


%% Use DBSCAN to automatically correct under-segmentation errors (i.e., two stems were segmented as one stem)
%%% At least 25 points in each cluster, and 2cm distance between clusters. Need to tune the parameter based on forest types and point density.
if 1  %%% 35min to run
    
     disp('DBSCAN to find stems that need further segmentation.....')     
     tic
    
     load('./TLidar_output/Petit_plateau_2cmsubsample_dist10cmseg_05_4m','labels_10cmdist_05_4m','trunk_pts_withZ_05_4m','height_threhold_05_15m','height_threhold_15_4m')     
     %%% trunk_pts_withZ_05_4m: x, y, z, bin_id, 
     %%% labels_10cmdist_05_4m: stem_id
    
     stem_id=unique(labels_10cmdist_05_4m);
     
     %%% Group the stems based on their ID. Much faster than indexing the ID within a for loop
     groups = splitapply( @(x){x}, trunk_pts_withZ_05_4m, labels_10cmdist_05_4m); %%%%% 'groups' is a cell_array, each cell contains the points of a stem (x,y,z,bin_id)

     counter_dbscan_15_4m=0; %%%%% count how many stems needs further segmentation

     stem_dbscan_15_4m=[];    
     stemID_15_4m_need_DBSCAN=[];
     
     existing_id_max=max(stem_id); %%%% will give the newly found stems (by DBSCAN) a new ID
     
     for i3=1:length(groups)
         
        if mod(i3,500)==0
            disp(strcat(num2str(100*i3/length(groups)),'%.....'))
        end
        
        %%%%%%%%%%%%%%%%%%%%%%% 90% time is consumed here, indexing large dataset is super time-consuming %%%%%%%%%%%%%%%%%%%%%%
%         pts_stem=trunk_pts_withZ_05_4m(labels_10cmdist_05_4m==stem_id(i3),:);

        %%%% much faster
        pts_stem=groups{i3,1};
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [class,type]=dbscan(pts_stem(:,1:2),25,0.02); % at least 25 points in each cluster, and 2cm distance between clusters

        if max(class)>1 %%% class of -1 means noise points, 1 means one cluster, 2 means the second cluster....
            
            stemID_15_4m_need_DBSCAN=[stemID_15_4m_need_DBSCAN;stem_id(i3)];
            
            counter_dbscan_15_4m=counter_dbscan_15_4m+1;
   
            stem_height=range(pts_stem(:,3));      
            dbscan_new=[];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check each component %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i=1:max(class)  

                dbscan_component_xyzbinID=pts_stem(class==i,:);

                dbscan_component_height=range(dbscan_component_xyzbinID(:,3));

                flag_height_component=dbscan_component_height<stem_height*2/3; %%%%% filter out the DBCAN detected clusters with a small z range (not likely to be a true stem)

                if ~flag_height_component

                    existing_id_max=existing_id_max+1;

                    dbscan_componentxyz_newid=[dbscan_component_xyzbinID ones(size(dbscan_component_xyzbinID,1),1)*existing_id_max];

                    dbscan_new=[dbscan_new;dbscan_componentxyz_newid]; %%% x, y, z, bin_id, stem_id
                    
                end

            end

        
            stem_dbscan_15_4m=[stem_dbscan_15_4m;dbscan_new];  
            
        end
        
        

     end
        
     
%      disp(counter_dbscan_15_4m)
     
    %%%%% replace the stems with newly detected stems by DBSCAN 
    stem_05_4m_xyz_binID_stemID=[trunk_pts_withZ_05_4m labels_10cmdist_05_4m]; %%% x, y, z, bin_id, stem_id
    clear trunk_pts_withZ_05_4m
    stem_05_4m_xyz_binID_stemID(ismember(stem_05_4m_xyz_binID_stemID(:,5),stemID_15_4m_need_DBSCAN),:)=[]; %%%%% remove dbscan detected trees    
    stem_05_4m_xyz_binID_stemID=[stem_05_4m_xyz_binID_stemID;stem_dbscan_15_4m]; %%%%% add dbscan corrected trees


    %%% relabel the stem ID from 1 to N consectively.    
    [~,~,id_new]=unique(stem_05_4m_xyz_binID_stemID(:,5)); 
    stem_05_4m_xyz_binID_stemID(:,5)=id_new; clear id_new


    disp('DBH estimation.....')
    groups = splitapply( @(x){x}, stem_05_4m_xyz_binID_stemID(:,1:4), stem_05_4m_xyz_binID_stemID(:,5)); %%%%% 'groups' is a cell_array, each cell contains the points of a stem (x,y,z,bin_id)
    groups_withDBH=DBH_TLS(groups);   
    stem_05_4m_xyz_binID_stemID_DBH=cell2mat(groups_withDBH); clear groups_withDBH groups


    save('./TLidar_output/Petit_plateau_2cmsubsample_dist10cmseg_05_4m_DBSCAN','stem_05_4m_xyz_binID_stemID','stem_05_4m_xyz_binID_stemID_DBH','height_threhold_05_15m','height_threhold_15_4m','-v7.3')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fileID = fopen(strcat('./TLidar_output/all_stems_DBSCAN_corrected.txt'),'a');   
    str{1}='x';
    str{2}='y';
    str{3}='z';
    str{4}='bin_id';
    str{5}='stem_id';
    str{6}='DBH';

    fprintf(fileID,'%s,',str{1:5});
    fprintf(fileID,'%s\n',str{6});

    fprintf(fileID,'%.3f,%.3f,%.3f,%i,%i,%.3f\n',stem_05_4m_xyz_binID_stemID_DBH(:,1:6)'); %%% x, y, z, bin_id, stem_id, DBH
    fclose(fileID);

    toc


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% visually check the stems with a DBH > 50cm to find and correct potential over-segmentation errors %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
end


%% Further detect bad quality stems. To be improved.......
if 0  %%%% ~5min to run
    
    disp('Refine the stems...')
    tic
    
    load('./TLidar_output/Petit_plateau_2cmsubsample_dist10cmseg_05_4m_DBSCAN','stem_05_4m_xyz_binID_stemID','height_threhold_05_15m','height_threhold_15_4m')
     
    size1=size(stem_05_4m_xyz_binID_stemID,1);
    %%% Group the stems based on their ID. Much faster than indexing the ID within a for loop
    groups = splitapply( @(x){x}, stem_05_4m_xyz_binID_stemID(:,1:4), stem_05_4m_xyz_binID_stemID(:,5)); %%%%% 'groups' is a cell_array, each cell contains the points of a stem (x,y,z,bin_id)

    %%%%% Get the bin ids of each stem.
    groups2 = splitapply( @(x){x}, stem_05_4m_xyz_binID_stemID(:,4), stem_05_4m_xyz_binID_stemID(:,5)); 
    stem_bin_composition = cellfun(@unique,groups2,'Uniform', 0);

    %%%%%% Get the id of the stems which have points only in the top-most bin(0.5-1m) and the lowest bin (3.5-4m). They meight not be real stems %%%%%
    head_feet_stem_id=[];
    for i=1:length(stem_bin_composition)

        if length(stem_bin_composition{i})==1
            disp('error')            
        end
            
        if (length(stem_bin_composition{i})==2 & diff(stem_bin_composition{i})>=3) | (length(stem_bin_composition{i})==1)           
            head_feet_stem_id=[head_feet_stem_id;i];
        end

    end
    clear groups2 stem_bin_composition
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    stem_higher_bin_level_nice_id=[];  %%% to store the output
    stem_higher_bin_level_cylinder_bad_id=[]; %%% to store the output    
    stem_higher_bin_level_cylinder_good_id=[]; %%% to store the output
    
    giant_stems_id=[]; %%% to store the output    

    stem_lower_bin_level_good_id=[];
    stem_lower_bin_level_bad_id=[]; %%% to store the output
    
    stem_higher_bin_level_cylinder_good_afercylinder=cell(length(groups),1);
    counter2=0;
        
    for i2=1:length(groups) 
        
        if mod(i2,500)==0
            disp(strcat(num2str(100*i2/length(groups)),'%.....'))
        end
        
        if ismember(i2,head_feet_stem_id) %%% two components stems
            continue
        end
        
        %%%%%%%%%%%%%%%%%%%%%%% 90% time is consumed here, indexing large dataset is super time-consuming %%%%%%%%%%%%%%%%%%%%%%
%         one_trees_xyz_binID_stemID=stem_05_4m_xyz_binID_stemID(stem_05_4m_xyz_binID_stemID(:,5)==stemid_unique_temp(i2),:); %%% x y z bin_id stem_id
        
        %%%% much faster version %%%%
        one_trees_xyz_binID_stemID=groups{i2,1};
        one_trees_xyz_binID_stemID=[one_trees_xyz_binID_stemID ones(size(one_trees_xyz_binID_stemID,1),1)*i2]; %%% add a colum of tree ID
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if range(one_trees_xyz_binID_stemID(:,3))< height_threhold_05_15m %%%%% make sure again no short stems
            continue
        end
        
        %% For stems detected from higher bins
        if range(one_trees_xyz_binID_stemID(:,3))>=height_threhold_15_4m  %%% 
                        
            [dbh_max_temp,dbh_mean_temp]=xyrange_dbh_bins_10cm(one_trees_xyz_binID_stemID(:,1:3));                             

            if max(dbh_max_temp)>1 %%%%%%% giant trees with large DBH. Better check them visually. Will write them into a txt file.          
                giant_stems_id=[giant_stems_id;i2];                
                continue;               
            end
            
            %%%%%% Define a quality indicator for the stems: calcualte a diameter every 10cm, check the difference of the diameters                     
            if median(dbh_max_temp)<0.15 %%%%%%% pass some small trees (max bin DBH<15, or 20cm,or 0). Din't influence the results.
                filter_strenghth=0;   
            else            
                filter_strenghth=max(dbh_max_temp)/min(dbh_max_temp); %%% quality indicator
            end

            if filter_strenghth>1.5 %%%%% A threshold in quality indicator
                
                referenceVector = [0,0,1];
                maxAngularDistance=10;
                
                %%%%% pcfitcylinder gives slightly different result each time
                [model,inlierIndices,~,~] = pcfitcylinder(pointCloud(one_trees_xyz_binID_stemID(:,1:3)),0.05,referenceVector,maxAngularDistance);

                if range(one_trees_xyz_binID_stemID(inlierIndices,3))< min(1,range(one_trees_xyz_binID_stemID(:,3))*2/3) %model.Radius>1.5*median(dbh_max_temp)                 
                    
                    stem_higher_bin_level_cylinder_bad_id=[stem_higher_bin_level_cylinder_bad_id;i2];

                else                    

                    stem_higher_bin_level_cylinder_good_id=[stem_higher_bin_level_cylinder_good_id;i2];
                    
                    %%% apply rangesearch (10cm) on the inlierIndices points, to get more stem points                   
                    NS = createns(one_trees_xyz_binID_stemID(:,1:3),'NSMethod','exhaustive');
                    idx = rangesearch(NS,one_trees_xyz_binID_stemID(inlierIndices,1:3),0.1);
                    trunk_id=[];
                    for i=1:size(idx,1)
                        trunk_id=[trunk_id idx{i,1}];
                        trunk_id=unique(trunk_id);
                    end

%                     stem_higher_bin_level_cylinder_good_afercylinder=[stem_higher_bin_level_cylinder_good_afercylinder;one_trees_xyz_binID_stemID(trunk_id,:)];           
                    counter2=counter2+1;
                    stem_higher_bin_level_cylinder_good_afercylinder{counter2}=one_trees_xyz_binID_stemID(trunk_id,:);

                end

            else
                
                stem_higher_bin_level_nice_id=[stem_higher_bin_level_nice_id;i2]; %% High quality stems
                
            end
            
        
        else 
            %% %%% for stems detected from lower bins (0.5-1.5m)
            [dbh_max_temp,~]=xyrange_dbh_bins_10cm(one_trees_xyz_binID_stemID(:,1:3));
            filter_strenghth=sum(dbh_max_temp>0.15); %%%%% count how many bins have dbh>10cm, the numbers will be a strenght of good or bad. Try relative strength also
            filter_strenghth2=sum(dbh_max_temp>0.25); %%% 0.3
            
            if filter_strenghth>=3 | filter_strenghth2>=1 %| max(dbh_max_temp)/min(dbh_max_temp)>1.5               
                stem_lower_bin_level_bad_id=[stem_lower_bin_level_bad_id;i2];                
            else                
                stem_lower_bin_level_good_id=[stem_lower_bin_level_good_id;i2];
            end
            
        end

    end

    head_feet_stem=stem_05_4m_xyz_binID_stemID(ismember(stem_05_4m_xyz_binID_stemID(:,5),head_feet_stem_id),:);
        
    giant_stems=stem_05_4m_xyz_binID_stemID(ismember(stem_05_4m_xyz_binID_stemID(:,5),giant_stems_id),:);    
    stem_lower_bin_level_good=stem_05_4m_xyz_binID_stemID(ismember(stem_05_4m_xyz_binID_stemID(:,5),stem_lower_bin_level_good_id),:);
    stem_lower_bin_level_bad=stem_05_4m_xyz_binID_stemID(ismember(stem_05_4m_xyz_binID_stemID(:,5),stem_lower_bin_level_bad_id),:);   
    stem_higher_bin_level_cylinder_good=stem_05_4m_xyz_binID_stemID(ismember(stem_05_4m_xyz_binID_stemID(:,5),stem_higher_bin_level_cylinder_good_id),:);
    stem_higher_bin_level_cylinder_bad=stem_05_4m_xyz_binID_stemID(ismember(stem_05_4m_xyz_binID_stemID(:,5),stem_higher_bin_level_cylinder_bad_id),:);    
    stem_higher_bin_level_nice=stem_05_4m_xyz_binID_stemID(ismember(stem_05_4m_xyz_binID_stemID(:,5),stem_higher_bin_level_nice_id),:);    


    stem_higher_bin_level_cylinder_good_afercylinder(counter2+1:end)=[];
    stem_higher_bin_level_cylinder_good_afercylinder=cell2mat(stem_higher_bin_level_cylinder_good_afercylinder);
    
    all_good_stem=[giant_stems;stem_higher_bin_level_nice;stem_higher_bin_level_cylinder_good_afercylinder;stem_lower_bin_level_good];
    all_potential_stem=[head_feet_stem;stem_higher_bin_level_cylinder_bad;stem_lower_bin_level_bad];

    
    
%     fileID = fopen('./TLidar_output/all_cylinder_good_insidecylinder.txt','w');
%     fprintf(fileID,'%.3f,%.3f,%.3f,%i,%i\n',stem_higher_bin_level_cylinder_good_afercylinder(:,1:5)');
%     fclose(fileID); 
%     
%     fileID = fopen('./TLidar_output/all_cylinder_good_nocylinder.txt','w');
%     fprintf(fileID,'%.3f,%.3f,%.3f,%i,%i\n',stem_higher_bin_level_cylinder_good(:,1:5)');
%     fclose(fileID); 
    

    %% %%%%%%% DBH  %%%%%%%    
    
    labels_stem_potential=all_potential_stem(:,5);
    %%% relabel the stem ID from 1 to N consectively.    
    [~,~,labels_final_new]=unique(labels_stem_potential); clear labels_stem_potential
    all_potential_stem(:,5)=labels_final_new; clear labels_final_new
    
     
    labels_stem_good=all_good_stem(:,5);
    %%% relabel the stem ID from 1 to N consectively.    
    [~,~,labels_final_new]=unique(labels_stem_good); clear labels_stem_good
    all_good_stem(:,5)=labels_final_new; clear labels_final_new
    

    
    disp('DBH estimation.....')
    groups = splitapply( @(x){x}, all_good_stem(:,1:4), all_good_stem(:,5)); %%%%% 'groups' is a cell_array, each cell contains the points of a stem (x,y,z,bin_id)

    
    %%% DBH estimation
    groups_withDBH=DBH_TLS(groups);

    %% Write results
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% High quality stems %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fileID = fopen(strcat('./TLidar_output/all_stems_DBSCAN_corrected_furtherRefined.txt'),'a');
    
    str{1}='x';
    str{2}='y';
    str{3}='z';
    str{4}='bin_id';
    str{5}='stem_id';
    str{6}='DBH';

    fprintf(fileID,'%s,',str{1:5});
    fprintf(fileID,'%s\n',str{6});

    all_good_stem_out=cell2mat(groups_withDBH); 
    fprintf(fileID,'%.3f,%.3f,%.3f,%i,%i,%.3f\n',all_good_stem_out(:,1:6)'); %%% x, y, z, bin_id, stem_id, DBH
    fclose(fileID);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Potential stems %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fileID = fopen(strcat('./TLidar_output/all_potential_stems_final.txt'),'a');   
    str{1}='x';
    str{2}='y';
    str{3}='z';
    str{4}='bin_id';
    str{5}='stem_id';

    fprintf(fileID,'%s,',str{1:4});
    fprintf(fileID,'%s\n',str{5});

    fprintf(fileID,'%.3f,%.3f,%.3f,%i,%i\n',all_potential_stem(:,1:5)'); %%% x, y, z, bin_id, stem_id
    fclose(fileID);
	
    toc
    
    
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% visually check the stems with a DBH > 50cm to find and correct potential over-segmentation errors %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
end


%% delete temp file
% delete ./TLidar_output/bin_statistics
% delete ./TLidar_output/Petit_plateau_2cmsubsample_dist10cmseg_05_4m
% delete ./TLidar_output/Petit_plateau_2cmsubsample_dist10cmseg_05_4m_DBSCAN
