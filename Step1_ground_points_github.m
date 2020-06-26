%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Ground point extraction from TLS data.
%%% The input is a txt file, with three columns (x, y, z) seperated by comma. Use lasground to get the input file.
%%% Shengli Tao. Juin. 2020. Toulouse.
%%% See Tao et al. 2020. Methods in Ecology and Evolution (submitted) for details.

clear
clc

%% Parameters  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% %%% Z value of any point won't be 'z_diff' higher than its neighbors within  'xy_range' meter in the xy dimension %%%%%%%%%

xy_range=0.5; % the neighborhood search range. Unit: m
z_diff=0.5; % Unit: m. The larger this value, the more local variation in the Z values of the ground points. Set to ~0.5m to capture the woody debris on the ground.

%% Input txt file. Obtained from 'lasground'. Three columns without header: x, y and z  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gound_noise=dlmread('E:\TLidar\2019\lowest_5cm_single_scans_merged_lowest_5cm_cutplot_tiles\peitiplateau_lowest5cm_lasground_tile_reversed_groundonly.txt',',');
gound_noise=dlmread('./peitiplateau_lowest5cm_lasground_tile_reversed_groundonly.txt',',');


%%
%%%% Time consumption: 4335 seconds for processing 12ha data (33535463 points). Less than 10 minutes for 1ha.

%% Get the ground points.
tic

gound_noise_xyz=gound_noise(:,1:3);
gound_noise_xy=gound_noise(:,1:2);

original_pts_id=1:length(gound_noise_xy);

minx=min(gound_noise_xy(:,1));
maxx=max(gound_noise_xy(:,1));
miny=min(gound_noise_xy(:,2));
maxy=max(gound_noise_xy(:,2));

%%%% tile into ~2*2m  tiles, with a buffer of xyrange (50cm). Make sure the input file is larger than 2*2m. 

nx = round(min(range(gound_noise(:,1))/2,range(gound_noise(:,2))));
ny = nx;
xedge = linspace(minx,maxx,nx+1); % divide the plot into tiles
yedge = linspace(miny,maxy,ny+1); % divide the plot into tiles

xbin_allpts = discretize(gound_noise_xy(:,1), xedge);  % xbin index for all the pts
ybin_allpts = discretize(gound_noise_xy(:,2), yedge);  % ybin index for all the pts


%%%% Because each point will search its neighbors within 'xyrange' distance, we need to generate a buffer of size 'xyrange' around each tile %%%%
buffer_size = xy_range; 

xbin_allpts_left = discretize(gound_noise_xy(:,1), xedge-buffer_size);  % xbin index for all the pts
xbin_allpts_right = discretize(gound_noise_xy(:,1), xedge+buffer_size);  % xbin index for all the pts

ybin_allpts_left = discretize(gound_noise_xy(:,2), yedge-buffer_size);  % xbin index for all the pts
ybin_allpts_right = discretize(gound_noise_xy(:,2), yedge+buffer_size);  % xbin index for all the pts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ground_ind_final=nan(length(gound_noise_xy),1);

for i=1:length(xedge)-1
    
    disp(100*i/length(length(xedge)-1))

    for j=1:length(yedge)-1
        
%         disp(j)
        
        %%%%%%%%%%%%%%%% slow version %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         tile_withbuffer_xind=gound_noise_xy(:,1)>=xedge(i)-0.51 & gound_noise_xy(:,1)<xedge(i+1)+0.51; %%% 3136s out of 10800s
%         tile_withbuffer_yind=gound_noise_xy(:,2)>=yedge(j)-0.51 & gound_noise_xy(:,2)<yedge(j+1)+0.51; %%% 3136s out of 10800s       
%         pts_withintile_buffer_ind = tile_withbuffer_xind & tile_withbuffer_yind;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pts_withintile_buffer_ind = (xbin_allpts==i | xbin_allpts_left==i | xbin_allpts_right==i) & (ybin_allpts==j | ybin_allpts_left==j | ybin_allpts_right==j);
        
        
        pts_withintile_ind = xbin_allpts==i & ybin_allpts ==j;
        
        if sum(pts_withintile_ind)==0
            continue
        end
        
        %%%%%% pts in the tiles should also be in the tiles with buffer %%%%%%    %%% 1000s out of 10800s
%         if sum(ismember(original_pts_id(pts_withintile_ind),original_pts_id(pts_withintile_buffer_ind)))~= sum(pts_withintile_ind)
%             disp('error match....')
%             disp(i)
%             disp(j)
%         end
        
        
        pts_xyz1=gound_noise_xyz(pts_withintile_buffer_ind,:);

        ind_ground_tile_buffer=filter_ground(pts_xyz1,xy_range,z_diff);   %%% 1500s out of 10800s
        
        ind_ground_tile=ind_ground_tile_buffer(ismember(original_pts_id(pts_withintile_buffer_ind),original_pts_id(pts_withintile_ind)));
        
        Ground_ind_final(original_pts_id(pts_withintile_ind))=ind_ground_tile;
                
    end 
    
end

toc


%% %%%%%%%%%%%%%%%%%%%%%%%%% write out ground pts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outtxt_name=strcat('ground_point_range_',numstr(xy_range),'_z_',num2str(z_diff),'.txt');
fileID = fopen(outtxt_name,'w');
outdata2=gound_noise_xyz(logical(Ground_ind_final),:);
fprintf(fileID,'%.3f,%.3f,%.3f\n',outdata2');
fclose(fileID);
clear outdata2



