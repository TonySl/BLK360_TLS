%% Ground point extraction from TLS data. Version 0.0.  

%%% Shengli Tao. Juin. 2020. Toulouse. Tested in Matlab 2018b
%%% See Tao et al. xxxxx for details.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    !!!!!!!!!!!!!!!!!!   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Input files format   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    !!!!!!!!!!!!!!!!!!   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% The input is a txt or csv file, with at least three columns (x, y, z) seperated by comma or space, with or without header line 
%%% The first three colums should be x, y, z, in order.
%%% It's better if the input file contains the lowest points in every 5cm grid (or 2cm, 10cm). Finer grid size means higher computation burden.
%%% Please make sure your data cover an area larger than 2 x 2 m.



%%%% Time consumption: 4335 seconds for processing 12-ha data (33,535,462 points). Less than 10 minutes for 1ha data.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     Please specify inputfile and its full path name       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FullFileName='./figshare/test_region/lowest_5cm_testregion.txt'; %'E:\TLidar\2019\lowest_5cm_single_scans_merged_lowest_5cm_cutplot_tiles\peitiplateau_lowest5cm_lasground_tile_reversed_groundonly.txt'




%%  Default Parameters   
%%%%%%%%%%%%%%% %%% Z value of any point won't be 'z_diff' higher than its neighbors within  'xy_range' meter in the xy dimension %%%%%%%%%

xy_range=0.5; % the neighborhood search range. Unit: m. Each point will search its neighbors within xy_range meter horizontally. 
%%% 0.5 was set as default. Values between 0.1 and 0.5 all work. Larger value will increase computation burden.  

z_diff=0.5; % Unit: m. The larger this value, the more local variation in the Z values of the ground points. 
%%% Set to ~0.5m to capture the woody debris on the ground.




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ready to start data processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Ground point extraction

clc

addpath('./util')


tic

gound_noise_struct = importdata(FullFileName);

if isa(gound_noise_struct,'struct')
	gound_noise = gound_noise_struct.data;
    user_header_line = gound_noise_struct.colheaders;
else
    gound_noise = gound_noise_struct;
    user_header_line = [];
end



original_pts_id=1:length(gound_noise);

minx=min(gound_noise(:,1));
maxx=max(gound_noise(:,1));
miny=min(gound_noise(:,2));
maxy=max(gound_noise(:,2));

%%%% tile into ~2*2m  tiles, with a buffer of xyrange (50cm). 

nx = round(min(range(gound_noise(:,1))/2,range(gound_noise(:,2))));
ny = nx;
xedge = linspace(minx,maxx,nx+1); % divide the plot into tiles
yedge = linspace(miny,maxy,ny+1); % divide the plot into tiles

xbin_allpts = discretize(gound_noise(:,1), xedge);  % xbin index for all the pts
ybin_allpts = discretize(gound_noise(:,2), yedge);  % ybin index for all the pts


%%%% Because each point will search its neighbors within 'xyrange' distance, we need to generate a buffer of size 'xyrange' around each tile %%%%
buffer_size = xy_range; 

xbin_allpts_left = discretize(gound_noise(:,1), xedge-buffer_size);  % xbin index for all the pts
xbin_allpts_right = discretize(gound_noise(:,1), xedge+buffer_size);  % xbin index for all the pts

ybin_allpts_left = discretize(gound_noise(:,2), yedge-buffer_size);  % xbin index for all the pts
ybin_allpts_right = discretize(gound_noise(:,2), yedge+buffer_size);  % xbin index for all the pts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




Ground_ind_final=nan(length(gound_noise),1);

for i=1:length(xedge)-1
    
    disp(strcat(num2str(floor(100*i/(length(xedge)-1))),'%'))

    for j=1:length(yedge)-1
        
%         disp(j)
        

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
        
        
        pts_xyz1=gound_noise(pts_withintile_buffer_ind,1:3);

        ind_ground_tile_buffer=filter_ground(pts_xyz1,xy_range,z_diff);   %%% 1500s out of 10800s
        
        ind_ground_tile=ind_ground_tile_buffer(ismember(original_pts_id(pts_withintile_buffer_ind),original_pts_id(pts_withintile_ind)));
        
        Ground_ind_final(original_pts_id(pts_withintile_ind))=ind_ground_tile;
                
    end 
    
end

toc


%% %%%%%%%%%%%%%%%%%%%%%%%%% write out ground pts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outdata2=gound_noise(logical(Ground_ind_final),:);
% outtxt_name=strcat('ground_points_range_',num2str(xy_range),'_z_',num2str(z_diff),'.txt');
outtxt_name='ground_points_z_diff_filter.txt';

disp('..............................................................................')
disp('..............................................................................')
disp(strcat('write ground into: ',outtxt_name))

number_colums_final=size(outdata2,2); %%%% the first three colums are x y z, the last three colums are bin_id, stem_id, DBH
fileID = fopen(outtxt_name,'w');
str{1}='x';
str{2}='y';
str{3}='z';
user_colum=4:number_colums_final;
for i=1:length(user_colum)

    if isempty(user_header_line)
        str{user_colum(i)}=strcat('user_colum',num2str(i)); %%%% if no colum header in the input txt, name it as "user_colum1", "user_colum2"... 
    else
        str{user_colum(i)}=user_header_line{user_colum(i)}; %%%% keep the user colum header 
    end
end


fprintf(fileID,'%s,',str{1:number_colums_final-1});
fprintf(fileID,'%s\n',str{number_colums_final});

%%%%% will write all data into float precision
tempstr=repmat('%.3f,',1,number_colums_final);

fprintf(fileID,strcat(tempstr(1:end-1),'\n'),outdata2');
fclose(fileID);
clear outdata2 gound_noise
 


