function [dbh_max,dbh_mean,bin_id]=xyrange_dbh_bins(pts)  


    %%% cut the stem points into 5 smaller vertical bins, and calculate a diameter for each bin

    %%% dbh_max: diamter calcaulted as maximum x- or y- range 
    %%% dbh_mean: diamter calcaulted as mean x- or y- range 
        

    [~,ind22]=histc(pts(:,3), linspace(min(pts(:,3))-0.01,max(pts(:,3))+0.01,6)); % 5 bins should be 6 here              

    if ismember(0,ind22)
        disp('bin error')
    end

    dbh_mean=[];
    dbh_max=[];

    bin_id=[];

    unique_ind22=unique(ind22);
    for jjj=1:length(unique_ind22)

        dbh_pts=pts(ind22==unique_ind22(jjj),:);

        if isempty(dbh_pts)
            continue
        end

        dbh_mean=[dbh_mean;mean([max(dbh_pts(:,1))-min(dbh_pts(:,1)) max(dbh_pts(:,2))-min(dbh_pts(:,2))])];
        dbh_max=[dbh_max;max([max(dbh_pts(:,1))-min(dbh_pts(:,1)) max(dbh_pts(:,2))-min(dbh_pts(:,2))])];

        bin_id=[bin_id;jjj*ones(size(dbh_pts,1),1)];

    end


    dbh_mean(dbh_mean==0)=[];
    dbh_max(dbh_max==0)=[];


end
