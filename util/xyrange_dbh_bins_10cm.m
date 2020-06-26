function [dbh_max,dbh_mean]=xyrange_dbh_bins_10cm(pts)  

    %%% calculate a diameter every 10cm along the z-axis

    %%% dbh_max: diamter calcaulted as maximum x- or y- range 
    %%% dbh_mean: diamter calcaulted as mean x- or y- range 
        

    binrange=(min(pts(:,3))-0.01):0.1:(max(pts(:,3))+0.01);
    
    if max(binrange)< max(pts(:,3))
        binrange=[binrange max(pts(:,3))+0.01];
    end
    [~,ind22]=histc(pts(:,3), binrange); 


    if ismember(0,ind22)
        disp('bin error')
    end

    dbh_mean=[];
    dbh_max=[];

    unique_ind22=unique(ind22);
    for jjj=1:length(unique_ind22)

        dbh_pts=pts(ind22==unique_ind22(jjj),:);

        if isempty(dbh_pts)
            continue
        end

        dbh_mean=[dbh_mean;mean([max(dbh_pts(:,1))-min(dbh_pts(:,1)) max(dbh_pts(:,2))-min(dbh_pts(:,2))])];
        dbh_max=[dbh_max;max([max(dbh_pts(:,1))-min(dbh_pts(:,1)) max(dbh_pts(:,2))-min(dbh_pts(:,2))])];

    end


    dbh_mean(dbh_mean==0)=[];
    dbh_max(dbh_max==0)=[];

        
end
