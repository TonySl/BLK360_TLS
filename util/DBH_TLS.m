%%
%%% DBH estimation for each stem 

%%
function group_withDBH=DBH_TLS(groups1)


%     DBH_circle=nan(length(groups),1);    
%     DBH_cylinder=nan(length(groups),1);

    for i2=1:length(groups1) 

%             if mod(i2,500)==0
%                 disp(strcat(num2str(100*i2/length(groups1)),'%.....'))
%             end

            %%%%%%%%%%%%%%%%%%%%%%% 90% time is consumed here, indexing large dataset is super time-consuming %%%%%%%%%%%%%%%%%%%%%%
            %%%% one_trees_xyz_binID_stemID=stem_05_4m_xyz_binID_stemID(stem_05_4m_xyz_binID_stemID(:,5)==stemid_unique_temp(i2),:);

            %%%% much faster %%%%
            one_trees_xyz_binID_stemID=groups1{i2,1};

            %%%%%%%%%%%% DBH xy range- bin 10cm %%%%%%%%%%%%
            [~,dbh_mean1]=xyrange_dbh_bins_10cm(one_trees_xyz_binID_stemID(:,1:3));

            groups1{i2,1}(:,5)=i2;  %%% add stem_id as fifth colum            
            groups1{i2,1}(:,6)=nanmedian(dbh_mean1);  %%% add dbh as sixth colum


            %%%%%%%%%%%% DBH circle fitting %%%%%%%%%%%%
        %         if size(one_trees_xyz_binID_stemID,1)<3 %%%% at least 3 pts for fitting a circle
        %             DBH_circle(i2)=NaN;
        %         else
        %             Par = CircleFitByTaubin(one_trees_xyz_binID_stemID(:,1:2));  % x y radius  
        %             DBH_circle(i2)=Par(3)*2;
        %         end        


            %%%%%%%%%%%% DBH cylinder %%%%%%%%%%%%                        
        %         referenceVector = [0,0,1];
        %         maxAngularDistance=10;
        %         [model,inlierIndices,~,~] = pcfitcylinder(pointCloud(one_trees_xyz_binID_stemID(:,1:3)),0.05,referenceVector,maxAngularDistance);
        %         DBH_cylinder(i2)=model.Radius*2;


    end
    
    
    group_withDBH=groups1;
    
end