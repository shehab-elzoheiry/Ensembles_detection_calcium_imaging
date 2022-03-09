% % %        'F:\Imaging\Imaging\165\NewEnsembles\Ensemble_analysis_165.mat' ...
% % %        'F:\Imaging\Imaging\165\Cleaned_traces_165.mat'...
files={'F:\Imaging\155\NewEnsembles\Ensemble_analysis_155.mat' ...
       'F:\Imaging\155\Cleaned_traces_155.mat'...
       'F:\Imaging\157\NewEnsembles\Ensemble_analysis_157.mat' ...
       'F:\Imaging\157\Cleaned_traces_157.mat'...
       'F:\Imaging\165\NewEnsembles\Ensemble_analysis_165.mat' ...
       'F:\Imaging\165\Cleaned_traces_165.mat'...
       'F:\Imaging\172\NewEnsembles\Ensemble_analysis_172.mat' ...
       'F:\Imaging\172\Cleaned_traces_172.mat'...
       'F:\Imaging\174\NewEnsembles\Ensemble_analysis_174.mat' ...
       'F:\Imaging\174\Cleaned_traces_174.mat'...
       'F:\Imaging\189\NewEnsembles\Ensemble_analysis_189.mat' ...
       'F:\Imaging\189\Cleaned_traces_189.mat'...
       'F:\Imaging\191\NewEnsembles\Ensemble_analysis_191.mat' ...
       'F:\Imaging\191\Cleaned_traces_191.mat'...
       'F:\Imaging\194\NewEnsembles\Ensemble_analysis_194.mat' ...
       'F:\Imaging\194\Cleaned_traces_194.mat'...
};
% % 
% % % files={'/media/agkann/Elements/Imaging/155/NewEnsembles/Ensemble_analysis_155.mat' ...
% % %        '/media/agkann/Elements/Imaging/155/Cleaned_traces_155.mat'...
% % %        '/media/agkann/Elements/Imaging/157/NewEnsembles/Ensemble_analysis_157.mat' ...
% % %        '/media/agkann/Elements/Imaging/157/Cleaned_traces_157.mat'...
% % %        '/media/agkann/Elements/Imaging/165/NewEnsembles/Ensemble_analysis_165.mat' ...
% % %        '/media/agkann/Elements/Imaging/165/Cleaned_traces_165.mat'...
% % %        '/media/agkann/Elements/Imaging/172/NewEnsembles/Ensemble_analysis_172.mat' ...
% % %        '/media/agkann/Elements/Imaging/172/Cleaned_traces_172.mat'...
% % %        '/media/agkann/Elements/Imaging/174/NewEnsembles/Ensemble_analysis_174.mat' ...
% % %        '/media/agkann/Elements/Imaging/174/Cleaned_traces_174.mat'...
% % %        '/media/agkann/Elements/Imaging/189/NewEnsembles/Ensemble_analysis_189.mat' ...
% % %        '/media/agkann/Elements/Imaging/189/Cleaned_traces_189.mat'...
% % %        '/media/agkann/Elements/Imaging/191/NewEnsembles/Ensemble_analysis_191.mat' ...
% % %        '/media/agkann/Elements/Imaging/191/Cleaned_traces_191.mat'...
% % %        '/media/agkann/Elements/Imaging/194/NewEnsembles/Ensemble_analysis_194.mat' ...
% % %        '/media/agkann/Elements/Imaging/194/Cleaned_traces_194.mat'...
% % % };
% % 
% 
% 
% 
for zew=[1 3 5 7 9 11 13 15]
    tic
fortitle=[155 157 165 172 174 189 191 194];
she=[1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8];
    
    load(files{zew},'EnsRecActIdPlotSt','EnsRecActStFrames','sSt','offsetNormdSt');
    load(files{zew+1},'mask_');
% load('F:\Imaging\155\NewEnsembles\Ensemble_analysis_155.mat','sSt','offsetNormdSt','EnsRecActIdPlotSt','EnsRecActStFrames')
% load('F:\Imaging\155\Cleaned_traces_155.mat','mask_')
mask_(:,:,end)=[];


xx=EnsRecActIdPlotSt;
xx(~isnan(EnsRecActIdPlotSt))=1;
xx(isnan(xx))=0;
NormdSt=xx;                                             % marked cells at times of interesting events

% add=zeros(150-size(NormdSt,1),size(NormdSt,2));   NormdSt=vertcat(NormdSt,add);
% NormdSt=horzcat(NormdSt,NormdSt);

[coeff,score,latent,~,explained]=pca(NormdSt');         % PCA of the events
figure,biplot(coeff(:,1:3),'scores',score(:,1:3));      % plotting of PCA
savefig(['Fig1 PCA_' num2str(fortitle(she(zew))) '.fig']);
close all

    D = pdist(score(:,1:3));                            % calculates eucledian distance (default) for the score 
    Z = squareform(D);                                  % making the matrix a square because this is how the 'dunns' function (lool down) work
    %no=[2,3,4,5,6,7,8];                                % possibilities for number of clusters
    no= 2: sqrt(size(EnsRecActStFrames,2))+1;           % possibilities for number of clusters based on how big the sample is. The 1 is to add a cluster for the 0 events (where there are no ensembles detected)

    %options=[2, 15, 0.001, 1];
    %options=[2, 10, 0.01, 1];
    %options=[1.5, 15, 0.01, 1];                        % fuzzines (default 2) , no of iterations, min improvement btwn iteration , show flags or not 
    options=[2, 100,  0.0001, 1]; 
    universal={};
    for i=1:length(no);
%         score(:,1)=score(:,1)/max(score(:,1));
%         score(:,2)=score(:,2)/max(score(:,2));
%         score(:,3)=score(:,3)/max(score(:,3));
        [centers,U] = fcm(score(:,1:4),no(i),options);  % does the fuzzy clustering, calculates the centers and the coordinates for the observations
        maxU = max(U);
        if ~any(~isnan(U(:))), break, flag=1; end
        for hh=1:size(U,2)                                          
            ind(hh)= find(U(:,hh)== max(U(:,hh))  );    % determines the the membership of each frame to its cluster (to which it has highest membership 'U')
        end
        test(i)=length(unique(ind));
        if length(unique(ind))~=no(i)                   % sometimes although you ask for certain number of clusters, fcm fails to give that number, it gives something less,
            break                                       % therefore break the loop here, because it already means that further clusters won't be better
        end
        %DI=dunns(no(i),D,ind);
        DI(i)=dunns(no(i),Z,ind);                       % determines the DI to the number of clusters set for this iteration in the for-loop
        universal{i,1}=centers;
        universal{i,2}=U;
        clear centers U
    end
    if flag == 1, continue, end

    
    [~, bestNOclust]=find(DI==max(DI));
    centers=universal{bestNOclust,1};
    U=      universal{bestNOclust,2};
    %[centers,U] = fcm(score(:,1:3),no(bestNOclust),options);   maxU = max(U);  % do the clustering again based on the 'best number of clusters' to be used
    maxU = max(U);
    
    % figure for plotting the clustered data
    figure , hold on, grid on
    if     bestNOclust==1                               % actually 2 clusters
            index1 = find(U(1,:) == maxU); 
            index2 = find(U(2,:) == maxU);
            plot3(score(index1,1),score(index1,2),score(index1,3),'ob')
            plot3(score(index2,1),score(index2,2),score(index2,3),'or')
            plot3(centers(1,1),centers(1,2),centers(1,3),'xb','MarkerSize',10,'LineWidth',3)
            plot3(centers(2,1),centers(2,2),centers(2,3),'xr','MarkerSize',10,'LineWidth',3)
            if isempty(index1) || isempty(index2); disp('One cluster is empty'); end
    elseif bestNOclust==2                               % actually 3 clusters
            index1 = find(U(1,:) == maxU);
            index2 = find(U(2,:) == maxU);
            index3 = find(U(3,:) == maxU);
            plot3(score(index1,1),score(index1,2),score(index1,3),'ob')
            plot3(score(index2,1),score(index2,2),score(index2,3),'or')
            plot3(score(index3,1),score(index3,2),score(index3,3),'ok')
            plot3(centers(1,1),centers(1,2),centers(1,3),'xb','MarkerSize',10,'LineWidth',3)
            plot3(centers(2,1),centers(2,2),centers(2,3),'xr','MarkerSize',10,'LineWidth',3)
            plot3(centers(3,1),centers(3,2),centers(3,3),'xk','MarkerSize',10,'LineWidth',3)
            if isempty(index1) || isempty(index2) || isempty(index3); disp('One cluster is empty'); end
    elseif bestNOclust==3                               % actually 4 clusters
            index1 = find(U(1,:) == maxU);
            index2 = find(U(2,:) == maxU);
            index3 = find(U(3,:) == maxU);
            index4 = find(U(4,:) == maxU);
            plot3(score(index1,1),score(index1,2),score(index1,3),'ob')
            plot3(score(index2,1),score(index2,2),score(index2,3),'or')
            plot3(score(index3,1),score(index3,2),score(index3,3),'ok')
            plot3(score(index4,1),score(index4,2),score(index4,3),'og')
            plot3(centers(1,1),centers(1,2),centers(1,3),'xb','MarkerSize',10,'LineWidth',3)
            plot3(centers(2,1),centers(2,2),centers(2,3),'xr','MarkerSize',10,'LineWidth',3)
            plot3(centers(3,1),centers(3,2),centers(3,3),'xk','MarkerSize',10,'LineWidth',3)
            plot3(centers(4,1),centers(4,2),centers(4,3),'xg','MarkerSize',10,'LineWidth',3)
            if isempty(index1) || isempty(index2) || isempty(index3) || isempty(index4); disp('One cluster is empty'); end
    elseif bestNOclust==4                               % actually 5 clusters
            index1 = find(U(1,:) == maxU);
            index2 = find(U(2,:) == maxU);
            index3 = find(U(3,:) == maxU);
            index4 = find(U(4,:) == maxU);
            index5 = find(U(5,:) == maxU);
            plot3(score(index1,1),score(index1,2),score(index1,3),'ob')
            plot3(score(index2,1),score(index2,2),score(index2,3),'or')
            plot3(score(index3,1),score(index3,2),score(index3,3),'ok')
            plot3(score(index4,1),score(index4,2),score(index4,3),'og')
            plot3(score(index5,1),score(index5,2),score(index5,3),'oy')
            plot3(centers(1,1),centers(1,2),centers(1,3),'xb','MarkerSize',10,'LineWidth',3)
            plot3(centers(2,1),centers(2,2),centers(2,3),'xr','MarkerSize',10,'LineWidth',3)
            plot3(centers(3,1),centers(3,2),centers(3,3),'xk','MarkerSize',10,'LineWidth',3)
            plot3(centers(4,1),centers(4,2),centers(4,3),'xg','MarkerSize',10,'LineWidth',3)
            plot3(centers(5,1),centers(5,2),centers(5,3),'xy','MarkerSize',10,'LineWidth',3)
            if isempty(index1) || isempty(index2) || isempty(index3) || isempty(index4) || isempty(index5); disp('One cluster is empty'); end
    elseif bestNOclust==5                               % actually 6 clusters
            index1 = find(U(1,:) == maxU);
            index2 = find(U(2,:) == maxU);
            index3 = find(U(3,:) == maxU);
            index4 = find(U(4,:) == maxU);
            index5 = find(U(5,:) == maxU);
            index6 = find(U(6,:) == maxU);
            plot3(score(index1,1),score(index1,2),score(index1,3),'ob')
            plot3(score(index2,1),score(index2,2),score(index2,3),'or')
            plot3(score(index3,1),score(index3,2),score(index3,3),'ok')
            plot3(score(index4,1),score(index4,2),score(index4,3),'og')
            plot3(score(index5,1),score(index5,2),score(index5,3),'oy')
            plot3(score(index6,1),score(index6,2),score(index6,3),'om')
            plot3(centers(1,1),centers(1,2),centers(1,3),'xb','MarkerSize',10,'LineWidth',3)
            plot3(centers(2,1),centers(2,2),centers(2,3),'xr','MarkerSize',10,'LineWidth',3)
            plot3(centers(3,1),centers(3,2),centers(3,3),'xk','MarkerSize',10,'LineWidth',3)
            plot3(centers(4,1),centers(4,2),centers(4,3),'xg','MarkerSize',10,'LineWidth',3)
            plot3(centers(5,1),centers(5,2),centers(5,3),'xy','MarkerSize',10,'LineWidth',3)
            plot3(centers(6,1),centers(6,2),centers(6,3),'xm','MarkerSize',10,'LineWidth',3)
            if isempty(index1) || isempty(index2) || isempty(index3) || isempty(index4) || isempty(index5) || isempty(index6); disp('One cluster is empty'); end
    elseif bestNOclust==6                               % actually 7 clusters
            index1 = find(U(1,:) == maxU);
            index2 = find(U(2,:) == maxU);
            index3 = find(U(3,:) == maxU);
            index4 = find(U(4,:) == maxU);
            index5 = find(U(5,:) == maxU);
            index6 = find(U(6,:) == maxU);
            index7 = find(U(7,:) == maxU);
            plot3(score(index1,1),score(index1,2),score(index1,3),'ob')
            plot3(score(index2,1),score(index2,2),score(index2,3),'or')
            plot3(score(index3,1),score(index3,2),score(index3,3),'ok')
            plot3(score(index4,1),score(index4,2),score(index4,3),'og')
            plot3(score(index5,1),score(index5,2),score(index5,3),'oy')
            plot3(score(index6,1),score(index6,2),score(index6,3),'om')
            plot3(score(index7,1),score(index7,2),score(index7,3),'r.')
            plot3(centers(1,1),centers(1,2),centers(1,3),'xb','MarkerSize',10,'LineWidth',3)
            plot3(centers(2,1),centers(2,2),centers(2,3),'xr','MarkerSize',10,'LineWidth',3)
            plot3(centers(3,1),centers(3,2),centers(3,3),'xk','MarkerSize',10,'LineWidth',3)
            plot3(centers(4,1),centers(4,2),centers(4,3),'xg','MarkerSize',10,'LineWidth',3)
            plot3(centers(5,1),centers(5,2),centers(5,3),'xy','MarkerSize',10,'LineWidth',3)
            plot3(centers(6,1),centers(6,2),centers(6,3),'xm','MarkerSize',10,'LineWidth',3)
            plot3(centers(7,1),centers(7,2),centers(7,3),'xr','MarkerSize',10,'LineWidth',3)
            if isempty(index1) || isempty(index2) || isempty(index3) || isempty(index4) || isempty(index5) || isempty(index6) || isempty(index7); disp('One cluster is empty'); end
    elseif bestNOclust==7                               % actually 8 clusters             
            index1 = find(U(1,:) == maxU);
            index2 = find(U(2,:) == maxU);
            index3 = find(U(3,:) == maxU);
            index4 = find(U(4,:) == maxU);
            index5 = find(U(5,:) == maxU);
            index6 = find(U(6,:) == maxU);
            index7 = find(U(7,:) == maxU);
            index8 = find(U(8,:) == maxU);            
            plot3(score(index1,1),score(index1,2),score(index1,3),'ob')
            plot3(score(index2,1),score(index2,2),score(index2,3),'or')
            plot3(score(index3,1),score(index3,2),score(index3,3),'ok')
            plot3(score(index4,1),score(index4,2),score(index4,3),'og')
            plot3(score(index5,1),score(index5,2),score(index5,3),'oy')
            plot3(score(index6,1),score(index6,2),score(index6,3),'om')
            plot3(score(index7,1),score(index7,2),score(index7,3),'r.')
            plot3(score(index8,1),score(index8,2),score(index8,3),'b.')
            plot3(centers(1,1),centers(1,2),centers(1,3),'xb','MarkerSize',10,'LineWidth',3)
            plot3(centers(2,1),centers(2,2),centers(2,3),'xr','MarkerSize',10,'LineWidth',3)
            plot3(centers(3,1),centers(3,2),centers(3,3),'xk','MarkerSize',10,'LineWidth',3)
            plot3(centers(4,1),centers(4,2),centers(4,3),'xg','MarkerSize',10,'LineWidth',3)
            plot3(centers(5,1),centers(5,2),centers(5,3),'xy','MarkerSize',10,'LineWidth',3)
            plot3(centers(6,1),centers(6,2),centers(6,3),'xm','MarkerSize',10,'LineWidth',3)
            plot3(centers(7,1),centers(7,2),centers(7,3),'xr','MarkerSize',10,'LineWidth',3)        
            plot3(centers(8,1),centers(8,2),centers(8,3),'xb','MarkerSize',10,'LineWidth',3)   
            if isempty(index1) || isempty(index2) || isempty(index3) || isempty(index4) || isempty(index5) || isempty(index6) || isempty(index7) || isempty(index8); disp('One cluster is empty'); end
    end
    xlabel('PCA 1'), ylabel ('PCA 2'), zlabel('PCA 3');
    title(num2str(fortitle(she(zew))));
    fig_info = get(gca);
    savefig(['Fig2 PCA_' num2str(fortitle(she(zew))) '_Ensembles.fig']); 
    

    
    
    summary={};
    summary{1}= 'Id of Events'; % how many frames (events)
    summary{2}= EnsRecActStFrames;  % id of frames
    
       
    if     no(bestNOclust)==2
%          summary{3}=ones(1,length(EnsRecActStFrames));
           all_ind=zeros(size(NormdSt,2),2);
           all_ind(1:length(index1),1)=index1;
           all_ind(1:length(index2),2)=index2;
           for k=1:length(EnsRecActStFrames), [~, clus_ID(k)]=find(all_ind==EnsRecActStFrames(k)); end
           leng_checker=[length(index1) length(index2)];
    elseif no(bestNOclust)==3
           all_ind=zeros(size(NormdSt,2),3);
           all_ind(1:length(index1),1)=index1;
           all_ind(1:length(index2),2)=index2;
           all_ind(1:length(index3),3)=index3;
           for k=1:length(EnsRecActStFrames), [~, clus_ID(k)]=find(all_ind==EnsRecActStFrames(k)); end
           leng_checker=[length(index1) length(index2) length(index3)];
    elseif no(bestNOclust)==4
           all_ind=zeros(size(NormdSt,2),4);
           all_ind(1:length(index1),1)=index1;
           all_ind(1:length(index2),2)=index2;
           all_ind(1:length(index3),3)=index3;           
           all_ind(1:length(index4),4)=index4;
           for k=1:length(EnsRecActStFrames), [~, clus_ID(k)]=find(all_ind==EnsRecActStFrames(k)); end
           leng_checker=[length(index1) length(index2) length(index3) length(index4)];
    elseif no(bestNOclust)==5
           all_ind=zeros(size(NormdSt,2),5);
           all_ind(1:length(index1),1)=index1;
           all_ind(1:length(index2),2)=index2;
           all_ind(1:length(index3),3)=index3;           
           all_ind(1:length(index4),4)=index4;           
           all_ind(1:length(index5),5)=index5;
           for k=1:length(EnsRecActStFrames), [~, clus_ID(k)]=find(all_ind==EnsRecActStFrames(k)); end
           leng_checker=[length(index1) length(index2) length(index3) length(index4) length(index5)];
    elseif no(bestNOclust)==6
           all_ind=zeros(size(NormdSt,2),6);
           all_ind(1:length(index1),1)=index1;
           all_ind(1:length(index2),2)=index2;
           all_ind(1:length(index3),3)=index3;           
           all_ind(1:length(index4),4)=index4;           
           all_ind(1:length(index5),5)=index5;
           all_ind(1:length(index6),6)=index6;
           for k=1:length(EnsRecActStFrames), [~, clus_ID(k)]=find(all_ind==EnsRecActStFrames(k)); end
           leng_checker=[length(index1) length(index2) length(index3) length(index4) length(index5) length(index6)];
    elseif no(bestNOclust)==7
           all_ind=zeros(size(NormdSt,2),7);
           all_ind(1:length(index1),1)=index1;
           all_ind(1:length(index2),2)=index2;
           all_ind(1:length(index3),3)=index3;           
           all_ind(1:length(index4),4)=index4;           
           all_ind(1:length(index5),5)=index5;
           all_ind(1:length(index6),6)=index6;
           all_ind(1:length(index7),7)=index7;
           for k=1:length(EnsRecActStFrames), [~, clus_ID(k)]=find(all_ind==EnsRecActStFrames(k)); end
           leng_checker=[length(index1) length(index2) length(index3) length(index4) length(index5) length(index6) length(index7)];    
    elseif no(bestNOclust)==8
           all_ind=zeros(size(NormdSt,2),8);
           all_ind(1:length(index1),1)=index1;
           all_ind(1:length(index2),2)=index2;
           all_ind(1:length(index3),3)=index3;           
           all_ind(1:length(index4),4)=index4;           
           all_ind(1:length(index5),5)=index5;
           all_ind(1:length(index6),6)=index6;
           all_ind(1:length(index7),7)=index7;
           all_ind(1:length(index8),8)=index8;
           for k=1:length(EnsRecActStFrames), [~, clus_ID(k)]=find(all_ind==EnsRecActStFrames(k)); end
           leng_checker=[length(index1) length(index2) length(index3) length(index4) length(index5) length(index6) length(index7) length(index8)];         
    end
summary{3}=clus_ID;

    for del=1:size(all_ind,2)
        if ~any(all_ind(:,del)) , %to_be_del(del)=del; 
            all_ind(:,del)=[];
        end
    end
   
    for z=1:size(all_ind,2)
        length_check(z:z)=length(find(all_ind(:,z))~=0);
    end
    del_aswell=find(length_check==max(length_check));
    %[~,del_aswell]=find(leng_checker==max(leng_checker));
    all_ind(:,del_aswell)=[];
    
    
% Identify the cells in each ensemble
for ens=1:size(all_ind,2)
    summary{ens,4}=find(sum(NormdSt(:,all_ind(1:max(find(all_ind(:,ens)~=0)),ens)),2)~=0); 
end               
                  

% plot traces with motifs labeled on them
whichfig=figure(3);
plot(sSt,offsetNormdSt,'color',[0.4 0.7 1])                             %%% plots normalised traces of all cells
hold on
for j=1:size(EnsRecActStFrames,2)
    color=fig_info.Children(  length(fig_info.Children) - summary{1,3}(j) +1 ).Color;
    plot(sSt(1,summary{1,2}(j)),(EnsRecActIdPlotSt(:,  summary{1,2}(j)  )),'.','Color',color,'MarkerSize',7)                       %%% marks (red dots) cells recruited in an ensemble
end
    title(num2str(fortitle(she(zew))));
    savefig(whichfig,['Fig3 Traces_' num2str(fortitle(she(zew))) '_Ensembles.fig']); close

  


% Find the center points of each ROI (used later for plotting the map and calculating mahalonobis dist)
for l=1:size(mask_,3)
measurements = regionprops(mask_(:,:,l), 'Centroid');
allCentroids(l,:) = [measurements.Centroid];
end

% plotting map of each ensemble with the correspondin color 
id_ens=unique(summary{1,3});
for no_ens=1:length(id_ens)
    figure, hold on
    color=fig_info.Children(  length(fig_info.Children) -id_ens(no_ens)+1 ).Color;
    incolor=color+0.8; incolor(incolor>1)=1;
    plot(allCentroids(summary{no_ens,4},1) ,allCentroids(summary{no_ens,4},2), 'o' ,'MarkerFaceColor',incolor,'Color', color )
    title(['Ensemble_' num2str(no_ens) '_' num2str(fortitle(she(zew)))]);
    savefig(['Fig4 Map_' num2str(no_ens) '_' num2str(fortitle(she(zew))) '_.fig']);
    close 
    
    % measure size of the ensemble
    meanX=mean(allCentroids(summary{no_ens,4},1));          % x-coordinate of the center point
    meanY=mean(allCentroids(summary{no_ens,4},2));          % y-coordinate of the center point
    ref(1:length(allCentroids(summary{no_ens,4})),1)=meanX; 
    ref(1:length(allCentroids(summary{no_ens,4})),2)=meanY;
    d2 = mahal(allCentroids(summary{no_ens,4},:),ref);      % d2 is squared Mahalanobis distance
    sizeofNetwork=mean(sqrt(d2));
    summary{1,5}= 'Network Size';                           % avg Mahalanobis distance
    summary{no_ens,6}=  sizeofNetwork;
end


% cells overlap summary
for p = 1:size(summary,1)
    x(p) = size(summary{p,4},1);  % size each ensemble (number of cells in each)
end

% ID of cells that participated in any ensemble activity
cells=nan(max(x),size(summary,1));
for m=1:size(summary,1)
    cells(1:length(summary{m,4}),m)=summary{m,4};
end

members(:,1)=unique(cells);                                     
members(isnan(members))=[];


for mm=1:size(members,1)
    members(mm,2)=length(find(cells==members(mm,1)));           % find how many times each cell participated in ensemble activity
end

summary{1,7}=  'cells_once';
summary{1,8}=  length(find(members(:,2)==1))/size(EnsRecActIdPlotSt,1)*100;
summary{1,9}=  'cells_several';
% summary{1,10}= length(find(members(:,2)>1 & members(:,2)< size(summary,1)))/size(EnsRecActIdPlotSt,1)*100;
summary{1,10}= length(find(members(:,2)>1 ))/size(EnsRecActIdPlotSt,1)*100;
summary{1,11}= 'cells_always';
summary{1,12}= length(find(members(:,2)==size(summary,1)))/size(EnsRecActIdPlotSt,1)*100;

summary{1,13}= 'no of cells';
summary{1,14}= (size(EnsRecActIdPlotSt,1));
summary{1,15}= 'no of events';
summary{1,16}= (length(summary{1,2}));
summary{1,17}= 'no of motifs';
summary{1,18}= (size(summary,1));

save(['Summary_' num2str(fortitle(she(zew))) '.mat'], 'summary');
save(['Variables_' num2str(fortitle(she(zew))) '.mat']);
clearvars -except files
close all
toc   

end



%% this part should be after the summing up of all summary files 
%correlating no of cells to no of events

summaries=dir('Summary*');
for n=1:length(summaries)
    load (summaries(n).name)
    noOFevents(n:n)=summary{1,16};
    noOFcells(n:n)=summary{1,14};
    noOFmotifs(n:n)=summary{1,18};
    
    for r=1:size(summary,1),x(r)= summary{r,6} ;end
    networksize(n:n)=mean(x);
    cells_once(n:n)=summary{1,8};
    cells_several(n:n)=summary{1,10};
    cells_always(n:n)=summary{1,12};
    clear summary
end



% correlations
        for i=1:8
            corrCells_Events(i,1)=noOFcells(i);
            corrCells_Events(i,2)=noOFevents(i);
        end
        [R_CE,P_CE] = corrcoef(corrCells_Events);   % _CE means Cells to Events
        figure, for i=1:8, plot(corrCells_Events(i,1) ,corrCells_Events(i,2),'kx','LineWidth',1.5,'MarkerSize',15); hold on; end
        xlabel('noOFcells'); ; ylabel('noOFevents');  xlim([0 max(noOFcells)+10]) ; ylim([0 max(noOFevents)+4]);
        title(['R = ' num2str(R_CE(1,2)) ', P = ' num2str(P_CE(1,2)) 'Cells_Events']);
        savefig(['Fig_Correlating_' 'Cells_vs_Events' '_.fig']); close

        % correlating no of cells to no of motifs
        for i=1:8
            corrCells_Motifs(i,1)=noOFcells(i);
            corrCells_Motifs(i,2)=noOFmotifs(i);
        end
        [R_CM,P_CM] = corrcoef(corrCells_Motifs);   % _CM means Cells to Motifs
        figure, for i=1:8, plot(corrCells_Motifs(i,1) ,corrCells_Motifs(i,2),'xk','LineWidth',1.5,'MarkerSize',15); hold on; end
        title(['R = ' num2str(R_CM(1,2)) ', P = ' num2str(P_CM(1,2)) 'Cells_Motifs']);
        xlabel('noOFcells'); ylabel('noOFmotifs');   xlim([0 max(noOFcells)+10]) ; ylim([0 max(noOFmotifs)+1]);
        savefig(['Fig_Correlating_' 'Cells_vs_Motifs' '_.fig']); close
        
        % correlating no of events to no of motifs
        for i=1:8
            corrEvents_Motifs(i,1)=noOFevents(i);
            corrEvents_Motifs(i,2)=noOFmotifs(i);
        end
        [R_EM,P_EM] = corrcoef(corrEvents_Motifs);   % _EM means Events to Motifs
        figure, for i=1:8, plot(corrEvents_Motifs(i,1) ,corrEvents_Motifs(i,2),'kx','LineWidth',1.5,'MarkerSize',15); hold on; end
        title(['R = ' num2str(R_EM(1,2)) ', P = ' num2str(P_EM(1,2)) 'Events_Motifs']);
        xlabel('noOFevents'); ylabel('noOFmotifs'); xlim([0 max(noOFevents)+4]) ; ylim([0 max(noOFmotifs)+1]);
        savefig(['Fig_Correlating_' 'Events_vs_Motifs' '_.fig']); close


%% The map of ensembles with their SD of activity 
    

    no=[155 157 165 172 174 189 191 194 ];   
    summaries=dir('Summary*');    
for zew=1:length(no)
    load(summaries(zew).name,'summary');
    %load(['F:\PhD\blowup ensembles\Final_Results_of_blowup\test\test'  '\Summary_' num2str(no(zew)) '.mat'],'summary');
    for i=1:size(summary,1)
        Ensembles{zew:zew,i:i}=summary{i,4}; 
    end
    clear summary
end
    gam_rot_ens=Ensembles;
    
    

no_gamma=[155     157     165     172     174     189     191     194    ];
no_roten=[    156     158     166     173     175     190     192     195];
for m=1:length(no_gamma)
    load(['F:\Imaging\' num2str(no_gamma(m)) '\Cleaned_traces_' num2str(no_gamma(m)) '.mat'],'clean_traces','mask_'); 
    load(['F:\Imaging\' num2str(no_gamma(m)) '\Cleaned_traces_' num2str(no_gamma(m)) '.mat'],'mask_'); 
    %load(['F:\Imaging\' num2str(no_gamma(m)) '\NewEnsembles\Ensemble_analysis_' num2str(no_gamma(m)) '.mat'],'NormdSt'); 
    preclean_traces_gam{m:m} = clean_traces'; %%% shehab
    masks{m:m} = mask_ ;
    clear clean_traces mask_
    load(['F:\Imaging\' num2str(no_roten(m)) '\Cleaned_traces_' num2str(no_roten(m)) '.mat'],'clean_traces'); 
    load(['F:\Imaging\' num2str(no_roten(m)) '\Cleaned_traces_' num2str(no_roten(m)) '.mat'],'mask_'); 
    %load(['F:\Imaging\' num2str(no_roten(m)) '\NewEnsembles\Ensemble_analysis_' num2str(no_roten(m)) '.mat'],'NormdSt'); 
    preclean_traces_rot{m:m} = clean_traces';  %%% shehab
    clear clean_traces
end


for i=1:size(gam_rot_ens,1)
    clean_traces_gam = preclean_traces_gam{i};
    clean_traces_rot = preclean_traces_rot{i};
    mask_= masks{i};
    for j=1:length(gam_rot_ens(i,:))
        for k=1:length(gam_rot_ens{i,j})
            gam_rot_ens{i,j}(k,2) = trapz(clean_traces_gam(gam_rot_ens{i,j}(k,1),:),2); % AUC for this cell during gamma 
            gam_rot_ens{i,j}(k,3) = trapz(clean_traces_rot(gam_rot_ens{i,j}(k,1),:),2); % AUC for this cell during rotenone 


            measurements = regionprops(mask_(:,:,gam_rot_ens{i,j}(k,1)), 'Centroid');
            allCentroids(1,:) = [measurements.Centroid];
            gam_rot_ens{i,j}(k,4:5)=allCentroids;
            gam_rot_ens{i,j}(1,6)=mean(gam_rot_ens{i,j}(:,2));
            gam_rot_ens{i,j}(1,7)=std(gam_rot_ens{i,j}(:,2));
            gam_rot_ens{i,j}(1,8)=std(gam_rot_ens{i,j}(:,3));
        end
    end
    clear clean_traces_gam clean_traces_rot
end
    
for i=1:size(gam_rot_ens,1)
    for j=1:size(gam_rot_ens,2)
        if isempty(gam_rot_ens{i,j}), break, end
        
        figure, xlim([0 672]); ylim([0 512]); 
        ratio=log10((gam_rot_ens{i,j}(:,3))./(gam_rot_ens{i,j}(:,2))); %%%%%%%%%%%%%% Maps are showing log the diff in ratio
        color=ratio;
        scatter(gam_rot_ens{i,j}(:,4),gam_rot_ens{i,j}(:,5),25,color,'filled') 
        colorbar
        
        % for all cells in all ensembles
            gam_ratio2avg=(gam_rot_ens{i,j}(:,2))./(gam_rot_ens{i,j}(1,6));
            rot_ratio2avg=(gam_rot_ens{i,j}(:,3))./(gam_rot_ens{i,j}(1,6));
            if i==1 && j==1
                    heatmap(1,:)=gam_ratio2avg';
                    heatmap(2,:)=rot_ratio2avg';
                    heatmap_log(1,:)=log10(gam_ratio2avg');
                    heatmap_log(2,:)=log10(rot_ratio2avg');
            else
                    ind=length(heatmap); % trick to avoid the changing size after the next line which screws the following line 
                    heatmap(1,ind+1:ind+length(gam_ratio2avg))=gam_ratio2avg';
                    heatmap(2,ind+1:ind+length(rot_ratio2avg))=rot_ratio2avg';
                    ind=length(heatmap_log); % trick to avoid the changing size after the next line which screws the following line 
                    heatmap_log(1,ind+1:ind+length(gam_ratio2avg))=log10(gam_ratio2avg');
                    heatmap_log(2,ind+1:ind+length(rot_ratio2avg))=log10(rot_ratio2avg');
            end

        % for ensembles as a whole
        gam2gamSD=1;
        rot2gamSD=(gam_rot_ens{i,j}(1,8))/(gam_rot_ens{i,j}(1,7));
        
        if   i==1 && j==1 , ens_SD=1; ens_SD(2,1)=rot2gamSD; 
        else
             ind=size(ens_SD,2);
             ens_SD(1,ind+1)=gam2gamSD;
             ens_SD(2,ind+1)=rot2gamSD;
        end
% % % % % % %         savefig(['Slice_' num2str(no_gamma(i)) 'Ensemble_' num2str(j) '_.fig']); close
        clear gam_ratio2avg rot_ratio2avg
    end
end
    
    % each cell compared to the average of its ensemble in gamma (Normal)
    figure
    imagesc(heatmap) , title('normal'); colorbar
    savefig('Ensembles_cellactivityTOaverageOfEnsemble_cartesian_.fig'); close
    
    % SORTED each cell compared to the average of its ensemble in gamma (Normal)
    newheatmap = sortrows(heatmap',1) ;
    figure
    subplot(3,1,2:3)
    imagesc(newheatmap') , title('sorted gamma normal'); colorbar
    subplot(3,1,1)
    histogram(newheatmap(:,1))
    savefig('Ensembles_SortGAMcellactivityTOaverageOfEnsemble_cartesian_.fig'); close
    
    newheatmap = sortrows(heatmap',2) ;
    figure
    subplot(3,1,1:2)
    imagesc(newheatmap') , title('sorted Rotenone normal'); colorbar
    subplot(3,1,3)
    histogram(newheatmap(:,2))
    savefig('Ensembles_SortROTcellactivityTOaverageOfEnsemble_cartesian_.fig'); close
    
    % each cell compared to the average of its ensemble in gamma (LOG)
    figure
    imagesc(heatmap_log), title('log'); colorbar
    savefig('Ensembles_cellactivityTOaverageOfEnsemble_log_.fig'); close

    % SORTED each cell compared to the average of its ensemble in gamma (LOG)
    newheatmap_log = sortrows(heatmap_log',1) ;
    figure
    subplot(3,1,2:3)
    imagesc(newheatmap_log'), title('sorted gamma log'); colorbar
    subplot(3,1,1)
    histogram(heatmap_log(1,:))
    savefig('Ensembles_SortGAMcellactivityTOaverageOfEnsemble_log_.fig'); close
    
    newheatmap_log = sortrows(heatmap_log',2) ;
    figure
    subplot(3,1,1:2)
    imagesc(newheatmap_log'), title('sorted Rotenone log'); colorbar
    subplot(3,1,3)
    histogram(heatmap_log(2,:))
    savefig('Ensembles_SortROTcellactivityTOaverageOfEnsemble_log_.fig'); close
    
    % all ensembles
    figure
    imagesc(ens_SD), title('Ensembles_SD_cartesian'); colorbar
    savefig('Ensembles_SD_log_.fig'); close
    
    figure
    imagesc(log10(ens_SD)), title('Ensembles_SD_log'); colorbar
    savefig('Ensembles_SD_log_.fig'); close

%% extract activity of all cells from 'gam_rot_ens'
% 
% for i=1:size(gam_rot_ens,1)
%     for j=1:size(gam_rot_ens,2)
%         if isempty(gam_rot_ens{i,j}), break, end
%         
%             gam_avg=gam_rot_ens{i,j}(:,2);
%             rot_avg=gam_rot_ens{i,j}(:,3);
%             if i==1 && j==1
%                     heatmap(1,:)=gam_avg';
%                     heatmap(2,:)=rot_avg';
%             else
%                     ind=length(heatmap); % trick to avoid the changing size after the next line which screws the following line 
%                     heatmap(1,ind+1:ind+length(gam_avg))=gam_avg';
%                     heatmap(2,ind+1:ind+length(rot_avg))=rot_avg';
%             end
% 
%        
%     end
% end

%% Correcting cells countings
files={'F:\Imaging\155\NewEnsembles\Ensemble_analysis_155.mat' ...
       'F:\Imaging\157\NewEnsembles\Ensemble_analysis_157.mat' ...
       'F:\Imaging\165\NewEnsembles\Ensemble_analysis_165.mat' ...
       'F:\Imaging\172\NewEnsembles\Ensemble_analysis_172.mat' ...
       'F:\Imaging\174\NewEnsembles\Ensemble_analysis_174.mat' ...
       'F:\Imaging\189\NewEnsembles\Ensemble_analysis_189.mat' ...
       'F:\Imaging\191\NewEnsembles\Ensemble_analysis_191.mat' ...
       'F:\Imaging\194\NewEnsembles\Ensemble_analysis_194.mat' ...
};
fortitle=[155 157 165 172 174 189 191 194];
summaries=dir('Summary*');
variables=dir('Variables*');
for n=1:length(summaries)
    load (summaries(n).name)
    load (variables(n).name, 'members')
    load (files{n},'EnsRecActIdPlotSt')
    summary{1,7}=  'cells_once';
    summary{1,8}=  length(find(members(:,2)==1))/size(EnsRecActIdPlotSt,1)*100;
    summary{1,9}=  'cells_several';
    summary{1,10}= length(find(members(:,2)>1 & members(:,2)< size(summary,1)))/size(EnsRecActIdPlotSt,1)*100;
    summary{1,11}= 'cells_always';
    summary{1,12}= length(find(members(:,2)==size(summary,1)))/size(EnsRecActIdPlotSt,1)*100;
    save(['Summary_' num2str(fortitle(n)) '.mat'], 'summary');
    clear summary EnsRecActIdPlotSt members
end


    
%%
% % correlation coefficient
% names=dir ('*.mat');
% corrMAT=[];
% corrMAT_sum=[];
% for ii=1:8
% load(names(ii).name);
% NormdSt=x;
% [coeff,score,latent,~,explained]=pca(NormdSt);
% corrMAT(ii,1)=explained(1);
% corrMAT(ii,2)=size(score,2);
% corrMAT_sum(ii,1)=sum(explained(1:3));
% corrMAT_sum(ii,2)=size(score,2);
% clearvars -except  names corrMAT corrMAT_sum
% end
% 
% [R,P] = corrcoef(corrMAT);
% [R_sum,P_sum] = corrcoef(corrMAT_sum);
% figure, for i=1:8, plot(corrMAT(i,1) ,corrMAT(i,2),'b.','MarkerSize',15); hold on; end
% figure, for i=1:8, plot(corrMAT_sum(i,1) ,corrMAT_sum(i,2) ,'r.','MarkerSize',15); hold on; end
% 
% 
% 
% 
% 
% 

% binned values
% [coeff,score,latent,~,explained]=pca(binNormdSt');
% figure,biplot(coeff(:,1:3),'scores',score(:,1:3))
% 
% 
% % 111111
% rawNormdSt=NormdSt;
% rawNormdSt(rawNormdSt>0)=1;
% [coeff,score,latent,~,explained]=pca(rawNormdSt);
% figure,biplot(coeff(:,1:2),'scores',score(:,1:2))
% 
% 
% 







%% the matrix is NOT circulating around itself
% rawNormdSt=NormdSt;
% %rawNormdSt(rawNormdSt>0)=1;
% bins=[1 2 4 8 12 16 20 40];
% timediff=[5 10 30 60 90 120];
% fr=4;
% preSImatrix=nan(1200,length(timediff),length(bins));
% 
% for j=1:length(bins) % for each bin width
%    
%     % create the binned matrix from the original matrix
%     if bins(j)>1                   
%         binNormdSt=zeros(floor(size(rawNormdSt,1)/bins(j)),size(rawNormdSt,2));
%         for k=1:size(binNormdSt,1)
%             binNormdSt(k:k,1:size(rawNormdSt,2))=sum(rawNormdSt( bins(j)*(k-1)+1 : bins(j)*(k)  ,:),1);   %%(bins(j)*k)-(bins(j)-1):bins(j)*k
%         end
%     else
%         binNormdSt=rawNormdSt;
%     end
%    
%     tempbinNormdSt=nan(size(binNormdSt,1),size(binNormdSt,2)) ;
%    
%     for i=1:length(timediff) % for each time difference range
%         for kk=1:size(binNormdSt,1) % for every bin of data
%                
%             ind = floor(1+timediff(i)*(fr/bins(j)));
%             
%             if kk  + ind <=    size(binNormdSt,1)
% 
%                     tempbinNormdSt(1:size(binNormdSt(kk:end,:),1),:)              =      binNormdSt(kk:end,:);     
%                     if kk>1 
%                        tempbinNormdSt(size(binNormdSt(kk:end,:),1)+1 : end,   :)  =      binNormdSt(1:kk-1,:)  ;
%                     end
%                     
%                     
%                           %ind = floor(1+timediff(i)*(fr/bins(j)));
%                     if    ind > size (tempbinNormdSt,1)
%           
%                           preSImatrix(kk,i,j) = (dot(tempbinNormdSt(1,:),tempbinNormdSt(end,:)))/...                     
%                                                ((dot(tempbinNormdSt(1,:),tempbinNormdSt(1,:))  +  (dot(tempbinNormdSt(end,:),tempbinNormdSt(end,:))))    /2); 
%                     else
%                           preSImatrix(kk,i,j) = (dot(tempbinNormdSt(1,:),tempbinNormdSt(ind,:)))/...                     
%                                                ((dot(tempbinNormdSt(1,:),tempbinNormdSt(1,:))  +  (dot(tempbinNormdSt(ind,:),tempbinNormdSt(ind,:))))    /2); 
%           
%                     end
%                     
%             else
%                 continue
%             end
%           
%         end
%         %clear tempbinNormdSt
%     end
% end
% 
% 
% final=[];
% for ii= 1:length(bins)
%     %for jj=1:length(timediff)
%         AVG_final(ii,:)=mean(preSImatrix(:,:,ii),1,'omitnan');
%         STD_final(ii,:)= std(preSImatrix(:,:,ii),1,'omitnan');
%         SEM_final(ii,:)= STD_final(ii,:)/sqrt(size(preSImatrix,1)) ;
%     %end
% end
% % 
% 
%