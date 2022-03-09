function DI=dunns(clusters_number,distM,ind)   
%%%Dunn's index for clustering compactness and separation measurement
% dunns(clusters_number,distM,ind)
% clusters_number = Number of clusters 
% distM = Dissimilarity matrix
% ind   = Indexes for each data point aka cluster to which each data point
% belongs
i=clusters_number;              % how many clusters
denominator=[];      
for i2=1:i                      % for every cluster
    indi=find(ind==i2);         % find all variables belonging to the current cluster  
    indj=find(ind~=i2);         % find the variables not belonging to it
    x=indi;                     % all those who belong
    y=indj;                     % all those who donot belong
    temp=distM(x,y);            % Eucledian distance between the variables in cluster and the rest (?)
    denominator=[denominator;temp(:)];
end
num=min(min(denominator));      % the lowest distance (numerator) 
neg_obs=zeros(size(distM,1),size(distM,2));
for ix=1:i                      % for each cluster
    indxs=find(ind==ix);        % find all variables belonging to the current cluster  
    neg_obs(indxs,indxs)=1;     % diagonal matrix with diagonal = not related variables to the current cluster
end
dem=neg_obs.*distM;             % 
dem=max(max(dem));              % SHEHAB 
DI=num/dem;
end