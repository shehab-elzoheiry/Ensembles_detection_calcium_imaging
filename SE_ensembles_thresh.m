function [thresh4clust,general_sim]=SE_ensembles_thresh(EnsRecActIdSt,iterations,quartile,RecActIdTreeSt)

EnsRecActIdSt_shift=EnsRecActIdSt;
nFrames=size(EnsRecActIdSt,2);                                                                              %%% number of frames that are going to be compared to each other
nCells=size(EnsRecActIdSt,1);                                                                               %%% number of cells 

dendro_shift=zeros(nFrames-1,iterations);                                                                   %%% preallocated matrix, collecting the results of all iterations at the end of for-loop

for iterate=1:iterations                                                                                          
    for i = 1:nFrames
            Shift = round((nCells)*(rand));                                                                 %%% randomly generated shifts                             
            EnsRecActIdSt_shift(:,i) = circshift(EnsRecActIdSt_shift(:,i),Shift,1);                         %%% applying the shifts within each frame (switching the cells that are active at each time point)
    end
    RecActIdTreeSt_shift = linkage(EnsRecActIdSt_shift','weighted','cosine');                               %%% heirarchical clustering 
    dendro_shift(:,iterate)=RecActIdTreeSt_shift(:,3);                                                      %%% gathering values of clustering for each iteration
end

dendro_shift_plot = [];
for ii = 1:iterations
    dendro_shift_plot = vertcat(dendro_shift_plot,dendro_shift(:,ii));                                      %%% concatenating all values vertically to determine the whatever quartile in next step
end

UpCutoff = quantile(dendro_shift_plot,quartile);                                                            %%% calculating the whatever quartile to be a uppercutoff for clustering
figure, histogram(dendro_shift,'NumBins', 50, 'BinWidth', 0.02,'Normalization', 'probability')              %%% plotting the bootstrap results    
        hold on
        histogram(RecActIdTreeSt(:,3),'NumBins', 50, 'BinWidth', 0.02,'Normalization', 'probability')       %%% plotting the original result
        line([UpCutoff UpCutoff], ylim,'color', [0 0 0]);                                                   %%% marking the threshold
        
thresh4clust=RecActIdTreeSt(min(find(RecActIdTreeSt(:,3)>UpCutoff)),3);                                     %%% the threshold is the first value above the uppercutoff

dendrogram(RecActIdTreeSt,0)
if   isempty(thresh4clust)
     display ('No threshold calculated :( ');
else line(xlim, [thresh4clust thresh4clust],'color', [1 0 0])
end
ylim([0 1])
%set(gca,'xtick',[])
title(['\fontsize{8} Diversity of ensemble activations'])
set(gca,'FontSize',8);
ylabel('SI')

general_sim=mean(RecActIdTreeSt(:,3));                                                                      %%% general similarity of the whole dendrogram
