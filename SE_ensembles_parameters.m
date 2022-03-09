%percentage of cells recruited in ensembles formation
function [percCellsRecru, RecruCellsID, percEnsDur] = SE_ensembles_parameters(EnsActStAll,EnsRecActStFrames,EnsRecActIdSt,dfoverf0St)

for i=1:size(EnsRecActIdSt,1)                                                               %%% look for cells that were active in any ensemble
    if   ~isempty(find(EnsRecActIdSt(i,:)))
         ind(i)=1;
    else ind(i)=0;
    end
end

percCellsRecru=((sum(ind)) / (length(ind))) *100;                                           %%% percentage of cells recruited in ensembles
RecruCellsID = find (ind);                                                                  %%% IDs of cells recruited in ensembles


for i=1:length(EnsRecActStFrames)                                                           %%% look for the duration of each ensemble event
    EnsActStAll(EnsRecActStFrames(i));
    endofevents = find (EnsActStAll(EnsRecActStFrames(i): length(EnsActStAll))==0);         %%% find when does the current event ends
    if isempty(endofevents)                                                                 %%% in case the current event ends with the end of the video
    DurOfeveEns(i)= length(EnsActStAll(EnsRecActStFrames(i): length(EnsActStAll)));         %%% the duration is till the end of the video 
    else DurOfeveEns(i)= endofevents(1)-1;                                                  %%% otherwise the duration is only until the current event ends
    end
    clear endofevents 
end
percEnsDur = (sum(DurOfeveEns)/size(dfoverf0St,1))*100;                                     %%% percentage of duration of ensembles to duration of video


