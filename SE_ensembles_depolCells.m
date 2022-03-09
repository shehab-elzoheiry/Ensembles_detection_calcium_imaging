function [percDepolCells , diffMaxMinEvent] = SE_ensembles_depolCells(NormdSt)
guide=nan(1,size(NormdSt,2));
for i=1:size(NormdSt,2)
    xx= find (NormdSt(:,i)==0);         %find traces that had zeros at any time point
    if isempty(xx)                      
        guide(i)=1;
    else guide(i)=0;
    end
    clear xx
end

percDepolCells=100*(sum(guide)/length(guide));

xx= find (guide==1);
for i=1:length(xx)
    if  isempty(xx)
        continue
    end
    maximum(i)=max(NormdSt(:,xx(i)));   %max. value
    minimum(i)=min(NormdSt(:,xx(i)));   %min. value
end
if isempty(xx)
diffMaxMinEvent = ['No cells that are constantly depolarised'];
else
diffMaxMinEvent = maximum-minimum;           %biggest difference
end 
% for i=1:nCells
% line([1 300],[i i],'LineWidth',0.01)
% end