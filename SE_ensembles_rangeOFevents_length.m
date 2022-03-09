% this is a function to determine the range of the lengths of ensemble activity events

[ind_test1]=find(EnsActSt);
[ind_test]=find(EnsActStAll);
x=zeros(length(ind_test1),1);

for i=2:length(ind_test1)
    x(i-1)=  ind_test(find(ind_test == ind_test1(i))-1)  -    ind_test1(i-1) +1 ;
end

shortest=min(x);
longest=max(x);