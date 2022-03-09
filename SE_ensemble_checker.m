

function [output] = SE_ensemble_checker(     q  ,   coactive_frames,   SIactIdSt, UpCutoffSIactSt)
    for i=1:length(q)
        z = coactive_frames(q(i), find(coactive_frames(q(i),:)));               % check how many elements in this member's group
        for ii=1:length(find(coactive_frames(q(i),:)))                          % for the number of elements in his group
            x = SIactIdSt(coactive_frames(q(i),ii),[z]) > UpCutoffSIactSt;      % check each element with the others
            if   all(x)                                                         % if all are fine
                 check(ii)=1;                                                   % then overall for this element is true
            else check(ii)=0;                                                   % otherwise false
            end
        end
        if   all(check)                                                         % if all are true then kick him out
             output(i)=nan; 
        else output(i)=q(i);
        end
    end



