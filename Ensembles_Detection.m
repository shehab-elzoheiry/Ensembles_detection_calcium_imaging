%% This script is based on bootstrapping apporach to detect above chance level coincident Ca++ imaging signals based on Hamm et al., Neuron 2017.
% This is combined effort from Shehabeldin Elzoheiry: shehab.elzohairy@gmail.com & Juan C. Boffi: boffi@ana.uni-heidelberg.de

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VARIABLES INDEX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dfoverf0St             = delta fluorescne over basal fluorescence. (n.m) n = length of recording & m = number of cells
% dt                     = sampling rate in miliseconds
% DownCutoffSt	          = lower cutoff, based on the 1st percentile of the randomised (bootstraped) data
% f0St 	                = value of basal fluorescence level of each cell. (1.m) m = number of cells
% FrameAvgSt             = average normalised activity of all cells at each time point. (n.1) n = length of recording
% EnsActIdSt 	          = logical, true at time points AND cells with activity above cutoff. (n.m) n = number of cells & m = number of frames above cutoff 
% EnsActIdShift          = randomly shifted version of EnsActIdSt. (n.m) n = number of cells & m = number of above-cutoff frames. 
% EnsActIdPlotSt         = mark of cells with above cutoff-activity at above-cutoff frames, for plotting red marks in figures. (n.m) n = number of cells & m = number of frames 
% EnsActIdStNonZero      = indices of frames with above-cutoff activity. (1.n) number of frames with above-cutoff activity
% EnsActSt	             = logical, true for first frames of each event contributing to above-cutoff frames. (n.1) n= number of frames
% EnsActStAll	          = logical, true for frames showing average (normalised to itself) cells activity above cutoff. (n.1) n = number of frames
% EnsCaActSt	          = traces of all cells at frames belonging to above-cutoff frames. (n.m) n = number of cells & m = number of frames
% EnsRecActIdSt	       = logical, true for frames that have SIs above upper cutoff of randomised SIs. (n.m) n = number of cells & m = number of frames showing ensembles activity
% EnsRecActStFrames      = index of frames that have SIs above upper cutoff of randomised SIs. (1.n) n = number of frames showing ensembles activity
% EnsRecActIdPlotSt      = for plotting, marks cells recruited in ensembles at times when ensembles are active. (n.m) n = number of cells & m = length of recording
% nCells 		          = number of cells analysed
% nFrames 	             = length of recording (number of frames) (sometimes this variable changes in different sections)
% NormdSt		          = fluorescence of each trace normalised to the maximum of itself. (n.m) n = length of recording & m = number of cells
% NormdFShift	          = RANDOMISED (bootstraped) fluorescence of each trace normalised to the maximum of itself. (n.m) n = length of recording & m = number of cells
% offsetNormdSt          = for plotting, all normalised traces are off-set according to their order. (n.m) n = length of recording & m = number of cells
% offsetEnsActSt         = for plotting, marking the beginings of above-cutoff events. (n.1) n = length of recording.
% PercFramesEnsActSt     = percentage of frames contributing to above-cutoff events out of all frames (whole recording)
% PercFramesEnsRecActSt  = percentage of frames contributing to ensemble's activity out of all frames (whole recording)
% rawSt   	             = raw calcium traces matrix. (n.m) n = length of recording & m = number of cells
% RandFrameAvg	          = average of the 10000 shifts (randomised) of all cells at each frame. (n.10000) n = length of recording (frames)
% RandFrameAvgSt         = concatenated randomised activity values of each cell, to calculate the upper and lower cutoffs. (n.1) n = sum of lengths of concatenated shifts (cells *frames)
% SIactIdSt   	          = similarity index between all above-cutoff frames. (n.n) n = number of above-cutoff frames.
% SIactIdStPlot          = ? vertically concatenated similarity index of all frames. (n.1) n = length of concatenated SIs   
% SIactIdShift	          = similarity index between all (randomly shifted) above-cutoff frames. (n.n) n = number of above-cutoff frames.
% SIactIdShiftPl	       = verticaly concatenated SIs of (randomly shifted) above-cutoff frames. (n.1) n = length of concatenated SIs   
% SIactIdShiftPlo	       = all (the 1000 or 10000) verticaly concatenated SIs of (randomly shifted) above-cutoff frames. (n.m) n = length of concatenated SIs & m = iterations for boostrap
% SIactIdShiftPlot       = vertical concatenation of the 1000 or 10000 vertical concatenations of the SIs of (randomly shifted) above-cutoff frames. (n.1) n = length of all bootstrapped SIs
% sSt 		             = time axis. (1.m) m = time point for each frame
% UpCutoffSt 	          = upper cutoff, based on the 99th percentile of the randomised (bootstraped) data
% UpCutoffSIactSt        = upper 99 percentile of the bootstrapped SIs

%% Automation (Batch Analysis)
   openfile={ '/media/Cleaned_traces_1', ...
              '/media/Cleaned_traces_2', ...
             };
   
   savedir={  '/media/NewEnsembles_1/', ...
              '/media/NewEnsembles_2/', ...
             };
             
   savefile={ '/media/Ensemble_analysis_1', ...
              '/media/Ensemble_analysis_2', ...
             };
%%             
   tic         
for files=1:length(openfile)
    load(openfile{files}, 'clean_traces'); %clearvars -except clean_traces openfile savedir savefile files
    display(files);
%% Import Raw F traces as numeric matrix. Columns = cells, Rows = frames. Name array 'rawSt'.

% Set your adquisition rate
dt = .250;                                                        %%% in milliseconds
if mod(files,2)==0                                                %%% if you are analyzing an even number file (rotenone) 
load(savefile{files-1}, 'rawSt');                                 %%% import the rawSt from the previous (baseline)
rawSt_former=rawSt; clear rawSt                                   %%% now clear it to avoid confusing with the rawSt that is belonging to the file being analysed in this loop
end
clean_traces(end,:)=[];                                           %%%%%%%%% In case I am using original traces from the begining
smooth_traces = movmean(clean_traces,(1/dt),2);                                              %%% smoothing the traces using sliding average
%rawSt=detrend(smooth_traces');                                                              %%% taking offsets away 
rawSt=smooth_traces';

%calculate F0 as 40% quantile
%%%f0St = ((median(rawSt)-min(rawSt))./1.25)+min(rawSt);          %%% choose how do you want to calculate F0 (median -min is 50 then divide by factor to get final threshold)
%f0St = 0.2 .* max(rawSt);                                        %%%In case I am using denoised traces from the begining 

if mod(files,2)==0                                                      %%% if you are analyzing an even number file (rotenone) 
    f0St = 0.2 .*(max(rawSt_former) - min(rawSt_former)) + min(rawSt_former);  %%% adjust it based on the minimum of the previously analyzed (baseline) traces
else
    f0St = 0.2 .*(max(rawSt) - min(rawSt)) + min(rawSt);              %%% In case I am using original traces from the begining
end

dfoverf0St = bsxfun(@minus, rawSt, f0St);                        
for i=1:length(f0St), dfoverf0St(dfoverf0St(:,i)<0,i)=0; end      
    
% Time axis in seconds
sSt = bsxfun(@times, dt, [1:size(dfoverf0St,1)]);

% Normalize dF/F0 to max response of each cell.
NormdSt = bsxfun(@rdivide,dfoverf0St,max(dfoverf0St));                                       %%% Normalise dF/F0 of each cell to the maximum F0 of that cell
NormdSt(isnan(NormdSt))=0;                                                                   %%% cells that became silent in the second phase of analysis (rotenone) will have negative values after using the 
FrameAvgSt = mean(NormdSt,2); %Network activity magnitude = to % of activity ceiling         %%% average activity in each frame to show frames with highest population activity

% Generate random distribution of chance level coactive events
NormdF = NormdSt;
nCells = size(dfoverf0St,2);                                                                 %%% number of cells in your stimulation
nFrames = size(dfoverf0St,1);                                                                %%% number of frames in your video

RandFrameAvg = zeros(nFrames,10000);          
NormdFShift = NormdF;

for ii = 1:10000                                                                             %%% number of iterations done for bootstrapping (to get random distribution)
     for i = 1:nCells                                                                        %%% the bootstraping is done for each cell
         Shift = round((nFrames)*(rand));                                                    %%% rand generates random number (between 0 & 1) every time it runs
         NormdFShift(:,i) = circshift(NormdF(:,i),Shift,1);                                  %%% random circular frame shift for each cell, each displacement is = to the scalar of 'shift' 
     end
     RandFrameAvg(:,ii) = mean(NormdFShift,2);                                               %%% average activity of each individaual frame, to show frames with highest population activity (for randomised population)
end

%Concatenate the 10000 calculations in one vector                                            %%% to enable using the quantile function below
RandFrameAvgSt = [];
for i = 1:nCells
    RandFrameAvgSt = vertcat(RandFrameAvgSt,RandFrameAvg(:,i));                              
end

UpCutoffSt = quantile(RandFrameAvgSt,0.99);                                                  % Determine 99% quantile as upper cutoff chance level coactivity
DownCutoffSt = quantile(RandFrameAvgSt,0.01);                                                % Determine 1% quantile as lower cutoff chance level coactivity


%% Plot chance level coactivity histogram vs data
figure(1)
histogram(RandFrameAvgSt,'NumBins', 50, 'BinWidth', 0.002,'Normalization', 'probability')    %%% blue histogram (randomised data from bootstrapping)
hold on
histogram(FrameAvgSt,'NumBins', 50, 'BinWidth', 0.002,'Normalization', 'probability')        %%% brown histogram (original data)
% xlim([0 0.5])
%plot cutoff lines
line([UpCutoffSt UpCutoffSt], ylim,'color', [0 0 0])                                         %%% vertical line representing the whatever quantile (determined in previous section) cutoff
% line([DownCutoffSt DownCutoffSt], ylim,'color', [0 0 0])
title(['\fontsize{8} Stimulated vs. resampled'])                                             %%% 'Stimulated' is the data under examination
set(gca,'FontSize',8);
xlabel('Frame Average Signal')                                                               
ylabel('Proportion of observed')
legend('show')
legend('Stim. Resampled', 'Stim.', '99% cutoff')
savefig([savedir{files} 'Fig1 PopActVsResampledSt.fig'])                                                          
hold off
close
%% Get frames above upper cutoff
EnsActStAll = FrameAvgSt > UpCutoffSt;                                                      %%% take the frames with average activity above the cutoff for further analysis
EnsActSt = [0; EnsActStAll];                                                                %%% this '0' is added at the beginning, to detect if the first frame active by using 'diff' function below
EnsActSt = diff(EnsActSt)==1;                                                               %%% determine the first frame of each event (composed of several frames), as the time point when the cell was active
PercFramesEnsActSt = ((sum(EnsActStAll))/nFrames)*100;                                      %%% percentage of frames that are involved in ensemble activity
EnsCaActSt = bsxfun(@times,NormdSt',EnsActSt');                                             %%% takes dF/F0 of ALL cells at the begining of each events during ensemble activity (frames above cutoff)

%Get ROI number of cells in above cutoff frames.
EnsActIdSt = zeros(nCells,nFrames);                                                         %%% identifying the cells recruited in ensembles and the time points when they were active
for i = 1:nCells                                                                            %%% for each cell
    for ii = 1:nFrames                                                                      %%% for every frame
        if  EnsCaActSt(i,ii) > UpCutoffSt                                                   %%% search for the cells participating in frames with activity above the cutoff
            EnsActIdSt(i,ii) = 1;
        else
            EnsActIdSt(i,ii) = 0;                                                           %%% excludes cells with activity below cutoff, so that they don't count as members of the ensemble
        end
    end
end

% Get frame n of non zero frames                                                            
EnsActIdStNonZero = (sum(EnsActIdSt,1)) > 1;                
EnsActIdStNonZero = find(EnsActIdStNonZero);                %% indices of frames involved in above-cutoff activity                 

% Remove columns (frames) with all zero elements
EnsActIdSt(:, ~any(EnsActIdSt,1)) = [];                                                     %%% ~isany finds zeros

EnsActIdPlotSt = zeros(nCells,nFrames);
for i = 1:nCells
    for ii = 1:nFrames
        if  EnsCaActSt(i,ii) > UpCutoffSt                                                   %%% search for the cells participating in frames with activity above the cutoff
            EnsActIdPlotSt(i,ii) = i;                                                       %%% mark cells active at frames above the cutoff (by their indices) in a matrix, to plot these points later (as red dots)
        else
            EnsActIdPlotSt(i,ii) = NaN;                                                     %%% empty, so as not to plot anything later at these timepoints for these cells
        end
    end
end

%% Get frames below lower cutoff
% EnsDeactSt = FrameAvgSt > 0; % >0 to ignore stim artifacts with closed shutter.
% EnsDeactSt = bsxfun(@times,FrameAvgSt,EnsDeactSt);
% for i = 1:nFrames
%         if EnsDeactSt(i,1) > 0;
%             EnsDeactSt(i,1) = 1;
%         else
%             EnsDeactSt(i,1) = NaN;
%         end
% end
% EnsDeactSt = bsxfun(@times,FrameAvgSt,EnsDeactSt);
% EnsDeactStAll = EnsDeactSt < DownCutoffSt;
% EnsDeactSt = [0; EnsDeactStAll];
% EnsDeactSt = diff(EnsDeactSt)==1; % Get first frame of events that span consecutive frames
% PercFramesEnsDeactSt = ((sum(EnsDeactStAll))/nFrames)*100;
% 
% % Get dFoverF0 values during EnsDeact
% EnsCaDeactSt = bsxfun(@times,NormdSt',EnsDeactSt');
% 
% %Get ROI number of cells in below cutoff frames.
% EnsDeactIdSt = zeros(nCells,nFrames);
% for i = 1:nCells
%     for ii = 1:nFrames
%         if EnsCaDeactSt(i,ii) < DownCutoffSt && EnsCaDeactSt(i,ii) > 0; %below threshold
%             EnsDeactIdSt(i,ii) = 1;
%         else
%             EnsDeactIdSt(i,ii) = 0;
%         end
%     end
% end
% 
% % Get frame n of non zero frames
% EnsDeactIdStNonZero = (sum(EnsDeactIdSt,1)) > 1;
% EnsDeactIdStNonZero = find(EnsDeactIdStNonZero);
% 
% % Remove columns (frames) with all zero elements
% EnsDeactIdSt(:, ~any(EnsDeactIdSt,1)) = [];
% 
% EnsDeactIdPlotSt = zeros(nCells,nFrames);
% for i = 1:nCells
%     for ii = 1:nFrames
%         if EnsCaDeactSt(i,ii) < DownCutoffSt && EnsCaDeactSt(i,ii) > 0; %below threshold
%             EnsDeactIdPlotSt(i,ii) = i;
%         else
%             EnsDeactIdPlotSt(i,ii) = NaN;
%         end
%     end
% end

%% Plot Ens act 

%Offset normalized Ca traces
offset = 1:nCells;
offsetNormdSt = bsxfun(@plus, NormdSt, offset);                                            %%% for the sake of plotting, each normalised trace is offset by its order to present each cell separately

offset = nCells+2;
offsetEnsActSt = bsxfun(@plus, EnsActSt, offset);                                          %%% generating an extra offset above the last trace, to mark the frame (above the cutoff) at which an event starts

for i = 1:nFrames                                                                          %%% this for-loop is a trick to help marking the frames above cutoff just after the trace of last cell
    if offsetEnsActSt(i) == offset
       offsetEnsActSt(i) = NaN;
    end
end
% offsetEnsDeactSt = bsxfun(@plus, EnsDeactSt, offset);
% for i = 1:nFrames
%     if offsetEnsDeactSt(i) == offset
%        offsetEnsDeactSt(i) = NaN;
%     end
% end

figure(2)
plot(sSt,offsetNormdSt,'color',[0.4 0.7 1])                                                 %%% plots normalised traces of all cells
hold on
plot(sSt,(EnsActIdPlotSt),'r.','MarkerSize',7)                                              %%% marks (red dots) cells contributing to an above-cutoff frame
% plot(sSt,(EnsDeactIdPlotSt),'b.','MarkerSize',4)
plot(sSt,(offsetEnsActSt),'r','Marker', '+')                                                %%% marks (red plus signs) frames that are above cutoff
% plot(sSt,(offsetEnsDeactSt),'b','Marker', '+')
xlabel('time (s)')
ylabel('cell #')
savefig([savedir{files} 'Fig2 EnsCaSt.fig'])
hold off
close
%% Calculate Similarity index (SI = cosine of the angle between the vectors) between each Id vector pair)    
                                                                                     %%%% is the heading correct about using the COSINE ??
%Act

nFrames = size(EnsActIdSt,2);                                                               %%% uses the frames that are above the cutoff
SIactIdSt = zeros(nFrames,nFrames);                                                         %%% creates a sqaure shape matrix for ploting SI between all pairs of frames (that are above the cutoff) 
parfor i = 1:nFrames
    for ii = 1:nFrames
        SIactIdSt(i,ii) = (dot(EnsActIdSt(:,i),EnsActIdSt(:,ii)))/...                       %%% this is the equation for calculating similarity index -->  SI= Ca.Cb / ((|Ca|²+|Cb|²) / 2)
            ( (dot(EnsActIdSt(:,i),EnsActIdSt(:,i))  +  (dot(EnsActIdSt(:,ii),EnsActIdSt(:,ii))))    /2);  
    end
end

%Concatenate the SI calculations in one vector                                        
SIactIdStPlot = [];                                                                 
for i = 1:nFrames                                                                    
    SIactIdStPlot = vertcat(SIactIdStPlot,SIactIdSt(i+1:nFrames,i));                  
end

% %Deact
% nFrames = size(EnsDeactIdSt,2);
% SIdeactIdSt = zeros(nFrames,nFrames);
% parfor i = 1:nFrames
%     for ii = 1:nFrames
%         SIdeactIdSt(i,ii) = (dot(EnsDeactIdSt(:,i),EnsDeactIdSt(:,ii)))/...
%             ((dot(EnsDeactIdSt(:,i),EnsDeactIdSt(:,i))+(dot(EnsDeactIdSt(:,ii),EnsDeactIdSt(:,ii))))/2);
%     end
% end
% 
% %Concatenate the SI calculations in one vector
% 
% SIdeactIdStPlot = [];
% for i = 1:nFrames
%     SIdeactIdStPlot = vertcat(SIdeactIdStPlot,SIdeactIdSt(i+1:nFrames,i));
% end


%% Determine SI of resampled random chance level ensembles                                       %%% SIactIdShiftPl --> SIactIdShiftPlo --> SIactIdShiftPlot

%Generate random shift along frames and evaluate SI of rearranged random ensembles              
%Act

EnsActIdShift = EnsActIdSt;
nFrames = size(EnsActIdSt,2);
SIactIdShiftPlo = zeros((((nFrames^2)-nFrames)/2),10000);                                         %%% can reduce to 1000. ((n^2)-n)/2 is an equation to get the right length of the 1st dim for concatenation below

for iii = 1:10000 % can be reduced to 1000, times works OK.
    for i = 1:nCells
        Shift = round((nFrames)*(rand));                                                         %%% creating random magnitudes of shifts
        EnsActIdShift(i,:) = circshift(EnsActIdShift(i,:),Shift,2);                                % for each cell, random circular frame shift by the magnitude of 'shift'  
    end
    
    % Calculate SI
    SIactIdShift = NaN(nFrames,nFrames);                                                         %%% creates a sqaure shape matrix for ploting SI between all pairs of frames (randomized above-cutoff frames) 
    for i = 1:nFrames
        for ii = i+1:nFrames                                                                     %%% this is a trick to generate half the SI plot w/o its mirror image, to reduce computational time
        SIactIdShift(i,ii) = (dot(EnsActIdShift(:,i),EnsActIdShift(:,ii)))/...                   %%% this is the equation for calculating similarity index -->  SI = Ca.Cb / ((|Ca|²+|Cb|²) / 2)
            (  (dot(EnsActIdShift(:,i),EnsActIdShift(:,i))   +   (dot(EnsActIdShift(:,ii),EnsActIdShift(:,ii))))     /2);
        end
    end
    
    SIactIdShiftPl = [];
    for i = 1:nFrames
        SIactIdShiftPl = vertcat(SIactIdShiftPl,SIactIdShift(i,((i+1):nFrames))');               %%% vertical concatenation of the outcome from one iteration
    end
    SIactIdShiftPlo(:,iii) = SIactIdShiftPl;                                                     %%% outcome of each bootstrap iteration will be stored through the 2nd dim
end

%Concatenate the SI calculations in one vector
SIactIdShiftPlot = [];
for i = 1:10000 % can be reduced to 1000 times, works OK.
    SIactIdShiftPlot = vertcat(SIactIdShiftPlot,SIactIdShiftPlo(:,i));                           %%% vertical concatenation of the 10000 iterations (bootstraping) to calculate the quantile below
end

% Determine 99% quantile as upper cutoff chance level recurrent coactivity
UpCutoffSIactSt = quantile(SIactIdShiftPlot,0.99);                                               %%% calculating the whatever quantile of the randomly distributed SI 

% %Deact
% EnsDeactIdShift = EnsDeactIdSt;
% nFrames = size(EnsDeactIdSt,2);
% SIdeactIdShiftPlo = zeros((((nFrames^2)-nFrames)/2),10000);
% 
% for iii = 1:10000
%     for i = 1:nCells
%         Shift = round((nFrames)*(rand));
%         EnsDeactIdShift(i,:) = circshift(EnsDeactIdShift(i,:),Shift,2); %random circular frame shift for each cell
%     end
%     % Calculate SI
%     SIdeactIdShift = NaN(nFrames,nFrames);
%     for i = 1:nFrames
%         for ii = i+1:nFrames
%         SIdeactIdShift(i,ii) = (dot(EnsDeactIdShift(:,i),EnsDeactIdShift(:,ii)))/...
%             ((dot(EnsDeactIdShift(:,i),EnsDeactIdShift(:,i))+(dot(EnsDeactIdShift(:,ii),EnsDeactIdShift(:,ii))))/2);
%         end
%     end
%     SIdeactIdShiftPl = [];
%     for i = 1:nFrames
%         SIdeactIdShiftPl = vertcat(SIdeactIdShiftPl,SIdeactIdShift(i,((i+1):nFrames))');
%     end
%     SIdeactIdShiftPlo(:,iii) = SIdeactIdShiftPl;
% end
% 
% %Concatenate the SI calculations in one vector
% SIdeactIdShiftPlot = [];
% for i = 1:10000
%     SIdeactIdShiftPlot = vertcat(SIdeactIdShiftPlot,SIdeactIdShiftPlo(:,i));
% end
% 
% % Determine 99% quantile as upper cutoff chance level recurrent coactivity
% UpCutoffSIdeactSt = quantile(SIdeactIdShiftPlot,0.99); 

%% Plot chance level recurrent coactivity histogram vs data
%Act
figure(3)
histogram(SIactIdShift,'NumBins', 50, 'BinWidth', 0.02,'Normalization', 'probability')              %%% histogram of randomised SIs (in blue)
hold on
histogram(SIactIdSt,'NumBins', 50, 'BinWidth', 0.02,'Normalization', 'probability')                 %%% histogram of non-randomised SIs (in brown)
xlim([0 1])
% ylim([0 0.05])
%plot cutoff lines
line([UpCutoffSIactSt UpCutoffSIactSt], ylim,'color', [0 0 0])                                      %%% vertical line representing the upper cutoff (based on randomised SIs)
title(['\fontsize{8} Stimulated vs. resampled'])
set(gca,'FontSize',8);
xlabel('Frame-Frame similarity')
ylabel('Proportion of observed')
legend('show')
legend('Stim. Resampled', 'Stim.', '99% cutoff')
savefig([savedir{files} 'Fig3 PopRecurActVsResampledSt.fig'])
hold off
close
% figure(4)
% histogram(SIdeactIdShift,'NumBins', 50, 'BinWidth', 0.02,'Normalization', 'probability')
% hold on
% histogram(SIdeactIdSt,'NumBins', 50, 'BinWidth', 0.02,'Normalization', 'probability')
% xlim([0 1])
% % ylim([0 0.07])
% %plot cutoff lines
% line([UpCutoffSIdeactSt UpCutoffSIdeactSt], ylim,'color', [0 0 0])
% title(['\fontsize{8} Stimulated vs. resampled'])
% set(gca,'FontSize',8);
% xlabel('Frame-Frame similarity')
% ylabel('Proportion of observed')
% legend('show')
% legend('Stim. Resampled', 'Stim.', '99% cutoff')
% savefig('PopRecurDeactVsResampledSt.fig')
% hold off

%% Get frames above upper SI cutoffs                            
%%% 'EnsRecActIdPlotSt' is a full-video-sized matrix with frames marked when ensembles are active // 'EnsRecActStFrames' is list of indices of frames when ensmbles are active  

% Percentage of frames above cutoff ACT

% Get Ens act above cutoff for plotting.
EnsRecActStFrames = SIactIdSt > UpCutoffSIactSt;                                      %%% getting co-active frames above the upper cutoff (based on randomised SIs)
EnsRecActStFrames = (sum(EnsRecActStFrames,1)) > 1;                                   %%% produce logical to identify frames having above-cutoff SIs
EnsRecActIdSt = bsxfun(@times,EnsActIdSt,EnsRecActStFrames);                          %%% getting the cells and time points related to co-active frames above the upper cutoff (based on randomised SIs)
EnsRecActIdSt(:,~any(EnsRecActIdSt,1)) = [];                                            % Remove non recurrent columns
EnsRecActStFrames = bsxfun(@times,EnsActIdStNonZero,EnsRecActStFrames);                 % switch from logical to frames indices (having above-cutoff SIs)
EnsRecActStFrames(:,~any(EnsRecActStFrames,1)) = [];                                    % Remove indices for frames not participating in ensemble activity

nFrames = size(dfoverf0St,1);
EnsRecActIdPlotSt = NaN(nCells,nFrames); 
for i = 1:nFrames % 19916
    for ii = 1:size(EnsRecActStFrames,2) %%409
        if i == EnsRecActStFrames(ii)                                                 %%% if the frame (in the full video) belongs to a time point when the ensemble was active, proceed
           EnsRecActIdPlotSt(:,i) = EnsActIdPlotSt(:,EnsRecActStFrames(ii));          %%% mark time points when SIs are > cutoff (based on randomised SIs) from cells at frames > 1st boostrap cutoff ('EnsActIdPlotSt')
        end
    end
end

PercFramesEnsRecActSt = ((sum((sum(EnsRecActIdPlotSt,'omitnan')) > 0))/nFrames)*100;  %%% ????????? why not just the size(EnsRecActStFrames,2)
% 
% % Percentage of frames above cutoff DEACT
% 
% % Get Ens deact above cutoff for plotting.
% EnsRecDeactStFrames = SIdeactIdSt > UpCutoffSIdeactSt;
% EnsRecDeactStFrames = (sum(EnsRecDeactStFrames,1)) > 1;
% EnsRecDeactIdSt = bsxfun(@times,EnsDeactIdSt,EnsRecDeactStFrames);
% EnsRecDeactIdSt(:,~any(EnsRecDeactIdSt,1)) = []; % Remove non recurrent columns
% EnsRecDeactStFrames = bsxfun(@times,EnsDeactIdStNonZero,EnsRecDeactStFrames);
% EnsRecDeactStFrames(:,~any(EnsRecDeactStFrames,1)) = []; % Remove non recurrent columns
% 
% nFrames = size(dfoverf0St,1);
% EnsRecDeactIdPlotSt = NaN(nCells,nFrames);
% for i = 1:nFrames
%     for ii = 1:size(EnsRecDeactStFrames,2)
%         if i == EnsRecDeactStFrames(ii)
%            EnsRecDeactIdPlotSt(:,i) = EnsDeactIdPlotSt(:,EnsRecDeactStFrames(ii));
%         end
%     end
% end
% 
% PercFramesEnsRecDeactSt = ((sum((sum(EnsRecDeactIdPlotSt,'omitnan')) > 0))/nFrames)*100;
% 
%% Plot Recurrent Ens act

%Offset normalized Ca traces
offset = 1:nCells;
offsetNormdSt = bsxfun(@plus, NormdSt, offset);                         %%% for the sake of plotting, each normalised trace is offset by its order to present each cell separately

offset = nCells+2;
offsetEnsActSt = bsxfun(@plus, EnsActSt, offset);                       %%% generating an extra offset above the last trace, to mark the frames (above the cutoff) at which an event starts    
for i = 1:nFrames                                                       %%% this for-loop is a trick to help marking the frames above cutoff just after the trace of last cell
    if offsetEnsActSt(i) == offset
       offsetEnsActSt(i) = NaN;
    end
end
% offsetEnsDeactSt = bsxfun(@plus, EnsDeactSt, offset);
% for i = 1:nFrames
%     if offsetEnsDeactSt(i) == offset
%        offsetEnsDeactSt(i) = NaN;
%     end
% end

figure(5)
plot(sSt,offsetNormdSt,'color',[0.4 0.7 1])                             %%% plots normalised traces of all cells
hold on
plot(sSt,(EnsRecActIdPlotSt),'r.','MarkerSize',7)                       %%% marks (red dots) cells recruited in an ensemble
% plot(sSt,(EnsRecDeactIdPlotSt),'b.','MarkerSize',4)
plot(sSt,(offsetEnsActSt),'r','Marker', '+')                            %%% marks (red plus signs) frames when ensembles are active  ?????? what is the difference then between this and previous figure ?????????????????
% plot(sSt,(offsetEnsDeactSt),'b','Marker', '+')
xlabel('time (s)')
ylabel('cell #')
savefig([savedir{files} 'Fig5 EnsRecCaSt.fig'])
hold off
close
%% Plot SI matrices

%Act
figure(6)
imagesc(SIactIdSt, [0 1])
axis image
axis square
axis vis3d
title(['\fontsize{8} SI of ensemble activations'])
set(gca,'FontSize',8);
xlabel('Co-active event #')
ylabel('Co-active event #')
savefig([savedir{files} 'Fig6 SIactSt.fig'])
close

figure(7)
imagesc(SIactIdSt, [UpCutoffSIactSt 1])
axis image
axis square
axis vis3d
title(['\fontsize{8} SI of recurrent ensemble activations'])
set(gca,'FontSize',8);
xlabel('Co-active event #')
ylabel('Co-active event #')
savefig([savedir{files} 'Fig7 SIrecurrentActSt.fig'])
close
% %Deact
% figure(8)
% imagesc(SIdeactIdSt, [0 1])
% axis image
% axis square
% axis vis3d
% title(['\fontsize{8} SI of ensemble deactivations'])
% set(gca,'FontSize',8);
% xlabel('Co-inactive event #')
% ylabel('Co-inactive event #')
% savefig('SIdeactSt.fig')
% 
% figure(9)
% imagesc(SIdeactIdSt, [UpCutoffSIdeactSt 1])
% axis image
% axis square
% axis vis3d
% title(['\fontsize{8} SI of recurrent ensemble deactivations'])
% set(gca,'FontSize',8);
% xlabel('Co-inactive event #')
% ylabel('Co-inactive event #')
% savefig('SIrecurrentDeactSt.fig')

%%
%% Determine Ensemble diversity


if  isempty(EnsRecActIdSt)
    RelNumClustEnsActSt=0;
else
    %Act
    % Hierarchical tree
    RecActIdTreeSt = linkage(EnsRecActIdSt', 'weighted','cosine');

    %Cluster data and get relative diversity of ensembles
    ClustersActSt = cluster(RecActIdTreeSt,'criterion','distance','cutoff',UpCutoffSIactSt);
    RelNumClustEnsActSt = (max(ClustersActSt))/(size(EnsRecActIdSt,2));

    % Plot dendrogram
    figure(10)
    dendrogram(RecActIdTreeSt,0)
    line(xlim, [UpCutoffSIactSt UpCutoffSIactSt],'color', [1 0 0])
    ylim([0 1])
    set(gca,'xtick',[])
    title(['\fontsize{8} Diversity of ensemble activations'])
    set(gca,'FontSize',8);
    ylabel('SI')
    savefig([savedir{files} 'Fig10 DendrogramEnsActSt.fig'])
    close

    % %Deact
    % % Hierarchical tree
    % RecDeactIdTreeSt = linkage(EnsRecDeactIdSt', 'weighted','cosine');
    % 
    % %Cluster data and get relative diversity of ensembles
    % ClustersDeactSt = cluster(RecDeactIdTreeSt,'criterion','distance','cutoff',UpCutoffSIdeactSt);
    % RelNumClustEnsDeactSt = (max(ClustersDeactSt))/(size(EnsRecDeactIdSt,2));
    % 
    % % Plot dendrogram
    % figure(11)
    % dendrogram(RecDeactIdTreeSt,0)
    % line(xlim, [UpCutoffSIdeactSt UpCutoffSIdeactSt],'color', [1 0 0])
    % set(gca,'xtick',[])
    % title(['\fontsize{8} Diversity of ensemble deactivations'])
    % set(gca,'FontSize',8);
    % ylabel('SI')
    % savefig('DendrogramEnsDeactSt.fig')
end
%% Result vector with:
% % of frames with higher than chance level co-active events
% % of frames with higher than chance level recurrent co-active events (Ensembles)
% Diversity index of ensembles

%%%%%%%%%%%%%%%%%%% 
if  isempty(EnsRecActIdSt)
    percCellsRecru=0; RecruCellsID=[]; percEnsDur=0;general_sim=0;thresh4clust=[];
else
[percCellsRecru, RecruCellsID, percEnsDur] = SE_ensembles_parameters(EnsActStAll,EnsRecActStFrames,EnsRecActIdSt,dfoverf0St);
[thresh4clust,general_sim] =                 SE_ensembles_thresh(EnsRecActIdSt,1000,0.95,RecActIdTreeSt);                           %%% SE_ensembles_thresh(EnsRecActIdSt,iterations,quartile,RecActIdTreeSt)
end
%AUC=sum(NormdSt,1); meanAUCs=mean(sum(NormdSt,1)); medianAUCs=median(sum(NormdSt,1));                                               %%% old way calculating area under the curve of the normalised traces, additionally the mean and median
AUC=trapz(clean_traces,2); meanAUCs=mean(trapz(clean_traces,2)); medianAUCs=median(trapz(clean_traces,2));                                          %%% new way calculating AUC
EnsAUC=sum(NormdSt(:,RecruCellsID),1); meanEnsAUCs=mean(sum(NormdSt(:,RecruCellsID),1)); medianEnsAUCs=median(sum(NormdSt(:,RecruCellsID),1));  %%% same but for cells recruited in ensembles
[percDepolCells , diffMaxMinEvent] = SE_ensembles_depolCells(NormdSt);
%%%%%%%%%%%%%%%%%%%

Result(1,:)  =    {'PercFramesEnsRecActSt', PercFramesEnsRecActSt};
Result(2,:)  =    {'RelNumClustEnsActSt', RelNumClustEnsActSt};
Result(3,:)  =    {'percCellsRecru', percCellsRecru};
Result(4,:)  =    {'RecruCellsID', RecruCellsID};
Result(5,:)  =    {'general_sim', general_sim};
Result(6,:)  =    {'thresh4clust', thresh4clust};
Result(7,:)  =    {'percEnsDur', percEnsDur};
Result(8,:)  =    {'AUC', AUC};
Result(9,:)  =    {'meanAUCs', meanAUCs };
Result(10,:) =    {'medianAUCs', medianAUCs};
Result(11,:) =    {'EnsAUC', EnsAUC };
Result(12,:) =    {'meanEnsAUCs', meanEnsAUCs };
Result(13,:) =    {'medianEnsAUCs', medianEnsAUCs};
Result(14,:) =    {'percDepolCells', percDepolCells};
Result(15,:) =    {'diffMaxMinEvent', diffMaxMinEvent};

%%
save(savefile{files});
display(['Video number ' num2str(files) ' is analysed and ensembles are saved ;)']);
close all
clearvars -except openfile savedir savefile
end
toc
