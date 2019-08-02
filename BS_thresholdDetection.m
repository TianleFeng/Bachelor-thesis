dirData = dir('*.nmd'); 
numfiles = numel(dirData);
allresults = strings([0,1]);
for k = 1:numfiles
    matFileName = sprintf('verified_goldstd_%d.nmd', k);
	if exist(matFileName, 'file')
        load(matFileName, '-mat')
	else
		fprintf('File %s does not exist.\n', matFileName);
        continue
    end
    
    combieeg=[];                              %concatenate EEG vectors
    for v=1:length(rcd.eeg)
        if ~isempty(rcd.eeg{v})
            combieeg=[combieeg;rcd.eeg{v}(:,1)];
        end
    end
     
    combigs=[];                               %concatenate gold standard vectors
    for v=1:length(rcd.verified)
        if size(rcd.verified{v}, 1) > 1
            combigs=[combigs;rcd.verified{v}(:,1)];
        end
    end
%%
counter=1;
for a=6                                       %amplitude threshold
    threshold=BS_threshold(combieeg,200,2,a); %BS_threshold(eegdata,sr,locutoff,amp_th)
    amp_th_uV(counter,1)=a;
    
    mcompare=threshold-combigs;               %substract gold standard from detection vector
    FP=find(mcompare==1);                     %false positive
    FN=find(mcompare==-1);                    %false negative
    FPratio(counter,1)=length(FP)/length(combigs);
    FNratio(counter,1)=length(FN)/length(combigs);
    
    pcompare=threshold+combigs;               %add gold standard and detection vector together
    TP=find(pcompare==2);                     %true positive
    TN=find(pcompare==0);                     %true negative
    sensi(counter,1)=length(TP)/(length(TP)+length(FN));
    speci(counter,1)=length(TN)/(length(TN)+length(FP));
    qual=sensi+speci;                         %quality indicator
    
    evatable_th = [amp_th_uV FPratio FNratio sensi speci qual];
    
    counter=counter+1;
end
%%
col=cellstr(["amp_th_uV" "FPratio" "FNratio" "sensi" "speci" "qual"]);
T=array2table(evatable_th, 'VariableNames', col);
[~,maxidx]=max(T.qual);
%create a table that contains the best qual with amplitude threshold, sensi and speci
opt_thresh_metrics=T(maxidx, :)
opt_thresh=opt_thresh_metrics.amp_th_uV;
allresults=[allresults;'Optimal threshold for file ' matFileName ' was set at ' num2str(opt_thresh) ' uV. Accuracy was ' num2str(opt_thresh_metrics.qual*50, 3) '%. Sensitivity ' num2str(opt_thresh_metrics.sensi) ' , specificity ' num2str(opt_thresh_metrics.speci) '. ']
end