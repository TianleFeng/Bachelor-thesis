load('verified_goldstd_9.nmd','-mat')    
    combieeg=[];                                %concatenate EEG vectors
    for v=1:length(rcd.eeg)
        if ~isempty(rcd.eeg{v})
            combieeg=[combieeg;rcd.eeg{v}(:,1)];
        end
    end
    
    combigs=[];                                 %concatenate gold standard vectors
    for v=1:length(rcd.verified)
        if size(rcd.verified{v}, 1) > 1
            combigs=[combigs;rcd.verified{v}(:,1)];
        end
    end
    
    fcombieeg=eegfilt(combieeg',200,8,47);      %eeg filtering

counter=1;
for a=1050                                      %Amplitude threshold
    nleo=NLEO(fcombieeg,100);                   %nleo=NLEO(EEGdata,halfsr)
    nleo(abs(nleo)>a)=nan;                          
    nleo(1)=nan;
    nsupp=find(isnan(nleo));
    diffns=diff(nsupp);
    suppsec=find(diffns>=200);                  %if the length between 2 NaNs >=1 sec will be detected as suppression
    suppvect=zeros(length(fcombieeg),1);
    for i=1:length(suppsec)
        start=nsupp(suppsec(i));
        stop=nsupp(suppsec(i)+1);
        suppvect(start:stop)=1;
    end

    amp_th_uV(counter,1)=a;
    
    mcompare=suppvect-combigs;                  %substract gold standard from detection vector
    FP=find(mcompare==1);                       %false positive
    FN=find(mcompare==-1);                      %false negative
    FPratio(counter,1)=length(FP)/length(combigs);
    FNratio(counter,1)=length(FN)/length(combigs);
    
    pcompare=suppvect+combigs;                  %add gold standard and detection vector together
    TP=find(pcompare==2);                       %true positive
    TN=find(pcompare==0);                       %true negative
    sensi(counter,1)=length(TP)/(length(TP)+length(FN));
    speci(counter,1)=length(TN)/(length(TN)+length(FP));
    qual=sensi+speci;                           %quality indicator
    
    evatable_nleo = [amp_th_uV FPratio FNratio sensi speci qual];
    
    counter=counter+1;
end

col=cellstr(["amp_th_uV" "FPratio" "FNratio" "sensi" "speci" "qual"]);
T=array2table(evatable_nleo, 'VariableNames', col);
[~,maxidx]=max(T.qual);                         
%create a table that contains the best qual with amplitude threshold, sensi and speci
opt_thresh_metrics=T(maxidx, :)
opt_thresh=opt_thresh_metrics.amp_th_uV;
results=['Optimal threshold was set at ' num2str(opt_thresh) ' uV. Accuracy was ' num2str(opt_thresh_metrics.qual*50, 3) '%. Sensitivity ' num2str(opt_thresh_metrics.sensi) ' , specificity ' num2str(opt_thresh_metrics.speci) '. ']
