function REDs=repetitiveIED_v02(MTABS,tmax,Tmin)
% Detector of repetitive epiletiform discharges (REDs)
%
% ---INPUTS---
% MTABS... matrix of absolute time idex correspond to matrix of detection
%          example: MTABS=datenum(0,0,0,0,0,discharges.MP)+tabs(1), where tabs(1) is start of
%          signal in absalute datenum format and discharges.MP is output of
%          spike_detector_hilbert_XX
%          
% tmax ... maximal time betwwen two IEDs (DEFAULT 1 SECOND)
% Tmin ... minimal duration of IED sequence (DEFAULT 10 SECONDS)
%
% ---OUTPUT---
% REDS... structure of detections
%       REDs.N number of detections (qEEG)
%       REDs.meanDUR mean duration of REDs sequence
%       REDs.meanRATE mean IED rate per second during REDs epoch
%       REDs.pos_dur{ch,1} CELL contain absolute time marker of detections,
%                          REDs.pos_dur{ch,1}(:,[start stop])
%   


if exist('tmax','var')==0
    tmax=1; % 1 seconds
end

if exist('Tmin','var')==0
    Tmin=10; % 10 seconds
end


Nmin=round(Tmin/tmax);
if Nmin<2; erroe('Tmin<2*tmin'); end

[~,idx]=sort(min(MTABS,[],2));
MTABS=MTABS(idx,:);


% MW=MW(idx,:);
% MP=MP(idx,:);
% MT=MT(idx,:);
% Q=sum(MW,1)/(total_time/60);

for ch=1:size(MTABS,2)
    
    ied_idx=find(~isnan(MTABS(:,ch)));
    t_abs=MTABS(ied_idx,ch);
    
    DT=diff(t_abs)/datenum(0,0,0,0,0,1); % èas mezi výboji v sekuindách

    bDT=DT<tmax;
    bDT=erodeN(bDT,Nmin-1);
    bDT=dilateN(bDT,Nmin-1);
    
    sidx=[];
    sidx=[find(diff([0; bDT(:)])>0) find(diff([bDT(:);0])<0)];
    
    slength=zeros(size(sidx,1),1);
    if ~isempty(sidx)
        for i=1:size(sidx,1)
%             slength(i,1)=diff(t_abs(sidx(i,:)))/datenum(0,0,0,0,0,1);
            slength(i,1)=sum(DT(sidx(i,1):sidx(i,2)));
        end
        
        idx=find(slength<Tmin);
        slength(idx)=[];
        sidx(idx,:)=[];
    end

    if ~isempty(sidx)
        REDs.N(ch,1)=length(slength);
        REDs.meanDUR(ch,1)=mean(slength);
        REDs.meanRATE(ch,1)=sum(diff(sidx,[],2))/sum(slength);
        REDs.pos_dur{ch,1}=[t_abs(sidx(:,1)) t_abs(sidx(:,2))];
    else
        REDs.N(ch,1)=0;
        REDs.meanDUR(ch,1)=0;
        REDs.meanRATE(ch,1)=0;
        REDs.pos_dur{ch,1}=[];
    end
end

