function clustering=clustering_main_fun(EVT,duration)
amp='ma';
pca_centering=0.4;

ch_err=(sum(isnan(EVT.MA),1)==0); % used channels=1
HV=[];
LV=[];
for ch=1:size(EVT.MA,2)
    Q3=quantile(EVT.MA(EVT.MV(:,ch)>0,ch),0.75);
    Q1=quantile(EVT.MA(EVT.MV(:,ch)>0,ch),0.25);
    HV(1,ch)=Q3+3*(Q3-Q1); 
    LV(1,ch)=Q1-3*(Q3-Q1); 
end
HV=repmat(HV,size(EVT.MA,1),1);
LV=repmat(LV,size(EVT.MA,1),1);
EVT.MA(EVT.MA>HV)=HV(EVT.MA>HV);
EVT.MA(EVT.MA<LV)=LV(EVT.MA<LV);

removed_lines=sort(find(sum(EVT.MA(:,ch_err),2)==0));

EVT.MA(removed_lines,:)=[];
EVT.MV(removed_lines,:)=[];
EVT.MW(removed_lines,:)=[];
EVT.MP(removed_lines,:)=[];



% klastrování s rekurzí -------------------------------------------------------------
 

% potential profile choice
switch amp
    case 'ma'
        CL=spike_cluster_v12(EVT.MA(:,ch_err),EVT.MW(:,ch_err).*double(EVT.MV(:,ch_err)>0),pca_centering); % merging 95%
    case 'raw'
        CL=spike_cluster_v12(EVT.MRAW(:,ch_err),EVT.MW(:,ch_err).*double(EVT.MV(:,ch_err)>0),pca_centering); % merging 95%
end



if ~isempty(removed_lines)
    for i=1:length(removed_lines)
        CL.class=[CL.class(1:removed_lines(i)-1); 0; CL.class(1:removed_lines(i))];
        CL.weight=[CL.weight(1:removed_lines(i)-1); 0; CL.weight(1:removed_lines(i))];
        CL.area=[CL.area(1:removed_lines(i)-1); 0; CL.area(1:removed_lines(i))];
    end
end

c_unique=sort(unique(CL.class(:))); c_unique(c_unique==0)=[]; c_unique=[c_unique; 0];

reds_tmin=1; reds_Tmin=10;

for c_n=1:length(c_unique);
    if sum(CL.class==c_unique(c_n))>0
        MQ=EVT.MW(CL.class==c_unique(c_n),:).*double((EVT.MV(CL.class==c_unique(c_n),:))>0);
        clustering.qIED(c_n,:)=sum(MQ,1)/(duration/60); % 1/min
        clustering.qSINGLE(c_n,:)=sum(MQ((sum(MQ>0,2)==1),:),1)/(duration/60); % 1/min
        clustering.qAR(c_n,:)=AR_subfunction(clustering.qIED(c_n,:)); % [0 0.5 1]
        
        REDs=repetitiveIED_v02(EVT.MTABS(CL.class==c_unique(c_n),:),reds_tmin,reds_Tmin);
        clustering.qREDs(c_n,:)=REDs.N/(duration/3600); % 1/h
        
        clustering.qDELAY(c_n,:)=DELAY_subfun(clustering.qAR(c_n,:),EVT.MP(CL.class==c_unique(c_n),:))*1000; % ms
        
        clustering.cluster_num(c_n,:)=[sum(CL.class==c_unique(c_n)) sum(CL.class>0)];
        clustering.cluster_percent(c_n,:)=100*sum(CL.class==c_unique(c_n))/sum(CL.class>0); % procento klastrovaných
        clustering.evt_num(c_n,:)=[sum(CL.class==c_unique(c_n)) length(CL.class)];
        clustering.evt_percent(c_n,:)=100*sum(CL.class==c_unique(c_n))/length(CL.class);
        clustering.evt_ied_MSM(c_n,:)=msm(sum(MQ>0,2));
        clustering.total_number_evt(c_n,:)=size(MQ,1);
        clustering.total_number_ied(c_n,:)=sum(MQ(:)>0);
        
        
        if ~isempty(EVT.MRAW)
            clustering.qAMPL(c_n,:)=median(EVT.MA(CL.class==c_unique(c_n),:),1);
            clustering.qAMPL_raw(c_n,:)=median(EVT.MRAW(CL.class==c_unique(c_n),:),1);
        else
            clustering.qAMPL(c_n,:)=median(EVT.MA(CL.class==c_unique(c_n),:),1);
            clustering.qAMPL_raw(c_n,:)=median(EVT.MA(CL.class==c_unique(c_n),:),1);
        end
    else
        clustering.qIED(c_n,:)=zeros(1,size(EVT.MA,2));     %qIED(c_n,sum(isnan(EVT.MA),1)>0)=NaN;
        clustering.qSINGLE(c_n,:)=zeros(1,size(EVT.MA,2));  %qIED(c_n,sum(isnan(EVT.MA),1)>0)=NaN;
        clustering.qAMPL(c_n,:)=zeros(1,size(EVT.MA,2));
        clustering.qAMPL_raw(c_n,:)=zeros(1,size(EVT.MA,2));    %qIED(c_n,sum(isnan(EVT.MA),1)>0)=NaN;
        clustering.qAR(c_n,:)=zeros(1,size(EVT.MA,2));      %qIED(c_n,sum(isnan(EVT.MA),1)>0)=NaN;
        clustering.qREDs(c_n,:)=zeros(1,size(EVT.MA,2));    %qIED(c_n,sum(isnan(EVT.MA),1)>0)=NaN;
        clustering.qDELAY(c_n,:)=zeros(1,size(EVT.MA,2));   %qIED(c_n,sum(isnan(EVT.MA),1)>0)=NaN;
        
        clustering.cluster_num(c_n,:)=[0 0];
        clustering.cluster_percent(c_n,:)=0;
        clustering.evt_num(c_n,:)=[0 0];
        clustering.evt_percent(c_n,:)=0;
        clustering.evt_ied_MSM(c_n,:)=0;
        clustering.total_number_evt(c_n,:)=0;
        clustering.total_number_ied(c_n,:)=0;
    end
end


clustering.total_time=duration;
clustering.class=CL.class;
clustering.cpca=pca_centering;


end
