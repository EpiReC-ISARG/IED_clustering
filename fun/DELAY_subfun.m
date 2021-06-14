
function qDELAY=DELAY_subfun(qAR,MP)

MP=MP-repmat(min(MP,[],2),1,size(MP,2));
MP(:,qAR==0)=NaN;
MP(sum(~isnan(MP),2)<2,:)=[];

% qDELAY=trimmean(MP,25,'weighted',1);
Q95=quantile(MP,0.95,1);
MP(MP>repmat(Q95,size(MP,1),1))=NaN;
MP(sum(~isnan(MP),2)<2,:)=[];

qDELAY=nanmean(MP,1);

end