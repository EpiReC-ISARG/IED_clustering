function MAIN_fun_standalone(inDATA)

if nargin==0
    [file,path] = uigetfile('*.mat');
    data=load([path file]);
else
    data=load(inDATA);
end



[d,labels]=ref2bip_v4(double(data.d),data.labels);

[~, discharges]=spike_detector_hilbert_v24(d,data.fs);
discharges.MTABS=datenum(0,0,0,0,0,discharges.MP)+data.tabs(1);

%% klastrování
duration=(size(d,1)/data.fs/60);
clustering=clustering_main_fun(discharges,duration);

%% figures
figure(1); clf
cidx=find(clustering.evt_percent>5);
n_line=ceil(length(cidx)/3);

% IED
qEEG=sum(discharges.MW.*(discharges.MV>0)/duration);
[idx,C]=kmeans(qEEG(:),2,'replicates',10);
[~,cpoz]=max(C);
subplot(n_line+1,3,[1 3])
bar(qEEG);hold on;
set(gca,'XTick',1:length(qEEG),'XTickLabel',labels); xtickangle(45);
bar(find(idx==cpoz),qEEG(idx==cpoz),'r')
title('total IED rate'); ylabel('IED/min')



% clusters
for i=1:length(cidx)
    subplot(n_line+1,3,3+i)
    cied=clustering.qIED(cidx(i),:);
    bar(cied);hold on; 
    ylim([0 max(qEEG)])
    set(gca,'XTickLabel',labels(1:2:end)); xtickangle(45);
    [idx,C]=kmeans(cied(:),2,'replicates',10);
    [~,cpoz]=max(C);
    bar(find(idx==cpoz),cied(idx==cpoz),'r')
    title(['Cluster #' num2str(i) ' - ' num2str(clustering.evt_percent(i),'%.01f') '%'])
end


%% originators
[IEDmax,IEDpoz]=max(clustering.qIED(cidx),[],2);
origin=unique(IEDpoz);

origin_weight=[];
for i=1:length(origin)
    origin_weight(i,1)=sum(IEDmax(IEDpoz==origin(i)));
end

figure(1)
subplot(n_line+1,3,[1 3])
plot(origin,origin_weight,'ok','MarkerFaceColor','c');
% set(gca,'XTick',origin,'XTickLabel',labels(origin)); xtickangle(45); xlim([0 length(qEEG)+1])
legend('background','irritative zone','origin')


tab=cell(size(clustering.qIED,1)+2,size(clustering.qIED,2)+1);
tab(1,:)=[{'Bip. channels'},labels(:,1)'];
tab(2,:)=[{'toal IED'},num2cell(qEEG)];
for i=1:size(clustering.qIED,1)
    tab(2+i,:)=[{['cluster #' num2str(i)]},num2cell(clustering.qIED(i,:))];    
end

writetable(cell2table(tab),'results.csv','Delimiter',';','WriteVariableNames',0)
