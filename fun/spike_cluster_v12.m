function cluster=spike_cluster_v12(MA,MW,pca_centering)

% Clustering algorithm based on principal component analysis (PCA) finds patterns in dataset. 
% Realizations with simmilar component are assigned to same class. Non-significant points has zero class.
% 
% 
% INPUT:
%     MA/MW... double amplitude/full_weight metrix [n x dim], n-number of realization, dim-dimension
%     MW.*(MV>0)... double condition metrix [n x dim], n-number of realization, dim-dimension
%     pca_centering... set centering of dataset MA: MA=MA-repmat(pca_centering*mean(MA,1),size(MA,1),1);
%                  ... 0 - no centered, higher merging/ lower separation
%                  ... 1 - full centered
%         
% OUTPUTS:
%     cluster.class... matrix [n x 1], cluster assigment (>0). Zero value marks unclassified realizations
%     cluster.weight ... matrix [n x 1], percentage of classes
%     cluster.area ... matrix [n x 1], number of realization in same classes
%
% Made by ISARG 2014, Radek Janca (jancarad@fel.cvut.cz)
%
% rozpracov�no !!!!!!!!!!!!!!!!!!!!!!!!!


if ~exist('pca_centering','var')
    pca_centering=0;
end

% init
idx_ch=1:size(MA,2);
cluster.class=zeros(size(MA,1),1);
cluster.weight=zeros(size(MA,1),1);
cluster.area=zeros(size(MA,1),1);
% z-score norm.

max_events=0.25e6;



fprintf(1,'events passed ')
stack=true;
while stack;
    
    zero_length=sum(cluster.class==0);
    
    fprintf(1,[num2str(100*(1-zero_length/size(MA,1)),'% 04.0f') '%%'])
    
    if exist('MW')>0
        warning off all
        
        qEEG_MW=sum(MW(cluster.class==0,:),1);
        try
            [cm,centr]=kmeans(qEEG_MW,2,'replicates',20);
            [~,poz]=max(centr);
            significant_ch=find(cm==poz);
        
        catch
            try
                [cm,centr]=kmeans(qEEG_MW,2,'start',[min(qEEG_MW) ;max(qEEG_MW)]);
                [~,poz]=max(centr);
                significant_ch=find(cm==poz);
            catch
                significant_ch=1:size(MW,2);
            end
        end
        
         
        idx_ch=[];
        
        [~,IX_ch]=sort(qEEG_MW,'descend');
        try
            idx_ch=sort(IX_ch(1:2*length(significant_ch)));
        catch
            idx_ch=1:size(MW,2);
        end
    end
    
    class_itt=max(cluster.class);
    
    idx_evt=cluster.class==0;
%     disp([num2str(sum(idx_evt)) ' events'])
    
    if sum(idx_evt)>max_events
        poz=find(idx_evt);
        idx_evt=false(size(idx_evt));
        idx_evt(poz(randperm(length(poz),round(0.75*max_events))))=true; 
    end

    
    [cluster.class(idx_evt),~,cluster.area(idx_evt)]=spike_cluster_sub(MA(idx_evt,idx_ch),pca_centering);
    if (zero_length-sum(cluster.class==0))==0 
        [cluster.class(idx_evt),~,cluster.area(idx_evt)]=spike_cluster_sub(MA(idx_evt,:),pca_centering);
    end
    cluster.class(idx_evt & cluster.class>0)=cluster.class(idx_evt & cluster.class>0)+class_itt;
    
    if sum(cluster.class==0)==0 || (zero_length-sum(cluster.class==0))==0 
       stack=0;
    end
    
    fprintf(1,repmat('\b',1,4))
end
fprintf(1,' done\n')


cluster=cluster_merging(cluster,MW);

end







function [class,w,area]=spike_cluster_sub(data,pca_centering)


if size(data,1)==1 % only one realization
   class=1;
   w=100;
   area=1;
   return
end

ch_err=(sum(isnan(data),1)==0); % removing NaN values
data=data(:,ch_err);


data=data-repmat(pca_centering*mean(data,1),size(data,1),1); % CENTERING

econ=true;
opt = statset('pca'); opt.MaxIter = 200;

try
    [~,D_pca,eigval] = princomp(data','econ'); % PCA
catch
    try
        [~,D_pca,eigval] = pca(data','Centered',false,'Economy',econ,'Algorithm','svd','Options',opt);
    catch
        [~,D_pca,eigval] = pca(data','Centered',false,'Economy',econ,'Algorithm','eig','Options',opt);
    end
%     % requires PCA toolbox for MATLAB 1.3 (Milano) 
%     model = pca_model(data',size(data,2),'none');
%     eigval=model.E;
%     D_pca=model.T;
    
end


parfor i=1:100
%   disp(['PCA threshold: ' num2str(i) '/100'])
    rand_eigval=[];
    if size(data,1)>0.5e5
        sub_data=data(randperm(size(data,1),0.5e5),:);
        sub_data=sub_data(reshape(randperm(numel(sub_data)),size(sub_data)))';
        try
            [~,~,rand_eigval] = princomp(sub_data,'econ');
        catch
            [~,~,rand_eigval] = pca(sub_data,'Centered',false,'Economy',econ,'Algorithm','svd','Options',opt);
        end
    else
        try
            [~,~,rand_eigval] = princomp(data(reshape(randperm(numel(data)),size(data)))','econ'); % PCA
        catch
            [~,~,rand_eigval] = pca(data(reshape(randperm(numel(data)),size(data)))','Centered',false,'Economy',econ,'Algorithm','svd','Options',opt); % PCA
        end
    end
    re(i)=max(rand_eigval/sum(rand_eigval));
end


% eigval=eigval/sum(eigval); % eigenvalue normalization to 100%
% n_cluster=sum(eigval>mean(re)); % significant eigenvalue treshold

% n_cluster=sum((eigval/sum(eigval))>mean(re)); % significant eigenvalue treshold
n_cluster=sum((eigval/sum(eigval))>quantile(re,0.95)); % significant eigenvalue treshold

C=corr(data',D_pca); % eigenvectors correlation to data

% positive
[m1,c1]=max(C,[],2); % maximal correlation of each positive eigenvectors

% negative
[m2,c2]=max(-C,[],2); % maximal correlation of each neagative eigenvectors

if sum(m1(c1<=n_cluster))>sum(m2(c2<n_cluster)) % use sign of better correlated eigenvectors
    c=c1;
else
    c=c2;
end
% ------------------


% cluster info
c(c>n_cluster)=0;
class=zeros(size(c));
area=zeros(size(c));


if n_cluster==0
    area=size(data,1);
    w=100*ones(size(class));
else
    
    for i=1:n_cluster
        W(i)=sum(c==i);
    end
    [~,w_idx]=sort(W,'descend');
    
    
    for i=1:n_cluster
        class(c==w_idx(i))=i;
        area(class==i)=sum(class==i);
    end
    area(class==0)=sum(class==0);
    w=100*area/size(data,1);
end
end



function cluster=cluster_merging(cluster,MW)

while true
    if sum(cluster.class>0)==0; break; end
    
    qeeg=[];
    for i=1:max(cluster.class);
        qeeg(:,i)=sum(MW(cluster.class==i,:),1)';
    end
    
    C=corr(qeeg);
    C(logical(eye(size(C))))=-Inf;
    [row,col]=find(C>0.9);
    
    if isempty(row); break; end
%     disp('correlation merging')
    cluster.class(cluster.class==col(1))=row(1);
    
    c_unique=unique(cluster.class);
    c_unique(c_unique==0)=[];
    
    class=zeros(size(cluster.class));
    area=zeros(size(cluster.area));
    for i=1:length(c_unique)
        class(cluster.class==c_unique(i))=i;
        area(cluster.class==c_unique(i))=sum(cluster.class==c_unique(i));
    end
    cluster.class=class;
    cluster.area=area;
end


cluster.weight=100*cluster.area/size(MW,1);

[c_uniq,c_idx]=unique(cluster.class);
c_weight=cluster.weight(c_idx(c_uniq>0));
c_uniq(c_uniq==0)=[];

[~,s_idx]=sort(c_weight,'descend');

CLS=cluster.class;
for i=1:length(s_idx)
    cluster.class(CLS==c_uniq(s_idx(i)))=i;
end

end
