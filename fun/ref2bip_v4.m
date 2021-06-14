function [d,C_label,mapa_bip,Num_label]=ref2bip_v4(data,source_file,mapa)

if nargin>2
    if isstruct(mapa)
        mapa=mapa.mapa;
    end
end

if iscell(source_file)
    C_load{1,1}=source_file;
else
    fid=fopen(source_file,'r');
    C_load = textscan(fid, '%s','Delimiter',',');
    fclose(fid);
end

% odstranìní " ' "
idx=cellfun(@(x) strfind(x,char(39)),C_load{1},'UniformOutput',0);
C_load{1}=cellfun(@(x,y) x(setxor(1:length(x),y)),C_load{1},idx,'UniformOutput',0);

idx=cellfun(@(x) ~isempty(intersect(char(x),[0:47 58:64 91:96 123:255])),C_load{1});
if sum(idx)>0
    error(['Not supported char (0-9,a-Z) in name: "' C_load{1}{idx} '"'])
end

prf_idx=cellfun(@(x) find((x>=65 & x<=90) | (x>=97 & x<=122),1,'last'),C_load{1});
prf=cellfun(@(x,y) x(1:y),C_load{1},num2cell(prf_idx),'UniformOutput',0);
sfx=cellfun(@(x,y) str2double(x(y+1:end)),C_load{1},num2cell(prf_idx),'UniformOutput',1);

% unikátní prefixy
[unique_prf,idx]=unique(prf);
[~,idx]=sort(idx);
unique_prf=unique_prf(idx);

unused_el=false(size(C_load{1},1),1);
bip_list=[];
for i=1:length(unique_prf)
    idx=cellfun(@(x) strcmp(x,unique_prf(i)),prf);
    if sum(idx)==1
        unused_el(idx)=true;
    end
    
    idx=find(idx);
    [~,sfx_idx]=sort(sfx(idx));
    idx=idx(sfx_idx);
    
    bip_list=[bip_list;[idx(1:end-1),idx(2:end)]];
end

C_label=cell(size(bip_list,1),2);
Num_label=cell(size(bip_list,1),1);
for i=1:size(bip_list,1)
    C_label(i,1)={[prf{bip_list(i,1)} num2str(sfx(bip_list(i,1))) '-' num2str(sfx(bip_list(i,2)))]};
    C_label(i,2)={bip_list(i,:)};
    Num_label(i,1)={[num2str(sfx(bip_list(i,1))) '-' num2str(sfx(bip_list(i,2)))]};
end



% bipolární data ------------------------------------------------------
if ~isempty(data)
    d=zeros(size(data,1),size(bip_list,1));
    for i=1:size(bip_list,1)
        d(:,i)=data(:,bip_list(i,1))-data(:,bip_list(i,2));
    end
else
    d=[];
end

% bipolární mapa
if nargin>2
    mapa_bip=zeros(size(mapa));
    for i=1:size(bip_list,1)
        [r1,c1]=find(mapa==bip_list(i,1));
        [r2,c2]=find(mapa==bip_list(i,2));
        
        mapa_bip(floor(mean([r1 r2])),floor(mean([c1 c2])))=i;
    end
else
    mapa_bip=[];
end

