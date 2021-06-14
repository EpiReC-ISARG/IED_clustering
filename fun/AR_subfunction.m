function qAR=AR_subfunction(qEEG)

try
    [cm,centr]=kmeans(qEEG,2,'replicates',20);
    [~,poz]=max(centr);
    qAR=0.5*double(cm==poz);
catch
    try
        [cm,centr]=kmeans(qEEG,2,'start',[min(qEEG) ;max(qEEG)]);
        [~,poz]=max(centr);
        qAR=0.5*double(cm==poz);
    catch
        qAR=0.5*ones(1,size(qEEG,2));
    end
end

qAR(qEEG>0.5*max(qEEG))=1;
qAR=qAR(:)';
end