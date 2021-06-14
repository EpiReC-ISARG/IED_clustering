function IMGD=dilateN(IMG,MASK)
% IMG ... n-dim binary image (logical)
% MASK ... dilate n-dim mask (logical) or size of dilatation mask (pixels)
%
% Example:
% IMGD=dilateN(IMG,10) extends boundaries of 5 pixels 

if ~islogical(IMG)
   error('only for logical image') 
end    

if numel(MASK)==1
    MASK=single(ones(MASK*ones(1,length(size(IMG)))));
end

IMG=single(IMG);
MASK=single(MASK);

try
    if gpuDeviceCount>0
%         disp('convn on GPU')
        IMG=gpuArray(IMG);
        MASK=gpuArray(MASK);
        IMGD=convn(IMG,MASK,'same');
        IMGD=gather(IMGD)>0;
    else
        IMGD=convn(IMG,MASK,'same')>0;
    end
catch
    IMGD=convn(IMG,MASK,'same')>0;
end
