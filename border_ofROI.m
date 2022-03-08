%this function outputs the region surrounding a ROI, the values of pxls in
%each of surrounding borders, median and the average of these pixels

function [border, timeaverage_out, timemedian_out, pxls_traces_out]= shehab_border_ofROI(mask_x, Y)
newmask=zeros(size(Y,1),size(Y,2));
thinout= bwboundaries(mask_x,'noholes');
aa=thinout{1}; clear thinout
linearInd = sub2ind([size(Y,1),size(Y,2)], aa(:,1), aa(:,2));clear aa
newmask(linearInd)=1; clear aa , %figure, imagesc(newmask)

newermask=zeros(size(Y,1),size(Y,2));
thickout=boundarymask(newmask); %figure, imagesc(thickout)
thinthickout=bwboundaries(thickout,'noholes');
aa=thinthickout{1}; clear thinthickout 
linearInd = sub2ind([size(Y,1),size(Y,2)], aa(:,1), aa(:,2)); clear aa
newermask(linearInd)=1; % figure, imagesc(newermask)

thicknewermask=boundarymask(newermask);
%figure,imagesc(thicknewermask)

newestmask=zeros(size(Y,1),size(Y,2));
thinnewest=bwboundaries(thicknewermask,'noholes');
aa=thinnewest{1}; clear thinnewest
linearInd = sub2ind([size(Y,1),size(Y,2)], aa(:,1), aa(:,2)); clear aa
newestmask(linearInd)=1; % figure, imagesc(newestmask)

border=boundarymask(newestmask);
%figure, imagesc(border)


[Index]=find(border==1);
       for j=1:size(Y,3)               % time dimension
           a=Y(:,:,j);                    % take a certain frame out of the whole stack
           pic(:,j)=a(Index);             % takes the pixels that passes through the mask from the chosen frame
       end , clear j
       pxls_traces_out=pic;, clear pic a Index;
       timeaverage_out=mean(pxls_traces_out,1);   
       timemedian_out=median(pxls_traces_out,1);   
       
       
       
