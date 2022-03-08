%% Directory for tif files, ROIs (regions of interest) produced by imageJ and directory for saving output which is traces of calcium fluorescent signal
%of all ROIs 

videos={ '/media/file_1.tif', ...
         '/media/file_2.tif', ...
       };
 
ROIs={  '/media/denois_RoiSet_1/', ...
        '/media/denois_RoiSet_2/', ...
     };

Savedir={  '/media/Cleaned_traces_1', ...
           '/media/Cleaned_traces_2', ...
        };          

%%
for file=1:length(videos);    
    
% loading videos and ROIs
fr= 4;              %%frames per second of the recording%%
nam = videos{file};
fileinfo=imfinfo(nam);
for i=1:length(fileinfo)
    Y(:,:,i:i)=double(imread(nam,i));
end
clear fileinfo i 
if ~isa(Y,'double');    Y = double(Y);  end         % convert to single
[d1,d2,T] = size(Y);                                % dimensions of dataset
d = d1*d2;
clear nam T

%%%% Extracting ROIs (masks) %%%%
directory=dir(ROIs{file});
numOFrois=length(directory)-2;
cellOFrois=cell(numOFrois,1);
for i=1:numOFrois
    [sROI] = ReadImageJROI([ROIs{file} '/' num2str(directory(i+2).name)]); 
    cellOFrois{i}=sROI.mnCoordinates;
end

resiz_factor= 1; % Insert the resizing factor of Y (the image in the main script)
test(1:d,numOFrois)=zeros;                               

mask_=[];
for ii=1:numOFrois
    mask=poly2mask(round(cellOFrois{ii}(:,1)),round(cellOFrois{ii}(:,2)),d1/resiz_factor,d2/resiz_factor);
    masksmall=imresize(mask,resiz_factor);
    mask_(:,:,ii)=masksmall;
clear mask onedimMASK onedimINTENSITY
end
clear d d1 d2 i ii masksmall numOFrois sROI cellOFrois directory test resiz_factor

%% BG deduction
tic
interv_length=0.25*size(Y,3);
%determine which ROI represents the BG, based on its size
for i=1:size(mask_,3)                                                              
    sizes(i)=length(find(mask_(:,:,i)));                                            % find the roi with the largest size
end , [~,liveBG] = max(sizes); clear i sizes                                        % take its index 

%% determine whether background is the whole field (including no tissue) or only the live tissue
where_to_get_BG='only_live_tissue';
%where_to_get_BG='the_whole_field'

switch where_to_get_BG
      case 'only_live_tissue'
            % build the SD projection
            for x=1:size(Y,1)
                for y=1:size(Y,2)
                    if   mask_(x,y,liveBG)==1                                       % if you are inside the ROI, do the following    
                         SD_proj(x,y)=std(Y(x,y,:));                                % calculate the STD through time (3rd dimension of Y)     
                    else SD_proj(x,y)=0;                                            % else just label the pixel as 0
                    end
                end
            end , clear x y 

            %detect the value of the 25% index in a sorted array of the SD projection
            oneD_SD_BG=reshape(SD_proj,size(SD_proj,1)*size(SD_proj,2),1);          % convert the matrix to a one-dimension vector 
            sort_oneD_SD_BG=sort(oneD_SD_BG,'ascend');                              % sort the values in an ascending order
            strt_ind=find(sort_oneD_SD_BG);                                         % take the first index for non-0 value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHOSE THE PERCENTAGE THRESHOLD FOR THE BG NOISE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            IND=round(0.05*length(strt_ind))+strt_ind(1);                           % take the index of the 25% quartile (since the vector is full of 0s, index of first non-0 is added)
            thresh=sort_oneD_SD_BG(IND);                                            % take the value of this index

            %build SD_projection MASK
            for x=1:size(Y,1)
                for y=1:size(Y,2)
                    if  SD_proj(x,y)<thresh & SD_proj(x,y)~=0                       % take the pixels with SD below the threshold but also not equal to 0, to represent the BG noise
                        SD_BG(x,y)=1;                                               % label the pixels below the threshold SD as 1
                    else
                        SD_BG(x,y)=0;                                               % label the pixels above the threshold SD as 0
                    end
                end
            end , clear x y

            % apply the mask on the video
            [Index]=find(SD_BG==1);                                                 % create a mask for the BG pixels
            for j=1:size(Y,3)                                                       
                a=Y(:,:,j);     
                pic(:,j)=a(Index);
            end 
            BG_pxls_traces=pic; low_pxl_value=min(min(pic)); clear pic a j Index ;  % calculate also the lowest value inside the ROI representing BG , to be added later as an offset to avoid negative values
            BG_timeaverage=mean(BG_pxls_traces,1);                                  % produce a trace based on the average of the BG activity 

                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      case 'the_whole_field';   
            % build the SD projection
            SD_proj=std(Y,0,3);                                                     % produce the STD-projection                                                     

            %detect the value of the 25% index in a sorted array of the SD projection
            oneD_SD_BG=reshape(SD_proj,size(SD_proj,1)*size(SD_proj,2),1);          % convert the matrix to a one-dimension vector 
            sort_oneD_SD_BG=sort(oneD_SD_BG,'ascend');                              % sort the values in an ascending order
            IND=round(0.25*length(sort_oneD_SD_BG));                                % take the index of the 25% quartile
            thresh=sort_oneD_SD_BG(IND);                                            % take the value of this index

            %build SD_projection MASK
            for x=1:size(Y,1)
                for y=1:size(Y,2)
                    if  SD_proj(x,y)<thresh;                                        % take the pixels with SD below the threshold, to represent the BG later
                        SD_BG(x,y)=1;                                               % label them as 1
                    else
                        SD_BG(x,y)=0;           bruss                                    % else label them as 0
                    end
                end
            end , clear x y

            % apply the mask on the video
            [Index]=find(SD_BG==1);                                                 % create a mask for the BG pixels
            for j=1:size(Y,3)                                   
                a=Y(:,:,j);     
                pic(:,j)=a(Index);
            end 
            BG_pxls_traces=pic; low_pxl_value=min(min(pic)); clear pic a j Index ;  % calculate also the lowest value inside the ROI representing BG , to be added later as an offset to avoid negative values
            BG_timeaverage=mean(BG_pxls_traces,1);                                  % produce a trace based on the average of the BG activity 
end

% create best fit for the BG trace over time
     BG_yfit=fit_function(BG_timeaverage,interv_length);                            % taking the best fit to be used as baseline

% subtract BG fit from each pixel
old_Y=Y;
for x=1:size(Y,1)                                                                   % let it run for all pixels in X-dim
    for y=1:size(Y,2)                                                               % let it run for all pixels in Y-dim 
        Y(x,y,:)=squeeze(old_Y(x,y,:))-BG_yfit+low_pxl_value;                      % deduct the BG_yfit from each pixel's trace + add the lowest pixel value (within the live material), to avoid negative values
    end
end, clear x y
disp('Background is denoised');
toc

clear BG_pxls_traces BG_timeaverage BG_yfit IND liveBG low_pxl_value oneD_SD_BG SD_BG SD_proj sort_oneD_SD_BG strt_ind thresh where_to_get_BG

%% Calculating diff. variables
%%% In this section, variables are extracted from raw data, and a guideMAT is created, including when do POTENTIAL events (true or false) ROUGHLY
%%% take place (start and end points are not precise). Then a newguideMAT is created with the POTENTIAL events sorted as TRUE (1) or FALSE (2) events.
%%%%%%%%%%%%%%%%%%% Needed variables from previous section : 1) Y 2)mask_ 3)fr %%%%%%%%%%%%%%%%%%% 
checker=[];
checker_indMAT=[];

tic
interv_length=0.25*size(Y,3);                                                                          % choose length of interval used for baseline fitting

newguideMAT(1:size(mask_,3),1:size(Y,3))=zeros;                                                        % preallocate 'newguideMAT' to avoid 

for roi=1:size(mask_,3)
        [~, timeaverage_out, timemedian_out, pxls_traces_out]= border_ofROI(mask_(:,:,roi), Y); % this is a function that provides info (pixels and values of pxls) about the borders surrounding each ROI
        [Index]=find(mask_(:,:,roi)==1);
        for j=1:size(Y,3)
            a=Y(:,:,j);
            pic(:,j)=a(Index);
        end , clear j
        pxls_traces=pic; clear pic a Index;                                                            % the pxls' values of each ROI are produced here
        
        timeaverage=mean(pxls_traces,1);   globaltimeaverage(roi:roi,1:length(timeaverage))=timeaverage;   % timeaverage for each ROI is produced and stored in globaltimeaverage according to its order among other ROIs
                                           globaltimeaverage_out(roi,1:length(timeaverage_out))=timeaverage_out;     
        timemedian =median(pxls_traces,1); globaltimemedian (roi,1:length(timemedian)) =timemedian;    % timemedian for each ROI is produced and stored in globaltimemedian according to its order among other ROIs
        %generalstd=std(pxls_traces,1);                                                                 % calculating the std of the values of pxls within the examined ROI

%%%%THRESHOLD for number of clusters%%%%       
%         for i=1:length(generalstd)                                                  % forloop for setting the number of clusters (no_of_k) based on the std (generalstd)
%             if      generalstd(i)>=prctile(generalstd,99.7);no_of_k(roi,i)=3;           % threshold is 2 STD points (99.7%)
%             elseif  generalstd(i)>=prctile(generalstd,95)  ;no_of_k(roi,i)=2;           % threshold is 1 STD points (95%)
%             else    no_of_k(roi,i)=1;
%             end
%         end, clear i

%%%%THRESHOLD for baseline%%%%
        %time10prctile=prctile(pxls_traces,10,1);                                   % the value of certain percentile (whatever is decided) is calculated here 
        %yfit=expofit4CA(1:length(time10prctile),time10prctile);                    % taking the best fit to be used as a baseline  (exponential fitting)
        %yfit(1:length(globaltimeaverage(roi,:)),1)=median(globaltimeaverage(roi,:));  % yfit based on median, this fits (probably) when video is denoised from before (using gaussian blur in FIJI)  
        yfit=fit_function(globaltimeaverage(roi,:),interv_length);                % taking the best fit to be used as a baseline  (alternative method based on linear fitting between chosen points)
        globalyfit(roi,1:length(yfit))=yfit;                                        % saving the best fit of each ROI in globalyfit according to its order among other ROIs
        [~, d1ydx, ~, d2ydx]=first_second_derivative(1, timeaverage);               % calculates first & second derivative 
        globald1ydx(roi,1:length(d1ydx))=d1ydx;                                     % the first derivative for all ROIs
        global_low_d1ydx(roi,1:length(d1ydx))=butterlowpass(d1ydx,fr,2,8);          % low-pass-filtered first derivative for all ROIs
        globald2ydx(roi,1:length(d2ydx))=d2ydx;                                     % the second derivative for all ROIs

%%%%THRESHOLD defining what is an event, based on the prctile of pxls crossing the baseline%%%%        
        guideMAT=(prctile(pxls_traces,5,1)>yfit');                                  % a binary matrix with ones at time points when the lowest (whatever percentile) is higher than the baseline
        indMAT(:,1)=find(diff(guideMAT)==1);                                        % find the ROUGHLY STARTing point for potential events (time points when the lowest (whatever percentile) is higher than the baseline)

%%%%%%% XX is to deal with any ROI that is active from the first frame
            XX=find(diff(guideMAT)==-1);
            if   guideMAT(1,1)==1                                                   % **** if the trace starts with an active cell 
                 if guideMAT(1,end)==1, indMAT(:,2)=[ XX(2:end) length(guideMAT)]; else indMAT(:,2)=XX(2:end); end               % **** find the ENDing point for potential events (time points when the lowest (whatever percentile) is higher than the baseline) 
                 indMAT(2:size(indMAT,1)+1,:)=indMAT; indMAT(1,1)=0;indMAT(1,2)=XX(1);   % in case there is event in the begining of the video AND YOU WANT TO KEEP IT
            else if guideMAT(1,end)==1, indMAT(:,2)=[ XX(1:end) length(guideMAT)]; else indMAT(:,2)=find(diff(guideMAT)==-1); end    % find the ROUGHLY ENDing point for potential events (time points when the lowest (whatever percentile) is higher than the baseline)
            end, clear XX;
        indMAT(:,1)=indMAT(:,1)+1;                                                  % this is to compensate for the shifting effect due to usage of 'diff' function on the starting points of potential events  
        
        checker(roi,1:length(guideMAT))=guideMAT;
        checker_indMAT(roi,1:length(indMAT))=indMAT(:,1);
        %%% Creating a preliminary 'newguidMAT' (including time points; when does ROUGHLY an event start and end) %%%
        %newguideMAT=zeros(size(guideMAT,1),size(guideMAT,2)); % this seems to be bad coding, because this deletes all previously calculated ROIs
        
        for noOFevents=1:size(indMAT,1)                                             % number of potential events  
           %if  indMAT(noOFevents,2)-indMAT(noOFevents,1) >=  2             %%%%%%%%% arbitrary number of frames set to be the minimum temporal size for an event  ###### for low frame rate this is problem ########
               %bef_last_checker{noOFevents,roi}=1:size(indMAT,1);
                for frame=indMAT(noOFevents,1):indMAT(noOFevents,2)

%%%%THRESHOLD defining what is a TRUE event based on the average (or median) value of the border crossing certain prctile of pxls within a ROI%%%%  
                last_checker{noOFevents,roi}=indMAT(noOFevents,1):indMAT(noOFevents,2);
                   if  timeaverage_out(frame) < prctile(pxls_traces(:,frame),5) 
                       newguideMAT(roi,frame)=1;                                    % the average or median of the border is NOT MORE than the set (previously) prctile of the ROI's pxls' values
                   else
                       newguideMAT(roi,frame)=2;                                    % the average or median of the border is MORE than the set (previously) prctile of the ROI's pxls' values
                   end, clear frame
               end
           %end             %%%%%%%%% ###### for low frame rate this is problem ######## %%%%%%%%%
        end
        %newguideMAT(roi,indMAT(end,2)+1:length(timeaverage))=0;                     % this is to fill the remaining time points with zeros instead of having short matrices ending at the last event
        
        clear guideMAT  indMAT pxls_traces %timeaverage timemedian time10prctile yfit guideMAT initialMAT timeaverage_out timemedian_out pxls_traces_out noOFevents
        disp(['ROI ' num2str(roi) ' is done ;)' ]);

end
toc
%%
if size(newguideMAT,2) >= size(globald2ydx,2), newguideMAT(:,size(globald2ydx,2):end)=[]; end % this is to avoid having a matrix longer than the second derivative matrix ('globald2ydx')




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Need to check the remaining one-frame events remainig as an artifact at the end of this section

%% Re-adjusting the starting and ending points of TRUE and FALSE events in 'newguideMAT'
%%% In this section, time points of the assorted TRUE and FALSE events are adjusted precisely; an event starts when there is a peak in the 2nd derivative,
%%% and ends when there is a peak in the low-pass-filtered 1st derivative. Events that last for a duration of one frame are cancelled out.
%%%%%%%%%%%%%%%%%%% Needed variables from previous section : 1)newguideMAT 2)global_low_d1ydx 3)globald2ydx %%%%%%%%%%%%%%%%%%% 

 which_should_be_adjusted_first = 'True_bef_False';                                              % choose which type of events should be adjusted first. Based on this, one of the two following cases will be performed
%which_should_be_adjusted_first = 'False_bef_True';

switch which_should_be_adjusted_first
case 'True_bef_False'
        tic , endings_of1=[]; endings_of2=[]; newerguideMAT=newguideMAT;                        % 'newerguideMAT' is made to make checking changes at each step easier to see (compared to 'newguideMAT')
        %%%% adjusting the TRUE events' starting and ending points%%%%
        for roi=1:size(newerguideMAT,1)
          for frame= 2:size(newerguideMAT,2)-1

% % % % % % %                     % clearing any single points with a value different from both the preceding and the following points
% % % % % % %                     if  newerguideMAT(roi, frame-1)==newerguideMAT(roi, frame+1)                % if the value at the previous and the following frames are equal, 
% % % % % % %                         newerguideMAT(roi, frame)=newerguideMAT(roi, frame-1);                  % then the middle value should be equal to them as well.
% % % % % % %                     end

                            % adjust the starting point (for TRUE EVENTS) to the peak of the 2nd deriv.
                            if  (newerguideMAT(roi,frame)==1 && newerguideMAT(roi, frame-1)==0) || (newerguideMAT(roi, frame)==1 && newerguideMAT(roi, frame-1)==2) % this is the point at which an event begins
                                if  frame<11;                                                   % in case the event starts WITHIN the first 10 frames
                                    [~,peak]=max(globald2ydx(roi,1:frame));                     % 'peak' is the point at which the 2nd deriv. is max,, which should be the point at which the event starts
                                    newerguideMAT(roi,1+peak:frame)=1;                          % set the starting point according to the peak of the 2nd deriv.
                                else                                                            % in case the event starts AFTER the first 10 frames
                                    [~,peak]=max(globald2ydx(roi,frame-10:frame));              % 'peak' is the point at which the 2nd deriv. is max,, which should be the point at which the event starts
                                    newerguideMAT(roi,frame-10+peak:frame)=1;                   % set the starting point according to the peak of the 2nd deriv.
                                end
                            end, clear peak

                    % adjust the ending point (for TRUE EVENTS) to the first peak of low-passed 1st deriv after the trough
                    if  (newerguideMAT(roi,frame)==0 && newerguideMAT(roi, frame-1)==1) || (newerguideMAT(roi, frame)==2 && newerguideMAT(roi, frame-1)==1) % this is the point at which an event ends
                         endings_of1(roi,frame)=1;      
                    end, %clear peaks
          end, %clear frame
                
                %Setting the end point based on what occurs first; either the 1st positive value in the 1st derivative or the 1st peak in the 1st derivativ
                [~,Ind]=find(endings_of1);
                for i=1:length(Ind)  
                            [~,pos_value]=find(global_low_d1ydx(roi,Ind(i):end)>0);             % detects when is the low-pass-filtered 1st deriv. positive again 
                            if  length(global_low_d1ydx(roi,Ind(i):end)) < 3, peaks=[];         % this is to avoid errors with the next function ('findpeaks') which require more than three points
                            else[~,peaks]=findpeaks(global_low_d1ydx(roi,Ind(i):end));          % 'peak' is the point at which the low-passed 1st deriv after the trough is max,, which should be the point at which the event ends
                            end
                    
                    if      isempty(pos_value) && isempty(peaks);                               % in case both are empty, then fill till the end with what is available
                            newerguideMAT(roi,Ind(i):end)=1;
                    elseif  isempty(peaks) && ~isempty(pos_value);                              % to avoid the comming steps if there were nothing in 'pos_value'         
                            newerguideMAT(roi,Ind(i):Ind(i)+pos_value(1)-1)=1;                  % choose the index of this positive value
                    elseif  isempty(pos_value) && ~isempty(peaks);                              % to avoid the coming steps if there were nothing in 'pos_value'
                            newerguideMAT(roi,Ind(i):Ind(i)+peaks(1)-1)=1;                      % else, set the ending point according to the peak of the low-passed 1st deriv.
                    else
                        if      pos_value(1) < peaks(1)                                             % in case the 1st derivative is positive before a peak exists
                                newerguideMAT(roi,Ind(i):Ind(i)+pos_value(1)-1)=1;                  % choose the index of this positive value
                        else    newerguideMAT(roi,Ind(i):Ind(i)+peaks(1)-1)=1;                      % else, set the ending point according to the peak of the low-passed 1st deriv.
                        end
                    end
                end

        end, %clear roi

        %%%% adjusting the FALSE events' starting and ending points%%%%
        finalguideMAT=newerguideMAT;                                                            % 'finalguideMAT' is made to make checking changes at each step easier to see (compared to 'newerguideMAT')
        for roi=1:size(finalguideMAT,1)
          for frame= 2:size(finalguideMAT,2)-1

                            % adjust the starting point (for FALSE EVENTS) to the peak of the 2nd deriv.
                            if  (finalguideMAT(roi,frame)==2 && finalguideMAT(roi, frame-1)==0) %|| (newguideMAT(roi, frame)==2 && newguideMAT(roi, frame-1)==1)
                                % this is the point at which an event was beginning
                                if  frame<11;                                                   % in case the event starts WITHIN the first 10 frames
                                    [~,peak]=max(globald2ydx(roi,1:frame));                     % 'peak' is the point at which the 2nd deriv. is max,, which should be the point at which the event starts
                                    finalguideMAT(roi,1+peak:frame)=2;                          % set the starting point according to the peak of the 2nd deriv.
                                else                                                            % in case the event starts AFTER the first 10 frames
                                    [~,peak]=max(globald2ydx(roi,frame-10:frame));              % 'peak' is the point at which the 2nd deriv. is max,, which should be the point at which the event starts
                                    finalguideMAT(roi,frame-10+peak:frame)=2;                   % set the starting point according to the peak of the 2nd deriv.
                                end
                            end, clear peak

                    % adjust the ending point (for FALSE EVENTS) to the first peak of peak of low-passed 1st deriv after the trough
                    if  (finalguideMAT(roi,frame)==0 && finalguideMAT(roi, frame-1)==2) %|| (newguideMAT(roi, frame)==1 && newguideMAT(roi, frame-1)==2) % this is the point at which an event ends
                         endings_of2(roi,frame)=1;     
                    end, %clear peaks
          end

                %Setting the end point based on what occurs first; either the 1st positive value in the 1st derivative or the 1st peak in the 1st derivativ
                [~,Ind_2]=find(endings_of2);
                for i=1:length(Ind_2)  
                           [~,pos_value]=find(global_low_d1ydx(roi,Ind_2(i):end)>0);            % detects when is the low-pass-filtered 1st deriv. positive again 
                           if length(global_low_d1ydx(roi,Ind_2(i):end)) < 3 , peaks=[];        % this is to avoid errors with the next function ('findpeaks') which require more than three points
                           else [~,peaks]=findpeaks(global_low_d1ydx(roi,Ind_2(i):end));        % 'peak' is the point at which the low-passed 1st deriv after the trough is max,, which should be the point at which the event ends
                           end
                    
                    if      isempty(pos_value) && isempty(peaks);                               % in case both are empty, then fill till the end with what is available
                            finalguideMAT(roi,Ind_2(i):end)=2;
                    elseif  isempty(peaks) && ~isempty(pos_value);                              % to avoid the coming steps if there were nothing in 'pos_value'         
                            finalguideMAT(roi,Ind_2(i):Ind_2(i)+pos_value(1)-1)=2;              % choose the index of this positive value
                    elseif  isempty(pos_value) && ~isempty(peaks);                                                 % to avoid the coming steps if there were nothing in 'pos_value'       
                            finalguideMAT(roi,Ind_2(i):Ind_2(i)+peaks(1)-1)=2;                  % else, set the ending point according to the peak of the low-passed 1st deriv.
                    else
                        if     pos_value(1) < peaks(1)                                              % in case the 1st derivative is positive before a peak exists
                               finalguideMAT(roi,Ind_2(i):Ind_2(i)+pos_value(1)-1)=2;               % choose the index of this positive value
                        else   finalguideMAT(roi,Ind_2(i):Ind_2(i)+peaks(1)-1)=2;                   % else, set the ending point according to the peak of the low-passed 1st deriv.
                        end
                    end
                end
             disp(['ROI ' num2str(roi) ' is adjusted ;)' ]);   
        end
        toc, disp('You adjusted the True events before False ones, case 1')

    
case 'False_bef_True'
        tic , endings_of1=[]; endings_of2=[]; newerguideMAT=newguideMAT;                        % 'newerguideMAT' is made to make checking changes at each step easier to see (compared to 'newguideMAT')
        %%%% adjusting the FALSE events' starting and ending points%%%%

        for roi= 1:size(newerguideMAT,1)
          for frame= 2:size(newerguideMAT,2)-1

% % % % % % %                     % clearing any single points with a value different from both the preceding and the following points
% % % % % % %                     if  newerguideMAT(roi, frame-1)==newerguideMAT(roi, frame+1)                % if the value at the previous and the following frames are equal, 
% % % % % % %                         newerguideMAT(roi, frame)=newerguideMAT(roi, frame-1);                  % then the middle value should be equal to them as well.
                    end
                            % adjust the starting point (for FALSE EVENTS) to the peak of the 2nd deriv.
                            if  (newerguideMAT(roi,frame)==2 && newerguideMAT(roi, frame-1)==0) %|| (newguideMAT(roi, frame)==2 && newguideMAT(roi, frame-1)==1)
                                % this is the point at which an event was beginning
                                if  frame<11;                                                   % in case the event starts WITHIN the first 10 frames
                                    [~,peak]=max(globald2ydx(roi,1:frame));                     % 'peak' is the point at which the 2nd deriv. is max,, which should be the point at which the event starts
                                    newerguideMAT(roi,1+peak:frame)=2;                          % set the starting point according to the peak of the 2nd deriv.
                                else                                                            % in case the event starts AFTER the first 10 frames
                                    [~,peak]=max(globald2ydx(roi,frame-10:frame));              % 'peak' is the point at which the 2nd deriv. is max,, which should be the point at which the event starts
                                    newerguideMAT(roi,frame-10+peak:frame)=2;                   % set the starting point according to the peak of the 2nd deriv.
                                end
                            end, clear peak

                    % adjust the ending point (for FALSE EVENTS) to the first peak of peak of low-passed 1st deriv after the trough
                    if  (newerguideMAT(roi,frame)==0 && newerguideMAT(roi, frame-1)==2) %|| (newguideMAT(roi, frame)==1 && newguideMAT(roi, frame-1)==2) % this is the point at which an event ends
                         endings_of2(roi,frame)=1;     
                    end, %clear peaks
          end

                %Setting the end point based on what occurs first; either the 1st positive value in the 1st derivative or the 1st peak in the 1st derivativ
                [~,Ind_2]=find(endings_of2);
                for i=1:length(Ind_2)  
                           [~,pos_value]=find(global_low_d1ydx(roi,Ind_2(i):end)>0);            % detects when is the low-pass-filtered 1st deriv. positive again 
                           [~,peaks]=findpeaks(global_low_d1ydx(roi,Ind_2(i):end));             % 'peak' is the point at which the low-passed 1st deriv after the trough is max,, which should be the point at which the event ends
                    
                    if     isempty(pos_value) && isempty(peaks);                                % in case both are empty, then fill till the end with what is available
                           newerguideMAT(roi,Ind(i):end)=2;
                    elseif isempty(peaks) && ~isempty(pos_value);                               % to avoid the coming steps if there were nothing in 'pos_value'
                           newerguideMAT(roi,Ind_2(i):Ind_2(i)+pos_value(1)-1)=2;               % choose the index of this positive value
                    elseif isempty(pos_value) && ~isempty(peaks);                               % to avoid the coming steps if there were nothing in 'pos_value'
                           newerguideMAT(roi,Ind_2(i):Ind_2(i)+peaks(1)-1)=2;                   % else, set the ending point according to the peak of the low-passed 1st deriv.
                    else
                        if     pos_value(1) < peaks(1)                                              % in case the 1st derivative is positive before a peak exists
                               newerguideMAT(roi,Ind_2(i):Ind_2(i)+pos_value(1)-1)=2;               % choose the index of this positive value
                        else   newerguideMAT(roi,Ind_2(i):Ind_2(i)+peaks(1)-1)=2;                   % else, set the ending point according to the peak of the low-passed 1st deriv.
                        end
                    end
                end
        

        finalguideMAT=newerguideMAT;                                                            % 'finalguideMAT' is made to make checking changes at each step easier to see (compared to 'newerguideMAT')
        %%%% adjusting the TRUE events' starting and ending points%%%%
        for roi= 1:size(finalguideMAT,1)
            for frame= 2:size(finalguideMAT,2)-1

                            % adjust the starting point (for TRUE EVENTS) to the peak of the 2nd deriv.
                            if  (finalguideMAT(roi,frame)==1 && finalguideMAT(roi, frame-1)==0) || (finalguideMAT(roi, frame)==1 && finalguideMAT(roi, frame-1)==2) % this is the point at which an event begins
                                if  frame<11;                                                   % in case the event starts WITHIN the first 10 frames
                                    [~,peak]=max(globald2ydx(roi,1:frame));                     % 'peak' is the point at which the 2nd deriv. is max,, which should be the point at which the event starts
                                    finalguideMAT(roi,1+peak:frame)=1;                          % set the starting point according to the peak of the 2nd deriv.
                                else                                                            % in case the event starts AFTER the first 10 frames
                                    [~,peak]=max(globald2ydx(roi,frame-10:frame));              % 'peak' is the point at which the 2nd deriv. is max,, which should be the point at which the event starts
                                    finalguideMAT(roi,frame-10+peak:frame)=1;                   % set the starting point according to the peak of the 2nd deriv.
                                end
                            end, clear peak

                    % adjust the ending point (for TRUE EVENTS) to the first peak of low-passed 1st deriv after the trough
                    if  (finalguideMAT(roi,frame)==0 && finalguideMAT(roi, frame-1)==1) || (finalguideMAT(roi, frame)==2 && finalguideMAT(roi, frame-1)==1) % this is the point at which an event ends
                         endings_of1(roi,frame)=1;      
                    end, %clear peaks
          end, %clear frame

                %Setting the end point based on what occurs first; either the 1st positive value in the 1st derivative or the 1st peak in the 1st derivativ
                [~,Ind]=find(endings_of1);
                for i=1:length(Ind)  
                            [~,pos_value]=find(global_low_d1ydx(roi,Ind(i):end)>0);             % detects when is the low-pass-filtered 1st deriv. positive again 
                            [~,peaks]=findpeaks(global_low_d1ydx(roi,Ind(i):end));              % 'peak' is the point at which the low-passed 1st deriv after the trough is max,, which should be the point at which the event ends
                    
                    if      isempty(pos_value) && isempty(peaks);                               % in case both are empty, then fill till the end with what is available
                            finalguideMAT(roi,Ind(i):end)=1;
                    elseif  isempty(peaks) && ~isempty(pos_value);                              % to avoid the coming steps if there were nothing in 'pos_value'
                            finalguideMAT(roi,Ind(i):Ind(i)+pos_value(1)-1)=1;                  % choose the index of this positive value
                    elseif  isempty(pos_value) && ~isempty(peaks);                              % to avoid the coming steps if there were nothing in 'pos_value'
                            finalguideMAT(roi,Ind(i):Ind(i)+peaks(1)-1)=1;                      % else, set the ending point according to the peak of the low-passed 1st deriv.
                    else
                        if      pos_value(1) < peaks(1)                                             % in case the 1st derivative is positive before a peak exists
                                finalguideMAT(roi,Ind(i):Ind(i)+pos_value(1)-1)=1;                  % choose the index of this positive value
                        else    finalguideMAT(roi,Ind(i):Ind(i)+peaks(1)-1)=1;                      % else, set the ending point according to the peak of the low-passed 1st deriv.
                        end
                    end
                end
             disp(['ROI ' num2str(roi) ' is adjusted ;)' ]); 
        end, %clear roi
        toc, disp('You adjusted the False events before True ones, case 2')
end


%% Cleaning the traces by following the 'finalguideMAT'
%%%%%%%%%%%%%%%%%% Needed variables from previous section : 1)Y  2)finalguideMAT  3)mask_ %%%%%%%%%%%%%%%%%%% 
tic
clean_traces=[];
no_of_k=[];
for roi=1:size(finalguideMAT,1)
        
    % determine the pxls_traces for the current ROI 
        [Index]=find(mask_(:,:,roi)==1);
        for j=1:size(Y,3)
            a=Y(:,:,j);
            pic(:,j)=a(Index);
        end , clear j
        pxls_traces=pic; clear pic a Index;    
        
        [~ , ind1] = find (finalguideMAT(roi,:)==0);
        high       = prctile(pxls_traces(:,ind1),95);
        low        = prctile(pxls_traces(:,ind1),5);
        highratio  = high./globalyfit(roi,ind1)';
        lowratio   = low./globalyfit(roi,ind1)';
        indicator  = highratio./lowratio;
        [~ , ind2] = find(indicator > prctile(indicator,95));
        no_of_k(roi,1:length(ones(1,length(finalguideMAT))))=ones(1,length(finalguideMAT));
        no_of_k(roi,ind1(ind2))=2;
       
    for frame= 1:size(finalguideMAT,2)
        
        if      finalguideMAT(roi,frame)==0;                                                              %%%%% IN CASE THERE IS NO EVENT %%%%%
               
                %clean_traces(roi,frame)=globalyfit(roi,frame);                                           %% take the baseline (based on the yfit) 
                clean_traces(roi,frame)=baseline_cluster(pxls_traces, frame, no_of_k(roi,:));       %% take the baseline (based on clustering pxls' values)
                
                
        elseif  finalguideMAT(roi,frame)==1;                                                              %%%%% IN CASE THE EVENT IS TRUE %%%%%
               
                %cleantrace(roi,frame)=globaltimemedian (roi,frame);                                      %% take the median
                clean_traces(roi,frame)=globaltimeaverage(roi,frame);                                     %% take the average
                        
                
        elseif  finalguideMAT(roi,frame)==2;                                                              %%%%% IN CASE THE EVENT IS FALSE %%%%%
                if  frame == 1                                                                             %% this is a trick to avoid having problem with indexing in the next line (frame-1)
                        y=find(finalguideMAT(roi,frame:end)~=2);
                        if    isempty(y), xx=frame:size(finalguideMAT,2); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        else,             xx=frame:frame+y(1)-1;          
                        end
                        clean_traces(roi,frame:xx(end))=globalyfit(roi,frame:xx(end));             %% just in case the very first event is false , make it equal the baseline because no other option is available

                elseif       finalguideMAT(roi,frame-1)==2                                                %% this is a trick to let the adjustment take place only at starting points instead of being contineously repeated till 
                    continue                                                                              %% the end of the video
                else
                    y=find(finalguideMAT(roi,frame:end)~=2);
                        if    isempty(y), xx=frame:size(finalguideMAT,2); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        else,             xx=frame:frame+y(1)-1;          %length(find(newguideMAT(frame:frame+y)==2));
                        end
                    interp_line1=drawline([xx(1)-1,xx(end)], [globaltimeaverage(roi,xx(1)-1)                                    , globaltimeaverage(roi,xx(end))],                                    [xx(1)-1:xx(end)]);
                    interp_line2=drawline([xx(1)-1,xx(end)], [globaltimeaverage(roi,xx(1)-1)-globaltimeaverage_out(roi,xx(1)-1) , globaltimeaverage(roi,xx(end))-globaltimeaverage_out(roi,xx(end))], [xx(1)-1:xx(end)]);
                    diffinROIoutROI = globaltimeaverage(roi,xx(1)-1:xx(end))  - globaltimeaverage_out(roi,xx(1)-1:xx(end)); 
                    clean_traces(roi,frame-1:frame+length(interp_line1)-2)=  diffinROIoutROI  +  interp_line1  - interp_line2 ;                %% take interpolated trace 
                end
        end
    end, clear frame pxls_traces
    disp(['ROI ' num2str(roi) ' is cleaned ;)' ])
end, clear roi

toc
%% 
save(Savedir{file});
clearvars -except videos ROIs Savedir
end
