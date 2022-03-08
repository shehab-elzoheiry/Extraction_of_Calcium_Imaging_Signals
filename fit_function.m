function [old_yfit]=shehab_oldfit(timeaverage,interv_length)  
    int=interv_length;                                                                                                       % interval length (in data points)
    for jj=1:size(timeaverage,2)/int
        %[~,ind]=find(timeaverage( (jj-1)*int+1 : (jj)*int )== (min(timeaverage( (jj-1)*int+1 : (jj)*int )))   );            % choose the min within the interval
        if   jj== size(timeaverage,2)/int                                                                                    % this is the end point
             xindices(jj)= jj*int;
        else [~,ind]=find(timeaverage( (jj-1)*int+1 : (jj)*int +1)== (median(timeaverage( (jj-1)*int+1 : (jj)*int+1 )))   ); % choose the median within the interval (the extra 1 compared to the next line is to have a string of odd numbers so that to get a real medaian)
            %[~,ind]=find(timeaverage( (jj-1)*int+1 : (jj)*int )== (min(timeaverage( (jj-1)*int+1 : (jj)*int )))   );        % choose the min within the interval
             xindices(jj)= ind(1) + (jj-1)*int ;                                                                             % adjust the index of the begining according to where it should be in the full trace
        end
    end
    clear jj int ind
    %[ycorr,yfit] = bf((timeaverage)',xindices,'confirm');                                                                   % ploting the baseline 
    x=size(timeaverage,2);
    [ycorr,old_yfit] = bf((timeaverage)', [1, [0.2*x 0.4*x 0.6*x 0.8*x], size(timeaverage,2)] ,0.2*x,'pchip');                       %linear fitting between chosen points
    