function [cdata] = despike(datain, slope)
% despike.m -> de-spikes a time-series by calculating point-to-point slopes
% determining whether a maximum allowable slope is exceeded
%
%   call:  [cdata] = despike(datain, slope);
%
%  input:  datain -> input data (e.g.: time-series)
%          slope  -> dlope threshold
%
% output:  cdata  -> clean data
%
% example: [cdata] = despike(datain, slope);
%
% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     10-Jan-2008
% modified: 13-Aug-2010
%

data_holder(1) = datain(1);

ii = 1;
number_of_points = length(datain);

for i = 2:number_of_points
    new_slope = datain(i) - data_holder(ii);
    % If the slope is okay, let the data through:
    if ( abs(new_slope) <= abs(slope) )
        ii = ii+1;
        data_holder(ii) = datain(i);
    else
    % Slope not okay, so look for the next data point:
        n=0;
        while ( (abs(new_slope) > abs(slope)) & ( (i+n) < number_of_points) )
        n = n+1; % Index offset of test point
        numerator   = datain(i+n) - data_holder(ii);
        denominator = (i+n) - ii;
        new_slope = numerator/denominator;
        end
        % If we have a "good" slope, calculate new point using linear interpolation:
        % point = {[(ngp - lgp)/(deltax)]*(actual distance)} + lgp
        % ngp = next good point
        % lgp = last good point
        % actual distance = 1, the distance between the last lgp and the point we want to interpolate
        % Otherwise, let the value through (i.e. we've run out of good data):
        if ( (i+n) < number_of_points )
            the_point = new_slope + data_holder(ii);
            ii = ii+1;
            data_holder(ii) = the_point;
        else
            ii = ii+1;
            data_holder(ii) = datain(i);
        end
    end
end

cdata = data_holder(:);