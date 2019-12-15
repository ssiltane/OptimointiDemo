% This routine helps to draw arrows in the plane.
%
% Arguments:
% sp    starting point of the arrow (complex number)
% ep    end point of the arrow (complex number)
% len   number in [0,1], gives length of the desired arrow relative
%       to the distance between start and end points
% alen  number in [0,1], gives length of the arrowheads relative to |sp-ep|
%
% Returns:
% x     vector of x-coordinates 
% y     vector of y-coordinates
%
% Example:
% [ap1,ap2] = arrowpoints(-1,1+i,.98,.03);
% plot(ap1,ap2)
%
% Samuli Siltanen September 2011

function [x,y] = arrowpoints(sp,ep,len,alen)

% Initialize vector of points
vec = sp;

% Calculate true end point
tep = sp + len*(ep-sp);

% Calculate tips of the arrowheads
alength = alen*abs(ep-sp);
unitv1  = exp(i*3*pi/4)*(ep-sp)/abs(ep-sp);
unitv2  = exp(-i*3*pi/4)*(ep-sp)/abs(ep-sp);
ta1     = tep + alength*unitv1;
ta2     = tep + alength*unitv2;

% Collect the results in a complex vector
z = [sp;tep;ta1;tep;ta2];

% Construct the result as real and imaginary parts of z
x = real(z);
y = imag(z);



