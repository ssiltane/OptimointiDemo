% Find the point where the given function has the smallest value along a
% designated half-line. The search along the half-line happens in finite
% difference steps and is bounded by the box [-MAX,MAX]x[-MAX,MAX].
%
% Arguments:
% x0       x-coordinate of the starting point of half-line
% y0       y-coordinate of the starting point of half-line
% dirvec   direction of the half-line from the starting point (unit vector)
% f        handle to the function (of two real variables) to be minimized
% h        steplength in the primitive line search
% MAX      bounding box limit for ending the search
%
% Returns: 
% x,y      the coordinates of the location of minimal value
% xvec     vector of x-coordinates along the half-line
% yvec     vector of y-coordinates along the half-line
% funvals  function values along the half-line
% minind   location index of the minimal function value
%
% Samuli Siltanen Nov 2019

function     [x,y,xvec,yvec,funvals,minind] = findminimalpoint(x0,y0,dirvec,f,h,xMAX,yMAX)

% Construct evaluation points
xvec = x0;
yvec = y0;
counter = 1;
while max(abs(xvec))<xMAX & max(abs(yvec))<yMAX
    xvec = [xvec,x0+counter*h*dirvec(1)];
    yvec = [yvec,y0+counter*h*dirvec(2)];
    counter = counter+1;
end

% Evaluate the function along the half-line
funvals = f(xvec,yvec);

% Find minimal point
minind = min(find(funvals==min(funvals)));
x = xvec(minind);
y = yvec(minind);


