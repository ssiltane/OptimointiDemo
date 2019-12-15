% Function of two real variables to be minimized
%
% Samuli Siltanen Nov 2019

function res = hillyterrain(x,y)

% Check that x and y have the same size
[rowx,colx] = size(x);
[rowy,coly] = size(y);
if ~(rowx==rowy) | ~(colx==coly)
    error('hillyterrain: x and y must have same size')
end

% Construct pieces of the function
valley1 = 1+abs(.07*(x+1i*y)).^2;
theta1 = pi/9;
xx = cos(theta1)*x - sin(theta1)*y;
yy = sin(theta1)*x + cos(theta1)*y;
fun1 = 4-cos(.4*(xx)).*cos(.3*(yy)+.05*yy.^2);

% Combine the results 
res = valley1.*fun1;
