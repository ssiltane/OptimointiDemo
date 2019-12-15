% Plot Newton's method, showing only local patches of the objective
% function
%
% The routine Opt_Newton_patch.m needs to be computed before running 
% this file.
%
% Samuli Siltanen Dec 2019

% Number of iterations
Niter = 60;

% Relaxation parameter
relparam = .8;

% Pseudoinverse tolerance parameter
pinvparam = 0.003;

% Parameters for controlling the plot
lwidth = 3;
thinline = 1;
verythinline = 1;
msize = 25;
minimsize = 10;
fsize = 20;
colorgray = [.6 .6 .6];
colordark = [.3 .3 .3];
gradplotlen = 1.5;
Ncontours = 100;
Rpatch = 2.5;
Nfii = 64;
fii = linspace(0,2*pi,Nfii);

% Finite difference stepsize
h = .1;

% Tolerance for small gradient
gradtol = 5*1e-2;
gradnorm = 2*gradtol; % Initialization

% Load the function from file
load data/minimfun minimfun funMIN funMAX contourvec minimindex X Y t1MAX t2MAX Ncontours

% Initialize point matrix with the initial guess
itermat = zeros(2,Niter);
itermat(:,1) = [8;-1];
itermat(:,1) = [9;-1];

% Initialize jump direction matrix
descdirmat = zeros(2,Niter);

% Initialize image counter
imcounter = 1;

%% Loop over iterations

% Counter for iterations
iii = 0;

while (iii<Niter)&(gradnorm>gradtol)
    % Increment counter
    iii = iii+1;
    
    % Current iterate
    curx = itermat(1,iii);
    cury = itermat(2,iii);
    
    % Compute current gradient
    curgradx = (hillyterrain(curx+h,cury)-hillyterrain(curx,cury))/h;
    curgrady = (hillyterrain(curx,cury+h)-hillyterrain(curx,cury))/h;
    
    % Compute current Hessian
    curHessian11 = (hillyterrain(curx+h,cury)-2*hillyterrain(curx,cury)+hillyterrain(curx-h,cury))/h^2;
    curHessian22 = (hillyterrain(curx,cury+h)-2*hillyterrain(curx,cury)+hillyterrain(curx,cury-h))/h^2;
    curHessian12 = (hillyterrain(curx+h,cury+h)-hillyterrain(curx+h,cury)-hillyterrain(curx,cury+h)+hillyterrain(curx,cury))/h^2;
    curHessian21 = curHessian12;
    curHessian = [[curHessian11,curHessian12];[curHessian21,curHessian22]];
    
    % Determine descent direction (unit vector).
    gradnorm = norm([curgradx;curgrady]); % Checking for convergence
    curdescdir = -[curgradx;curgrady]/norm([curgradx;curgrady]);
    descdirmat(:,iii) = curdescdir;
    
    % Find the next iterate
    nextpoint = [curx;cury]-relparam*pinv(curHessian,pinvparam)*[curgradx;curgrady];
    itermat(:,iii+1) = nextpoint;
    
    % Create plot window
    figure(1)
    clf
    
    % Plot the function to be minimized; just a patch round current point visible
    patchind = abs(X+1i*Y-(curx+1i*cury))<Rpatch;
    minimfun_patch = minimfun;
    minimfun_patch(~patchind) = NaN;
    contour(X,Y,minimfun_patch,contourvec);
    set(gca,'CLim',[funMIN,funMAX])
    hold on
    plot(curx+Rpatch*cos(fii),cury+Rpatch*sin(fii),'k','linewidth',verythinline)
    
    % Plot history of iterates and descent directions
    for jjj = 1:(iii-1)
        p2 = plot([itermat(1,jjj),itermat(1,jjj+1)],[itermat(2,jjj),itermat(2,jjj+1)],'k','linewidth',lwidth);
        set(p2,'color',colorgray)
        p1 = plot(itermat(1,jjj),itermat(2,jjj),'k.','markersize',msize);
        set(p1,'color',colordark)
    end
    
    % Plot current iterate and current neighborhood analysis
    [V,D] = eig(curHessian);
    eigvec = diag(D);
    eigvec2 = eigvec/max(abs(eigvec));
    % Vertices determined by the Hessian
    v1x = curx-2*eigvec2(1)*V(1,1);
    v1y = cury-2*eigvec2(1)*V(2,1);
    v2x = curx+2*eigvec2(1)*V(1,1);
    v2y = cury+2*eigvec2(1)*V(2,1);
    v3x = curx-2*eigvec2(2)*V(1,2);
    v3y = cury-2*eigvec2(2)*V(2,2);
    v4x = curx+2*eigvec2(2)*V(1,2);
    v4y = cury+2*eigvec2(2)*V(2,2);
    plot([v1x,v2x],[v1y,v2y],'b','linewidth',thinline)
    plot([v1x,v3x],[v1y,v3y],'b','linewidth',thinline)
    plot([v1x,v4x],[v1y,v4y],'b','linewidth',thinline)
    plot([v3x,v4x],[v3y,v4y],'b','linewidth',thinline)
    plot([v2x,v4x],[v2y,v4y],'b','linewidth',thinline)
    plot([v2x,v3x],[v2y,v3y],'b','linewidth',thinline)
    %     plot(curx+[-1;1]*2*eigvec2(1)*V(1,1),cury+[-1;1]*2*eigvec2(1)*V(2,1),'b','linewidth',lwidth)
    %     plot(curx+[-1;1]*2*eigvec2(2)*V(1,2),cury+[-1;1]*2*eigvec2(2)*V(2,2),'b','linewidth',lwidth)
    [ap1,ap2] = arrowpoints(curx+1i*cury,curx+gradplotlen*curdescdir(1)+1i*(cury+gradplotlen*curdescdir(2)),gradplotlen,.4);
    plot(ap1,ap2,'r','linewidth',lwidth)
    plot(curx,cury,'r.','markersize',msize)
    
    % Plot minimizer
    %plot(X(minimindex),Y(minimindex),'b.','markersize',minimsize)
    
    % Axis settings
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    axis equal
    axis([-t1MAX t1MAX -t2MAX t2MAX])
    axis off
    %title('Newton method','fontsize',fsize)
    
    % Grab the plot into a color image matrix
    im1 = print('-r500','-RGBImage');
    [row1,col1] = size(im1(:,:,1));
    
    % Crop the image
    startrow = round(.2*row1);
    endrow   = round(.76*row1);
    startcol = round(.2*col1);
    endcol   = round(.8*col1);
    im2 = im1(startrow:endrow,startcol:endcol,:);
    
    % Adjust image size to 1080x1920
    im2 = imresize(im2, [1080 NaN]);
    [~,col2,~] = size(im2);
    im3 = uint8(5*ones(1080,1920,3));
    im3(:,round((1920-col2)/2)+[1:col2],:) = im2;
    
    % Save to file
    filename = ['frames_Newton_patchB/OptFrameNewtonB_',num2str(imcounter),'.png'];
    imwrite(uint8(im3),filename,'png');
    imcounter = imcounter+1;
    disp([iii Niter])
end


%% Draw final image

% Plot the function to be minimized
figure(1)
clf
contour(X,Y,minimfun,contourvec)
hold on

% Plot history of iterates and descent directions
for jjj = 1:(iii-1)
    p2 = plot([itermat(1,jjj),itermat(1,jjj+1)],[itermat(2,jjj),itermat(2,jjj+1)],'k','linewidth',lwidth);
    set(p2,'color',colorgray)
    p1 = plot(itermat(1,jjj),itermat(2,jjj),'k.','markersize',msize);
    set(p1,'color',colordark)
end

% Plot current iterate and current descent direction
plot(curx,cury,'r.','markersize',msize)

% Plot minimizer
plot(X(minimindex),Y(minimindex),'b.','markersize',minimsize)

% Axis settings
set(gca,'xtick',[])
set(gca,'ytick',[])
axis equal
axis([-t1MAX t1MAX -t2MAX t2MAX])
axis off

% Grab the plot into a color image matrix
im1 = print('-r500','-RGBImage');
[row1,col1] = size(im1(:,:,1));

% Crop the image
startrow = round(.2*row1);
endrow   = round(.76*row1);
startcol = round(.2*col1);
endcol   = round(.8*col1);
im2 = im1(startrow:endrow,startcol:endcol,:);

% Adjust image size to 1080x1920
im2 = imresize(im2, [1080 NaN]);
[~,col2,~] = size(im2);
im3 = uint8(5*ones(1080,1920,3));
im3(:,round((1920-col2)/2)+[1:col2],:) = im2;

% Save to file
filename = ['frames_Newton_patchB/OptFrameNewtonB_',num2str(imcounter),'.png'];
imwrite(uint8(im3),filename,'png');


%% Plot very first image

% Current iterate
curx = itermat(1,1);
cury = itermat(2,1);

% Create plot window
figure(1)
clf

% Plot the function to be minimized
contour(X,Y,minimfun,contourvec)
hold on

% Plot initial point
curx = itermat(1,1);
cury = itermat(2,1);
plot(curx,cury,'r.','markersize',msize)

% Plot minimizer
plot(X(minimindex),Y(minimindex),'b.','markersize',minimsize)

% Axis settings
set(gca,'xtick',[])
set(gca,'ytick',[])
axis equal
axis([-t1MAX t1MAX -t2MAX t2MAX])
axis off
%title('Newton method','fontsize',fsize)

% Grab the plot into a color image matrix
im1 = print('-r500','-RGBImage');
[row1,col1] = size(im1(:,:,1));

% Crop the image
startrow = round(.2*row1);
endrow   = round(.76*row1);
startcol = round(.2*col1);
endcol   = round(.8*col1);
im2 = im1(startrow:endrow,startcol:endcol,:);

% Adjust image size to 1080x1920
im2 = imresize(im2, [1080 NaN]);
[~,col2,~] = size(im2);
im3 = uint8(5*ones(1080,1920,3));
im3(:,round((1920-col2)/2)+[1:col2],:) = im2;

% Save to file
filename = ['frames_Newton_patchB/OptFrameNewtonB_0.png'];
imwrite(uint8(im3),filename,'png');
imcounter = imcounter+1;
