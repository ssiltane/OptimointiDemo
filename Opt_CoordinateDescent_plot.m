% Plot coordinate descent method.
%
% The routine Opt_Newton_patch.m needs to be computed before running 
% this file.
%
% Samuli Siltanen Dec 2019

% Number of iterations
Niter = 10;

% Parameters for controlling the plot
lwidth = 3;
msize = 25;
minimsize = 10;
fsize = 20;
colorgray = [.6 .6 .6];
colordark = [.3 .3 .3];
gradplotlen = 1.5;
Ncontours = 100;

% Finite difference stepsize
h = .01;

% Tolerance for small gradient
gradtol = 5*1e-3;
gradnorm = 2*gradtol; % Initialization

% Create grid for plotting
t2MAX = 10;
t1MAX = 1.7778*t2MAX;
t1 = linspace(-t1MAX,t1MAX,512);
t2 = linspace(-t2MAX,t2MAX,1.7778*512);
[X,Y] = meshgrid(t1,t2);

% Evaluate the function to be minimized
minimfun = hillyterrain(X,Y);
figure(2)
clf
mesh(X,Y,minimfun)

% Find minimizer
minimindex = min(find(minimfun==min(minimfun(:))));

% Initialize point matrix with the initial guess
itermat = zeros(2,Niter);
itermat(:,1) = [7;-1];
itermat(:,1) = [8;-1];
% itermat(:,1) = [-5.5;-5];
% itermat(:,1) = [-2;-5];
itermat(:,1) = [-11;-6.5];

% Initialize descent direction matrix
descdirmat = zeros(2,Niter);
descdirmat(1,1) = 1;

% Counter for iterations
iii = 0;

% Initialize image counter
imcounter = 1;

%% Loop over iterations
while (iii<Niter)&(gradnorm>gradtol)
    
    % Increment counter
    iii = iii+1;
    
    % Current iterate
    curx = itermat(1,iii);
    cury = itermat(2,iii);
    
    % Compute current gradient
    curdescdir = descdirmat(:,iii);
    curgradx   = (hillyterrain(curx+h*curdescdir(1),cury)-hillyterrain(curx,cury))/h;
    curgrady   = (hillyterrain(curx,cury+h*curdescdir(2))-hillyterrain(curx,cury))/h;
    
    % Update direction
    descdirmat(:,iii+1) = flipud(curdescdir);
    
    % Determine descent direction (unit vector).
    gradnorm = norm([curgradx;curgrady]); % Checking for convergence
    curdescdir = -[curgradx;curgrady]/norm([curgradx;curgrady]);
    descdirmat(:,iii) = curdescdir;
    
    % Find the minimal point along the negative gradient direction
    [x,y] = findminimalpoint(curx,cury,curdescdir,@hillyterrain,h,t1MAX,t2MAX);
    itermat(:,iii+1) = [x;y];
    
    % Plot the function to be minimized
    figure(1)
    clf
    contour(X,Y,minimfun,Ncontours)
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
    %plot(curx+[0 gradplotlen*curdescdir(1)],cury+[0 gradplotlen*curdescdir(2)],'r','linewidth',lwidth)
    [ap1,ap2] = arrowpoints(curx+1i*cury,curx+gradplotlen*curdescdir(1)+1i*(cury+gradplotlen*curdescdir(2)),gradplotlen,.4);
    plot(ap1,ap2,'r','linewidth',lwidth)
    
    % Plot minimizer
    plot(X(minimindex),Y(minimindex),'b.','markersize',minimsize)
    
    % Axis settings
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    axis equal
    axis([-t1MAX t1MAX -t2MAX t2MAX])
    axis off
    %title('Coordinate descent method','fontsize',fsize)
    
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
    filename = ['frames_Coord/OptFrameCoord_',num2str(imcounter),'.png'];
    imwrite(uint8(im3),filename,'png');
    imcounter = imcounter+1;
end


%% Draw final image

% Plot the function to be minimized
figure(1)
clf
contour(X,Y,minimfun,Ncontours)
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
%title('Coordinate descent method','fontsize',fsize)


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
filename = ['frames_Coord/OptFrameCoord_nopatch_',num2str(imcounter),'.png'];
imwrite(uint8(im3),filename,'png');

