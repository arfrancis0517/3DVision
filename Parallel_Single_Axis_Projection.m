%% Load the image date and Visualization
% Guo Tang
% 01.June 2021
% References: 
% Radermacher M. (2007) Weighted Back-projection Methods. 
% In: Frank J. (eds) Electron Tomography. Springer, New York, NY.
clc,clear


img3d= importdata('volume.mat'); 
[x,y,z]=size(img3d); % get the size of 3d image

% Visualize the 3D image using volshow
figure
intensity = [-3024,-16.45,641.38,3071];
alpha = [0, 0, 0.72, 0.72];
color = ([0 0 0; 186 65 77; 231 208 141; 255 255 255]) ./ 255;
queryPoints = linspace(min(intensity),max(intensity),256);
alphamap = interp1(intensity,alpha,queryPoints)';
colormap = interp1(intensity,color,queryPoints);
volshow(img3d,'Colormap',colormap,'Alphamap',alphamap)

% Save the figure
output_folder = '/Volumes/bigbigfranc/data';
filename = [output_folder filesep 'Original_img.fig'];
savefig(filename)
filename = [output_folder filesep 'Original_img.png'];
saveas(gcf,filename)

%% Projection Preparation From -60° to 60° with angle increment 0.12°
n = 1000; % Projection number, I have generated n+1 projection
ar = 120; % Angle range
an = ar / n;  % Angle increment of each projection
startang = - 60; % Start angle 
endang = ar + startang; % End angle
% Prepare to pad the image with zeros so we don't lose anything when we rotate.
diag = sqrt(x^2 + y^2);
lenpad = ceil(diag-x)+2;
l = x + lenpad;
widpad = ceil(diag-y)+2;
w = y + widpad;
imgpad = zeros(l,w); 
% Prepare a zero-valued 2D Projection, 1D Projection and a 2D grid
proj2d = zeros(l,z); 
proj1d = zeros(l,1);
m = linspace(-1,1,l);
[X1,Y1] = meshgrid(m,m);

%% Projection, Visualization and Save

for ang = startang : an : endang
    % for each angle ang, we loop all 2D image slice of
    % the original 3D volume along z axis, and rotate the
    % 2D image with the angle ang,  sum up the values along
    % the new y axis to get a vector column, combine all
    % 1D vector columns to build 2D Projection with angle ang.
    for i = 1:z
        subimg2d = img3d(:,:,i); % 2D image slice
        % Pad each sub-2D-image with zeros
        imgpad(ceil(lenpad/2):(ceil(lenpad/2)+x-1), ... 
               ceil(widpad/2):(ceil(widpad/2)+y-1)) = subimg2d;
        %loop over the number of angles, rotate theta degree, and then add up.
        theta = ang * pi / 180;
        X = cos(theta)*X1 + -sin(theta)*Y1;
        Y = sin(theta)*X1 + cos(theta)*Y1; 
        % perform linear interpolation on the rotating to get 1D projection of sub-2D-image.
        tmpimg = interp2(X1,Y1,imgpad,X,Y);
        tmpimg(isnan(tmpimg)) = 0; % set NaN to zero
        proj1d = (sum(tmpimg))';
        % Add each 1D projection from every sub-2D-image with same angle to build 2D projection.
        proj2d(:,i) = proj1d;
    end
    
    % Add Gaussian noise to 2D projection with amplitude
    amp = 200; % Amplitude
    noisy_proj2d = proj2d + amp * randn(size(proj2d)); % Noisy 2D projection
    
    dim = round(rand)+1; % Random orientation of image shifting 
    e = 8; % Shifting coefficient
    K = round(round(size(proj2d,2)*(rand-0.5))/e); % Shifting distance with random orientation
    shifted_noisy_proj2d = circshift(noisy_proj2d,K,dim); % Shifted noisy 2D projection
    
    output_folder = '/Volumes/bigbigfranc/data'; % Set the output folder 
    filename = [output_folder filesep 'Proj_' sprintf('%.2f',ang) '°.mat']; % Set the filename
	save(filename,'n','ar', 'startang', 'proj2d','noisy_proj2d','amp','shifted_noisy_proj2d','dim','K','ang');
    
    
    figure
    image(proj2d,'CDataMapping','scaled');
    title(['Projection of Object With ', sprintf('%.2f',ang),' Degree']);
    xlabel('X');
    ylabel('Y');
    axis equal
    xlim('tight')
    ylim('tight')
    set(gca,'FontSize',12);
    colormap jet
    colorbar
    

    filename = [output_folder filesep 'Proj_' sprintf('%.2f',ang) '°.fig'];
    savefig(filename)
    filename = [output_folder filesep 'Proj_' sprintf('%.2f',ang) '°.png'];
    saveas(gcf,filename)

    
    figure
    image(noisy_proj2d,'CDataMapping','scaled');
    title(['Noisy Projection of Object With ', sprintf('%.2f',ang),' Degree']);
    xlabel('X');
    ylabel('Y');
    axis equal
    xlim('tight')
    ylim('tight')
    set(gca,'FontSize',12);
    colormap jet
    colorbar
    
    filename = [output_folder filesep 'Noisy_proj_' sprintf('%.2f',ang) '°.fig'];
    savefig(filename)
    filename = [output_folder filesep 'Noisy_proj_' sprintf('%.2f',ang) '°.png'];
    saveas(gcf,filename)

    
    figure
    image(shifted_noisy_proj2d,'CDataMapping','scaled');
    title(['Shifted Noisy Projection of Object With ', sprintf('%.2f',ang),' Degree']);
    xlabel('X');
    ylabel('Y');
    axis equal
    xlim('tight')
    ylim('tight')
    set(gca,'FontSize',12);
    colormap jet
    colorbar
    
    
    filename = [output_folder filesep 'Shifted_noisy_proj_' sprintf('%.2f',ang) '°.fig'];
    savefig(filename)
    filename = [output_folder filesep 'Shifted_noisy_proj_' sprintf('%.2f',ang) '°.png'];
    saveas(gcf,filename)
    close all

end