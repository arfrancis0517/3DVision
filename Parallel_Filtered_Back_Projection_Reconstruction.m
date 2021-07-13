%% Reconstruction Preparation
% Guo Tang
% 01.June 2021
% References: 
% Radermacher M. (2007) Weighted Back-projection Methods. 
% In: Frank J. (eds) Electron Tomography. Springer, New York, NY.
clc,clear
tic

path ='/Volumes/bigbigfranc/data'; % Set input_folder
filename = [path filesep '*.mat']; % Define the file path
files=dir(filename); 
L = length(files)-1; % Get the number of the .mat projection files
all_theta = zeros(1,L); % Prepare a theta collector
load([path filesep files(2).name]); % Load a example file to check the data in the file
an = ar/n; % Calculate the angle increment
endang = ar + startang; % The end angle
[l,z] = size(shifted_noisy_proj2d); % Get the size of shifted noisy 2D projection.
sum_proj1d = zeros(l,L); % Prepare a zero-valued 1D projection collector.
N = 2*floor(l/(2*sqrt(2))); % the number of rows and columns in the reconstructed image.
img2d_new = zeros(N);        % Prepare a zero-valued 2D image.
img3d_new = zeros(N,N,z); % Prepare a zero-valued 3D image.
% Set the x & y axis for the reconstructed image.
dist = (1:N)-ceil(N/2); % the distance to the center of the reconstructed image
x = repmat(dist, N, 1);    % the x coordinates
y = rot90(x); % the y coordinates
ctrIdx = ceil(l/2);     % index of the center of the projections
imgDiag = 2*ceil(N/sqrt(2))+1;  % largest distance through image.

%% Design a filter in Fourier space
% First design a Ramp Filter
order = max(64,2^nextpow2(2*l)); % Speed up the FFT
filter = 2*( 0:(order/2) )./order;
w = 2*pi*(0:size(filter,2)-1)./order;   % frequency axis up to Nyquist
d = 1; % d is a scalar in the range (0,1] that modifies the filter by 
       % rescaling its frequency axis.  The default is 1. 

% The Ramp filter can be regulated by some window functions:

% The Shepp-Logan filter multiplies the ramp filter by a sinc function.
% filter(2:end) = filter(2:end) .* (sin(w(2:end)/(2*d))./(w(2:end)/(2*d)));

% The hamming filter with cosine function
% filter(2:end) = filter(2:end) .*(.54 + .46 * cos(w(2:end)/d));

% The hann filter with cosine function
% filter(2:end) = filter(2:end) .*(1+cos(w(2:end)./d)) / 2;

filter(w>pi*d) = 0; % Crop the frequency response, all frequencies above d are set to 0
filter = [filter' ; filter(end-1:-1:2)'];    % Symmetry of the filter

figure
plot(filter)
title('Filter');
xlabel('k_y');
ylabel('Factor');
xlim('tight')
ylim('tight')
set(gca,'FontSize',12);

%% Collecte the angles
j = 1;
    for ang = startang : an : endang
        filename = [path filesep 'Proj_' sprintf('%.2f',ang) '°.mat'];        
        load(filename);
        all_theta(1,j) = ang;
        j = j + 1;       
    end  
    
theta = pi * all_theta/180;
costheta = cos(theta);
sintheta = sin(theta);

%% Reconstruction

for i = 1:z    
    k=1; 
    % Collecte the 1D projection of same sub-2D-image of each angle
    for ang = startang : an : endang
        filename = [path filesep 'Proj_' sprintf('%.2f',ang) '°.mat'];        
        load(filename);
        subproj = shifted_noisy_proj2d(:,i);
        sum_proj1d(:,k) = subproj;  
        k = k+1;
    end
    proj_body = sum_proj1d; % sinogram
    proj_body(length(filter),1)=0; % Zero pad projections
    % Use 2D weighted back projection to get all new-sub-2D-images
    f = fft(proj_body);    % f holds fft of projections
    for a = 1:size(f,2)
        f(:,a) = f(:,a).* filter; % frequency domain filtering
    end
    f = real(ifft(f));   % f is the filtered projections
    f(l+1:end,:) = [];   % Truncate the filtered projections

    if size(f,1) < imgDiag 
        rz = imgDiag - size(f,1);  % how many rows of zeros
        f = [zeros(ceil(rz/2),size(f,2)); f; zeros(floor(rz/2),size(f,2))];
        ctrIdx = ctrIdx+ceil(rz/2);
    end

    % linear interpolation
    for c=1:length(theta)  
        proj = f(:,c);
        t = x.*costheta(c) + y.*sintheta(c); 
        tt = floor(t);  
        img2d_new = img2d_new + (t-tt).*proj(tt+1+ctrIdx) + (tt+1-t).*proj(tt+ctrIdx);
    end
    
    img2d_new = img2d_new*pi/(2*length(theta));
    
    % Build a new 3D image
    img3d_new(:,:,i) = img2d_new;        
end
% 
toc
%% Visualization and Save

% output_folder = '/Volumes/bigbigfranc/data'; % Set the output folder 
% filename = [output_folder filesep '3DImage_new.mat'];
% save(filename,'img3d_new','all_theta','ar');

figure
intensity = [-3024,-16.45,641.38,3071];
alpha = [0, 0, 0.72, 0.72];
color = ([0 0 0; 186 65 77; 231 208 141; 255 255 255]) ./ 255;
queryPoints = linspace(min(intensity),max(intensity),256);
alphamap = interp1(intensity,alpha,queryPoints)';
colormap = interp1(intensity,color,queryPoints);
volshow(img3d_new,'Colormap',colormap,'Alphamap',alphamap)

output_folder = '/Volumes/bigbigfranc/data';
filename = [output_folder filesep '3DImage_new.fig'];
savefig(filename)
filename = [output_folder filesep '3DImage_new.png'];
saveas(gcf,filename)