%% Code to construct 2D maps of Shear Wave Velocity (SWV)
%
% Applied to ARFI data from the IEEE UFFC-JS 2024 Challenge Algorithms for 
% Mapping Elastic Modulus (a-MEM) 
%
% Please refer to corresponding proceeding when using the code: 
% Cihan A., Segers P. and Caenen A., "A Robust Pipeline Including Directional 
% Filtering and Wave Propagation Zones for Mapping Shear Wave Velocity",
% IEEE International Ultrasonics Symposium (IUS), 2024.
%
% Organization:     Institute for Biomedical Engineering and Technology,
%                   Ghent University, Ghent, Belgium
%
% Authors:          Ariana Cihan, Annette Caenen
%
% Contact:          ariana.cihan@ugent.be
%
% Input data:       ARF_phantom_**.mat-file containing the following variables:
%                   - IQ_Data: 1x2 cell, containing 2 3D datasets of IQ data (z,x,time) 
%                   - frame_rate: frame rate [Hz]
%                   - X_IQ: vector depicting horizontal (x) axis [mm]
%                   - Z_IQ: vector depicting vertical (z) axis [mm]
%                   - X_push: X-coordinate of push beam [mm]
%                   - x_circle, z-circle: coordinates of inclusion in phantom
%
% Output data:      ***.mat- and ***.png-files containing the variable:
%                   - elastogram: 2D matrix depicting shear wave velocity (SWV) [m/s]
%
% Description:      Elastogram is calculated based on cross-correlation
%                   in 2D, after applying directional and bandpass
%                   filtering, for zones outside the push beam region and
%                   with significant wave propagation
%
% Last updated: 17th of October 2024
% Written in MATLAB 2023a
clear; clc; close all;

%% Defining directories
% Start directory
startdir = 'C:\Users\name\...\';
% Directory where data is stored
datadir = [startdir 'Data\'];
% Directory were results are stored
BFdir = [startdir 'a-MEM_Results\'];
% Subdirectory where results are stored
savedir  = [BFdir 'OverviewResults_date'];
if ~exist(savedir,'dir'); mkdir(savedir); end

%% Choosing phantom data
% Phantom data provided by the IEEE UFFC-JS 2024 a-MEM challenge
phantom = 3;

% Selecting filenames and inclusion boundaries (determined from visual inspection)
% x_circle, z_circle: coordinates of inclusion boundary [mm]
theta = linspace(0,2*pi,100); % angles for plotting inclusions [rad]
if phantom == 1
    % No inclusion is present in phantom 1
    file = 'ARF_phantom_I.mat';
    theta = 0; x_circle = 0; z_circle = 0;
elseif phantom == 2
    % Stiffer inclusion at top
    file = 'ARF_phantom_II.mat';
    x_circle = 8*cos(theta); z_circle = 16+8*sin(theta);
elseif phantom == 3
    % Softer inclusion at bottom
    file = 'ARF_phantom_III.mat';
    x_circle = 10*cos(theta); z_circle = 35+8*sin(theta);
else 
    error('Wrong phantom selected')
end 

%% Parameters
% Set size of figures
scrsize = get(0, 'ScreenSize');
figset1 = scrsize(3)/4;
figset2 = scrsize(4)/2;

% Switches (1=on, 0=off)
SHOW_BMODE              = 0; % to show B-mode
SAVE_MOVIE_BMODE        = 1; % to save B-mode movie
SHOW_ITERIM_STEPS       = 0; % to show debug figures
SHOW_TDI                = 0; % to show wave propagation movie
SAVE_MOVIE_TDI          = 1; % to save wave propagation movie
PUSH_SELECTION          = 1:2; % selection pushes of interest

% Settings particle motion
TDIcolorscale           = [-1.5 1.5]*10^-3; % colorscale for axial particle velocities [m/s]
SWVcolorscale           = [0 3]; % colorscale for shear wave velocity [m/s] 
Bmode_dBlevel           = [-70 -10]; % colorscale for B-mode [dB]
f_center                = 5e6; % center frequency of push [Hz]

% Temporal and spatial smoothing
FILT_TEMP_SPAT          = 1; % spatial and temporal smoothing filter (applied after one-lag autocorrelation technique)
N_spatialz              = 7; % size of smoothing filter (# pixels in z-direction)
N_spatialx              = 7; % size of smoothing filter (# pixels in x-direction)
N_temporal              = 1; % size of smoothing filter (# time frames)

% Shear wave speed estimation method
t_ind                   = 8:100; % indices of temporal range to track shear waves
numref                  = t_ind(1); % start of timeframe to consider
TEMP_UPSAMP             = 1; % temporal upsampling 
BANDPASS_FILT           = 1; % bandpass filtering
DIRECTIONAL_FILT        = 1; % directional filtering
margin_mm               = 3; % margin considered for push mask [mm]

% SWV map (elastogram)
w                       = 12; % # pixels in window when calculating SWV
p                       = 8; % # pixels in window when calculating SWV
Vmin                    = 0.5; % minimal SWV [m/s]
Vmax                    = 6.0; % maximal SWV [m/s]
SMOOTHING_MAP           = 1; % smoothing of elastogram

% Set spatial smoothing filter
if FILT_TEMP_SPAT 
    mk3Dwin = zeros(N_spatialz,N_spatialx,N_temporal);
    for kk=1:N_temporal
        mk3Dwin(:,:,kk) = ( gausswin(N_spatialz) * gausswin(N_spatialx)' );
    end
    for kk=1:N_spatialz
        for mm=1:N_spatialx
            mk3Dwin(kk,mm,:) = squeeze(mk3Dwin(kk,mm,:)) .* gausswin(N_temporal);
        end
    end
end

%% Start Analysis
load([datadir file])

% Create folder for processed results per measurement
idx_ = find(file == '.')-1;
filename = [file(1:idx_(end)) '_Results'];
PROC_DAT_DIR = fullfile(savedir, filename);
if ~exist(PROC_DAT_DIR,'file'); mkdir(PROC_DAT_DIR); end

% Create results directory and convert number of acquisitions
dirname=[BFdir  filename '\'];

% Parameters for TDI
PRF = frame_rate; % pulse repetition frequency
dt = 1/PRF;

%% Defining ROI's
% Used to calculate metrics
% ROI1 = inclusion mask
% ROI2 = mask for background around inclusion
Nz = length(Z_IQ); Nx = length(X_IQ);
[Xaxis,Zaxis] = meshgrid(X_IQ, Z_IQ);

if phantom == 1
    ROI1 = zeros(Nz,Nx,'logical');
elseif phantom == 2
    ROI1 = ((Xaxis).^2+(Zaxis-16).^2 <= 8^2);
elseif phantom == 3
    ROI1 = ((Xaxis/10).^2+((Zaxis-35)/8).^2 <= 1);
end 
ROI2 = ~ROI1;
ROI2([1:w end-w:end],:) = 0; % discard edges
ROI2(:,[1:w end-w:end]) = 0;

if SHOW_ITERIM_STEPS
    figure('Position', [50 100 figset1*2 figset2],'PaperPositionMode','auto')
    tiledlayout(1,2,TileSpacing="tight",Padding="tight");
    nexttile 
    imagesc(X_IQ,Z_IQ,ROI1)
    axis image
    colormap gray;
    title('ROI1 - Inclusion')
    xlabel('x [mm]',FontSize=13); ylabel('z [mm]',FontSize=13);
    nexttile
    imagesc(X_IQ,Z_IQ,ROI2)
    axis image
    colormap gray;
    title('ROI2 - Background')
    xlabel('x [mm]',FontSize=13); ylabel('z [mm]',FontSize=13);
    colorbar
end

%% Concatenate Data - Bmode
numAcqs = length(PUSH_SELECTION); % number of acquisitions
numFrames = size(IQ_Data{1},3); % number of frames per acquisition
dX = mean(diff(X_IQ)); % pixel size in x-direction for IQ data [mm]
dZ = mean(diff(Z_IQ)); % pixel size in z-direction for IQ data [mm]
taxis = 0:dt:(numFrames-1)*dt; % time axis [s]

Frames = zeros(length(Z_IQ),length(X_IQ),numFrames,numAcqs); % z,x,t,push_selection
for i = PUSH_SELECTION
    Frames(:,:,:,i) = IQ_Data{i};
end 
FramesCart = Frames;  %scanconvert Frames
max_Frames  = max(abs(FramesCart(:)));

%% Show Bmode
if SHOW_BMODE
    name_Bmodemovie = [PROC_DAT_DIR '\Bmode.gif'];
    gcf3 = figure('Color','white','Position', [250 250 400 400],'PaperPositionMode','auto');
    for j = 1:numAcqs
        for kk = 1:numFrames
            imagesc(X_IQ,Z_IQ,20*log10(abs(squeeze(FramesCart(:,:,kk,j)))./max_Frames))
            hold on 
            %plot(x_circle,z_circle,Color=0.5*[1 1 1],LineWidth=0.5)
            axis image
            title(sprintf('Push %.0f - Time %.2f ms',j,taxis(kk)*1e3))
            xlabel('x [mm]');ylabel('z [mm]');
            colorbar
            clim(Bmode_dBlevel)
            colormap gray

            if SAVE_MOVIE_BMODE
                frame = getframe(gcf3);
                im = frame2im(frame);
                [A,map] = rgb2ind(im,256);

                if kk == 1 && j == 1
                    imwrite(A,map,name_Bmodemovie,'gif','LoopCount',Inf,'DelayTime',0.15);
                else
                    imwrite(A,map,name_Bmodemovie,'gif','WriteMode','append','DelayTime',0.1);
                end
            end
        end
    end
end


%% Compute displacements/TDI - Kasai one lag autocorrelation
TDI_dat_Angles = zeros(size(Frames));
for j = 1:size(Frames,4)
    TDI_dat_Angles(:,:,2:size(Frames,3),j) = (Frames(:,:,2:end,j)) .* conj((Frames(:,:,1:end-1,j)));
end

% Temporal and spatial smoothing in polar domain
if FILT_TEMP_SPAT
    for i = 1:size(TDI_dat_Angles,4)
        TDI_dat_Angles(:,:,:,i) = convn(squeeze(TDI_dat_Angles(:,:,:,i)),mk3Dwin,'same');
    end
end

% From phase differences to displacements or TDI, sign-convention: positive means in the direction of the probe displacements
TDI_Angles= angle((TDI_dat_Angles))./pi*(SoS*PRF)/(4*f_center); % [m/s]
TDI = squeeze(TDI_Angles);

clear TDI_Angles TDI_dat_Angles TDICart Frames
TDICart = double(TDI);
clear  TDI

%% Bandpass Filtering TDI data
taxis_long = taxis; 
taxis = taxis(t_ind);
numFrames_long = numFrames;
numFrames = length(taxis);

if BANDPASS_FILT
    bp = designfilt('bandpassfir','FilterOrder',3,'CutoffFrequency1',50,'CutoffFrequency2',1000,'SampleRate',PRF); % bandpass filter
    TDI_temp = permute(TDICart(:,:,t_ind,:),[3 1 2 4]);
    TDI_filt = filtfilt(bp,TDI_temp);
    TDI_filt = permute(TDI_filt,[2 3 1 4]);
    if SHOW_ITERIM_STEPS
        % Example of the effect of the bandpass filter
        figure('InnerPosition', [50 100 300*2.5 250*2],'PaperPositionMode','auto')
        tcl = tiledlayout(2,1,TileSpacing="tight",Padding="tight");
        nexttile
        hold on
        plot(taxis,squeeze(TDICart(100,50,t_ind,1)),LineWidth=2)
        plot(taxis,squeeze(TDI_filt(100,50,:,1)),LineWidth=2)
        ylabel('v_{axial} (m/s)',FontSize=13)
        legend('Unfiltered','Bandpass')
        title(sprintf('z = %.1f mm, x = %.1f mm',Z_IQ(100),X_IQ(50)),FontSize=14)
        nexttile
        hold on 
        plot(taxis,squeeze(TDICart(100,100,t_ind,1)),LineWidth=2)
        plot(taxis,squeeze(TDI_filt(100,100,:,1)),LineWidth=2)
        xlabel('t (s)',FontSize=13)
        ylabel('v_{axial} (m/s)',FontSize=13)
        title(sprintf('z = %.1f mm, x = %.1f mm',Z_IQ(100),X_IQ(100)),FontSize=14)
        title(tcl,'Effect of bandpass in two points',FontSize=16)
    end
else
    TDI_filt = TDICart(:,:,t_ind,:);
end 

clear TDICart TDICart;


%% Temporal upsampling
if TEMP_UPSAMP
    interp_factor = 5;
    t_ind = t_ind(1)*interp_factor:t_ind(end-1)*interp_factor;
    taxis_before = taxis;
    taxis = interp(taxis(1:end-1), interp_factor);
    dt = mean(diff(taxis));
    PRF = 1/dt;
    numFrames = length(taxis);
    numref = numref*interp_factor;
    TDI_filt = permute(TDI_filt, [3 1 2 4]);
    TDI_filt = interp1(taxis_before, TDI_filt, taxis,'spline');
    TDI_filt = permute(TDI_filt, [2 3 1 4]);
else 
    interp_factor = 1;
end

%% Directional Filtering
% Based on: Lipman et al., IEEE TUFFC, 63 (8), 2016. 

if DIRECTIONAL_FILT 
    %% Fourier transform
	N = 2^9; % nearest power of 2 for all dimensions
	TDI_fft = zeros(N,N,N,2);			 
    for i = PUSH_SELECTION
        TDI_fft(:,:,:,i) = fftn(TDI_filt(:,:,:,i), [N N N]);
        TDI_fft(:,:,:,i) = fftshift(TDI_fft(:,:,:,i));
    end 
    f = linspace(0,PRF,size(TDI_fft,3))-PRF/2; % temporal frequencies
    kx = linspace(0,1/dX,size(TDI_fft,2))-1/dX/2; % lateral spatial frequencies
    kz = linspace(0,1/dZ,size(TDI_fft,1))-1/dZ/2; % axial spatial frequencies
    [KX, KZ, F] = meshgrid(kx,kz,f);

    if SHOW_ITERIM_STEPS
        % Showing tissue displacements in Fourier domain 
        push = 1;
        figure('InnerPosition', [50 100 300*4 250*1.5],'PaperPositionMode','auto')
        tcl = tiledlayout(1,4,TileSpacing="tight",Padding="tight");
        nexttile
        temp_max = max(abs(TDI_fft(:,:,:,push)),[],'all');
        h = slice(KX,KZ,F,abs(TDI_fft(:,:,:,push)/temp_max),0,0,0);
        xlabel('k_x [1/mm]',Fontsize=13); ylabel('k_z [1/mm]',Fontsize=13); zlabel('f [Hz]',Fontsize=13);
        set(h,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp')
        alpha('color') 

        nexttile
        h = slice(KX,KZ,F,abs(TDI_fft(:,:,:,push))/temp_max,0,0,0);
        view(0,0)
        axis square
        xlabel('k_x [1/mm]',Fontsize=13); ylabel('k_z [1/mm]',Fontsize=13); zlabel('f [Hz]',Fontsize=13);
        set(h,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp')
        alpha('color') 
        xlim([-0.5 0.5]); ylim([-0.5 0.5]); zlim([-1000 1000])
        zticks([-1000 0 1000]); xticks([-0.5 0 0.5]); yticks([-0.5 0 0.5])
        box on
        
        nexttile
        h = slice(KX,KZ,F,abs(TDI_fft(:,:,:,push))/temp_max,0,0,0);
        view(90,0)
        axis square
        xlabel('k_x [1/mm]',Fontsize=13); ylabel('k_z [1/mm]',Fontsize=13); zlabel('f [Hz]',Fontsize=13);
        set(h,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp')
        alpha('color') 
        xlim([-0.5 0.5]); ylim([-0.5 0.5]); zlim([-1000 1000])
        zticks([-1000 0 1000]); xticks([-0.5 0 0.5]); yticks([-0.5 0 0.5]);
        box on

        nexttile
        h = slice(KX,KZ,F,abs(TDI_fft(:,:,:,push))/temp_max,0,0,0);
        view(0,90)
        axis square
        xlabel('k_x [1/mm]',Fontsize=13); ylabel('k_z [1/mm]',Fontsize=13); zlabel('f [Hz]',Fontsize=13);
        set(h,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp')
        alpha('color') 
        xlim([-0.5 0.5]); ylim([-0.5 0.5]); zlim([-1000 1000])
        zticks([-1000 0 1000]); xticks([-0.5 0 0.5]); yticks([-0.5 0 0.5]);
        colorbar
        box on

        title(tcl,sprintf('Tissue displacements in Fourier domain (push %.0f) before directional filtering',push),FontSize=14)
    end

    %% First filter component (spatial variation)
    u1 = [0 1]; % [z x] Propagating in +x direction
    u2 = [0 -1]; % [z x] Propagating in -x direction
    q = 2; % typically 2 or 3
    
    filter11 = zeros(length(kz),length(kx));
    filter12 = zeros(length(kz),length(kx));
    for kxi = 1:length(kx)
        for kzi = 1:length(kz)
            k = [kz(kzi); kx(kxi)];
            if k == [0; 0] 
                filter11(kzi,kxi) = 1;
                filter12(kzi,kxi) = 1;
            else
                filter11(kzi,kxi) = (u1*k/(norm(u1)*norm(k)))^q;
                filter12(kzi,kxi) = (u2*k/(norm(u2)*norm(k)))^q;
            end 
        end
    end 

    if SHOW_ITERIM_STEPS
        % Show the first filter component
        figure('Position', [50 100 figset1 figset2],'PaperPositionMode','auto')
        tiledlayout(1,1,TileSpacing="tight",Padding="tight");
        nexttile
        imagesc(kx,kz,filter11)
        axis square
        xlabel('k_x [1/mm]',FontSize=13)
        ylabel('k_z [1/mm]',FontSize=13)
        title('Spatial filter',FontSize=16)
        colormap('gray')
        colorbar
        clim([0 1])
    end 
    filter11 = repmat(filter11,[1 1 length(f)]);
    filter12 = repmat(filter12,[1 1 length(f)]);

    %% Second filter component (space-time)
    if rem(length(f),2)
        midf = ceil(length(f)/2);
    else 
        midf = length(f)/2:length(f)/2+1;
    end 
    if rem(length(kx),2)
        midx = ceil(length(kx)/2);
    else 
        midx = length(kx)/2:length(kx)/2+1;
    end 
    filter21 = double(KX.*F<=0);
    filter21(:,midx,:) = 0.5;
    filter21(:,:,midf) = 0.5;
    filter22 = 1-filter21;

    if SHOW_ITERIM_STEPS
        % Show second filter component
        figure('Position', [50 100 figset1*2 figset2],'PaperPositionMode','auto')
        tiledlayout(1,2,TileSpacing="tight",Padding="tight");
        nexttile
        imagesc(kx,f,squeeze(filter21(1,:,:))')
        axis square
        xlabel('k_x [1/mm]',Fontsize=13)
        ylabel('f [Hz]',Fontsize=13)
        title('Space-time filter (push 1)',Fontsize=16)
        colormap gray
        nexttile 
        imagesc(kx,f,squeeze(filter22(1,:,:))')
        axis square
        xlabel('k_x [1/mm]',Fontsize=13)
        ylabel('f [Hz]',Fontsize=13)
        title('Space-time filter (push 2)',Fontsize=16)
        colormap gray
        colorbar;
    end

    %% Complete directional filter
    filter = zeros(size(TDI_fft));
    filter(:,:,:,1) = filter11.*filter21;
    filter(:,:,:,2) = filter12.*filter22;
    
    if SHOW_ITERIM_STEPS
        % Show the directional filter
        figure('InnerPosition', [50 100 600 500],'PaperPositionMode','auto')
        h = slice(KX,KZ,F,filter(:,:,:,1),[],[],f(1:5:end));
        xlabel('k_x [1/mm]',FontSize=13)
        ylabel('k_z [1/mm]',FontSize=13)
        zlabel('f [Hz]',FontSize=13)
        title('Directional filter',FontSize=14)
        set(h,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp')
        alpha('color') 
        alphamap('rampup')
        alphamap('increase',.05)
        colormap('turbo')
        colorbar
    end

    %% Filter data & Inverse Fourier transform
    TDI_fft_filt = TDI_fft.*filter;
	TDI_filt2 = zeros(size(TDI_filt));				 
    for i = PUSH_SELECTION
        if SHOW_ITERIM_STEPS
            figure('InnerPosition', [50 100 300*4 250*1.5],'PaperPositionMode','auto')
            tcl = tiledlayout(1,4,TileSpacing="tight",Padding="tight");
            nexttile
            temp_max = max(abs(TDI_fft_filt(:,:,:,i)),[],'all');
            h = slice(KX,KZ,F,abs(TDI_fft(:,:,:,i)/temp_max),0,0,0);
            xlabel('k_x [1/mm]',Fontsize=13); ylabel('k_z [1/mm]',Fontsize=13); zlabel('f [Hz]',Fontsize=13);
            set(h,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp')
            alpha('color') 
    
            nexttile
            h = slice(KX,KZ,F,abs(TDI_fft_filt(:,:,:,i))/temp_max,0,0,0);
            view(0,0)
            axis square
            xlabel('k_x [1/mm]',Fontsize=13); ylabel('k_z [1/mm]',Fontsize=13); zlabel('f [Hz]',Fontsize=13);
            set(h,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp')
            alpha('color') 
            xlim([-0.5 0.5]); ylim([-0.5 0.5]); zlim([-1000 1000])
            zticks([-1000 0 1000]); xticks([-0.5 0 0.5]); yticks([-0.5 0 0.5]);
            box on
            
            nexttile
            h = slice(KX,KZ,F,abs(TDI_fft_filt(:,:,:,i))/temp_max,0,0,0);
            view(90,0)
            axis square
            xlabel('k_x [1/mm]',Fontsize=13); ylabel('k_z [1/mm]',Fontsize=13); zlabel('f [Hz]',Fontsize=13);   
            set(h,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp')
            alpha('color') 
            xlim([-0.5 0.5]); ylim([-0.5 0.5]); zlim([-1000 1000])
            zticks([-1000 0 1000]); xticks([-0.5 0 0.5]); yticks([-0.5 0 0.5]);
            box on
    
            nexttile
            h = slice(KX,KZ,F,abs(TDI_fft_filt(:,:,:,i))/temp_max,0,0,0);
            view(0,90)
            axis square
            xlabel('k_x [1/mm]',Fontsize=13); ylabel('k_z [1/mm]',Fontsize=13); zlabel('f [Hz]',Fontsize=13);
            set(h,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp')
            alpha('color') 
            xlim([-0.5 0.5]); ylim([-0.5 0.5]); zlim([-1000 1000])
            zticks([-1000 0 1000]); xticks([-0.5 0 0.5]); yticks([-0.5 0 0.5]);
            colorbar
            box on
            title(tcl,sprintf('Tissue displacements in Fourier domain (push %.0f) after directional filtering',i),FontSize=14)
        end
        TDI_fft_filt(:,:,:,i) = ifftshift(squeeze(TDI_fft_filt(:,:,:,i)));
        temp(:,:,:,i) = real(ifftn(TDI_fft_filt(:,:,:,i)));
		TDI_filt2(:,:,:,i) = temp(1:size(TDI_filt,1),1:size(TDI_filt,2),1:size(TDI_filt,3),i);
    end 

    TDI_filt_old = TDI_filt;
    TDI_filt = TDI_filt2;
    clear TDI_filt2
end 

%% Showing TDI data
if SHOW_TDI
    gcf4 = figure('Color','white','Position', [50 20 figset1 figset2],'PaperPositionMode','auto');
    for j = 1:size(TDI_filt,4)
        if DIRECTIONAL_FILT
            pushdir = [PROC_DAT_DIR '\Push'  num2str(j) '_DirFilt'];
        else 
            pushdir = [PROC_DAT_DIR '\Push'  num2str(j) '_NoDirFilt'];
        end 
        if SAVE_MOVIE_TDI
            if exist(pushdir,'dir'); rmdir(pushdir,'s'); end
            if ~exist(pushdir,'dir'); mkdir(pushdir); end
            name_TDImovie = [pushdir '\TDI_push' num2str(j) '.gif'];
        end 

        for kk = 1:interp_factor:size(TDI_filt,3)
            imagesc(X_IQ,Z_IQ,squeeze(TDI_filt(:,:,kk,j)));
            axis image
            hold on
            plot(x_circle,z_circle,Color=0.5*[1 1 1],LineWidth=0.5)         
            xlabel('x [mm]',Fontsize=13);ylabel('z [mm]',Fontsize=13);
            cbar = colorbar('Location','eastoutside');
            clim(TDIcolorscale)
            colormap('parula')
            ylabel(cbar, 'Axial velocities [m/s]',Fontsize=13);
            title([num2str((taxis(kk))*10^3,'%.1f') ' ms'],Fontsize=14)
            pause(0.005)

            if SAVE_MOVIE_TDI
                frame = getframe(gcf4);
                im = frame2im(frame);
                [A,map] = rgb2ind(im,256);

                if kk == 1
                    imwrite(A,map,name_TDImovie,'gif','LoopCount',Inf,'DelayTime',0.1);
                else
                    imwrite(A,map,name_TDImovie,'gif','WriteMode','append','DelayTime',0.1);
                end
            end
            print('-dpng','-r200',[pushdir '\SWpropagation_frame' num2str(kk) '.png']);
        end
    end
end
toc

Bmode = squeeze(FramesCart(:,:,2,:))/max_Frames;
clear FramesCart

%% Amplitude and Push Mask
% Amplitude mask: consider region close to first push
ampl_tot = prctile(TDI_filt(:,:,1,:),0.5,'all')*0.05; % find minimal value (w/o outliers)
ampl_mask = squeeze(min(TDI_filt,[],3)) < ampl_tot;

% Push Mask: consider region with significant wave propagation 
margin_p = ceil(margin_mm/dX);

% Push mask based on given location of push
push_mask1 = Xaxis>=(X_push(1)+margin_mm);
push_mask2 = Xaxis<=(X_push(2)-margin_mm);

% Adapting push mask 
for iz = 1:size(TDI_filt,1)
    [~, locs1] = findpeaks(-TDI_filt(iz,:,1,1),MinPeakProminence=1e-3);
    [~, locs2] = findpeaks(-TDI_filt(iz,:,1,2),MinPeakProminence=1e-3);
    
    if isempty(locs1) 
        locs1 = 1; 
    elseif locs1(end)+margin_p>size(TDI_filt,2)
        locs1 = locs1(end);
    else 
        locs1 = locs1(end)+margin_p;
    end
    
    if isempty(locs2) 
        locs2 = size(TDI_filt,2); 
    elseif locs2(1)-margin_p<1
        locs2 = locs2(1);
    else 
        locs2 = locs2(1)-margin_p;
    end

    push_mask1(iz,1:locs1) = 0;
    push_mask2(iz,locs2:end) = 0;
end
push_mask = cat(3,push_mask1,push_mask2); % concatenate
tot_mask = push_mask.*ampl_mask;

if SHOW_ITERIM_STEPS
    % Show amplitude and push masks for selected push
    push = 2;
    figure('Position', [50 100 figset1*3 figset2],'PaperPositionMode','auto')
    tcl = tiledlayout(1,4,TileSpacing="tight",Padding="tight");
    
    n1 = nexttile;
    imagesc(X_IQ,Z_IQ,squeeze(TDI_filt(:,:,1,push)))
    axis image
    xlabel('x [mm]',FontSize=13); ylabel('z [mm]',FontSize=13);
    colorbar
    title(sprintf('TDI - %.1f ms',taxis(1)*1e3),FontSize=14)
    colormap(n1,"parula")

    n2 = nexttile;
    imagesc(X_IQ,Z_IQ,ampl_mask(:,:,push))
    axis image
    xlabel('x [mm]',FontSize=13); ylabel('z [mm]',FontSize=13);
    title('Amplitude mask',FontSize=14)
    colormap(n2,"gray")

    n3 = nexttile;
    imagesc(X_IQ,Z_IQ,push_mask(:,:,push))
    axis image
    xlabel('x [mm]',FontSize=13); ylabel('z [mm]',FontSize=13);
    title('Push mask',FontSize=14)
    colormap(n3,"gray")

    n4 = nexttile;
    imagesc(X_IQ,Z_IQ,tot_mask(:,:,push))
    axis image
    xlabel('x [mm]',Fontsize=13); ylabel('z [mm]',Fontsize=13);
    title('Total mask',FontSize=14)
    colormap(n4,"gray")
    cb = colorbar;
    cb.Ticks = [0 1];
    
    title(tcl,sprintf('Push %.0f',push),FontSize=16)
end 

%% Elastogram
% Based on: Song et al., Ultras. Med. Biol., 40 (6), 2014.
% Algorithm considers a 2D window

h = round(w-p)/2;
Lac = p*dX*1e-3;
Lab = p*dZ*1e-3;
maxlagx = ceil(p*dX*1e-3/(Vmin*dt)); % max lag considering minimal velocity
maxlagz = ceil(p*dZ*1e-3/(Vmin*dt));
Vmap = zeros(size(TDI_filt,1),size(TDI_filt,2),size(TDI_filt,4));

tic
for push = PUSH_SELECTION
    for m = w/2+h+1:size(TDI_filt,2)-w/2-h % in lateral (x) dimension
        for n = w/2+h+1:size(TDI_filt,1)-w/2-h % in axial (z) dimension
            if tot_mask(n,m,push) 
                i = m-h:m+h;
                j = n-h:n+h;
                CCx_p = zeros(2*h+1,2*h+1);
                CCz_p = zeros(2*h+1,2*h+1);
                tac_p = zeros(2*h+1,2*h+1);
                tab_p = zeros(2*h+1,2*h+1);
                r_ij = ones(2*h+1,2*h+1);
                for ii = 1:length(i)
                    for jj = 1:length(j)
                        % Reciprocal of distance to center pixel of window
                        if ~(i(ii)==m && j(jj)==n)
                            r_ij(ii,jj) = sqrt((i(ii)-m)^2+(j(jj)-n)^2);
                        end 
                        % Calculate timelag and correlation
                        [CCx,lagx] = xcorr(squeeze(TDI_filt(j(jj),i(ii)+p/2,:,push)),squeeze(TDI_filt(j(jj),i(ii)-p/2,:,push)),'normalized',maxlagx);
                        [CCz,lagz] = xcorr(squeeze(TDI_filt(j(jj)+p/2,i(ii),:,push)),squeeze(TDI_filt(j(jj)-p/2,i(ii),:,push)),'normalized',maxlagz);
                        [CCx_temp, tdx] = max(CCx); 
                        [CCz_temp, tdz] = max(CCz);
                        CCx_p(ii,jj) = CCx_temp; 
                        CCz_p(ii,jj) = CCz_temp; 
                        tac_p(ii,jj) = lagx(tdx)/PRF;
                        tab_p(ii,jj) = lagz(tdz)/PRF;
                    end 
                end 
                % Normalize the correlations
                small_correlation = 0;
                if mean(CCx_p,'all')<0.5
                    small_correlation = 1;
                end 
                CCx_p = (CCx_p.^2./r_ij)./sum(CCx_p.^2./r_ij,'all'); 
                CCz_p = (CCz_p.^2./r_ij)./sum(CCz_p.^2./r_ij,'all');
                % Calculate velocity
                Vmap(n,m,push) = sum(Lac*Lab*CCx_p.*CCz_p./sqrt((Lac*tab_p.*CCx_p).^2+(Lab*tac_p.*CCz_p).^2),'all');
                if abs(Vmap(n,m,push)) > Vmax || small_correlation
                    Vmap(n,m,push) = NaN;
                end
            else 
                Vmap(n,m,push) = NaN;
            end
        end
    end
end 

% Combining SWV maps of pushes
elastogram = mean(Vmap,3,'omitnan');
elastogram = fillmissing(elastogram,'knn');
if SMOOTHING_MAP; elastogram = medfilt2(elastogram,[5 5]); end
    
% Quality metrics
SWV_ROI1 = median(elastogram(ROI1),'all','omitnan');
SWV_ROI2 = median(elastogram(ROI2),'all','omitnan');
IQR_ROI1 = iqr(elastogram(ROI1),'all');
IQR_ROI2 = iqr(elastogram(ROI2),'all');
metrics_txt = sprintf('I: %.2f $\\pm$ %.2f m/s, B: %.2f $\\pm$ %.2f m/s',SWV_ROI1,IQR_ROI1,SWV_ROI2,IQR_ROI2);
if phantom==1; metrics_txt = sprintf('%.2f $\\pm$ %.2f m/s',SWV_ROI2,IQR_ROI2); end 

%% Plotting SWV maps
figure('Position', [50 100 figset1*2.7 figset2],'PaperPositionMode','auto')
tcl = tiledlayout(1,3,TileSpacing="tight",Padding="tight");
nexttile 
imagesc(X_IQ,Z_IQ,Vmap(:,:,1),'AlphaData',~isnan(Vmap(:,:,1)))
axis image
title('Push 1',Fontsize=13)
subtitle(' ',Fontsize=13)
xlabel('x [mm]',Fontsize=13); ylabel('z [mm]',Fontsize=13);
clim(SWVcolorscale)

nexttile 
imagesc(X_IQ,Z_IQ,Vmap(:,:,2),'AlphaData',~isnan(Vmap(:,:,2)))
axis image
title('Push 2',Fontsize=13)
subtitle(' ',Fontsize=13)
xlabel('x [mm]',Fontsize=13); ylabel('z [mm]',Fontsize=13);
clim(SWVcolorscale)

nexttile 
imagesc(X_IQ,Z_IQ,elastogram)
hold on
plot(x_circle,z_circle,Color=0.5*[1 1 1],LineWidth=0.5)
axis image
title('SWV map',Fontsize=13)
subtitle(metrics_txt,Fontsize=13,Interpreter="latex")
xlabel('x [mm]',Fontsize=13); ylabel('z [mm]',Fontsize=13);
cbar = colorbar;
clim(SWVcolorscale)
ylabel(cbar, 'SWV [m/s]',Fontsize=13);
cbar.Layout.Tile = "east";
toc

% Saving results
if DIRECTIONAL_FILT
    elastogramdir = 'Elastogram_DirFilt';
    SWVsPushes_dir = 'SWVsPushes_DirFilt';
else 
    elastogramdir = 'Elastogram3_NoDirFilt';
    SWVsPushes_dir = 'SWVsPushes_NoDirFilt';
end 
print('-dpng','-r200',[PROC_DAT_DIR '\' elastogramdir '.png']);
save([PROC_DAT_DIR '\' file(1:end-4) '_' elastogramdir '.mat'],'elastogram');
save([PROC_DAT_DIR '\' file(1:end-4) '_' SWVsPushes_dir '.mat'],'elastogram');