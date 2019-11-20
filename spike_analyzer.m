function [nuis,cfg] = spike_analyzer(P,rpfile,cfg)
%
% USAGE: [nuis,cfg] = spike_analyzer(P,rpfile,cfg);
%
% ARGS:
% P       -- cell array with preprocessed EPI file names (use the set of
%            images that are included in the statistical analysis), e.g. 
%            selected with P = cellstr(spm_select(inf,'image')
% rpfile  -- text file with realignment parameters
% cfg     -- struct with configuration options with the following fields:
%  .voxselect -- method of voxel selection, with the following options:
%                'globmean':  similar to the algorithm in spm_global
%                             (faster, but prone to outliers)
%                'brainmask': use spm_segment to compute a brainmask for
%                             the EPI images, see cfg.bmask below for
%                             futher options (slow,but recommended)
%                'roi':       specify an anatomical ROI image and select
%                             all voxels within the ROI
%                             (for evaluating spikes in a target ROI only)
%  .roiimg  -- full path name to an anatomical ROI image [only used if
%              cfg.voxselect = 'roi'], make sure that the EPI series is
%              warped to the standard space in which the ROI is defined and
%              overlay the ROI image on the EPI using CheckReg and Blobs
%  .bmask   -- use either the gray or withe matter segment only or the
%              complete brain mask ('gray','white','both') [only used if
%              cfg.voxselect = 'brainmask']
%  .globthr -- global mean threshold in SD (defined on 1st derivative)
%              (default 0.75)
%  .winsize -- size of sliding window to predict global mean from movement
%              parameters (default: 20 images)
%  .movethr -- threshold in mm/TR for fast movement detection (default 0.2)
%  .ralignthr-- threshold for translational and rotational movement
%              parameters (conservative: transl. = 1*voxel size, rot. = 1
%              degree, liberal: transl. = 2*voxel size, rot. = 2 degrees)
%              (default = 1)
%  .rsqrthr -- threshold of R-Suare between movement params and global mean
%              (default 0.8)
%  .pflag   -- if true, the Graphics window will be plotted as Postscript to 
%              the file spike_analyzer.ps in the input directory (default:
%              false)
%
% OUTPUT:
% nuis    -- struct with matrices of dummy variables based on the differen
%            thresholds
%            .gm      -- global mean threshold
%            .m       -- movement velocity threshold
%            .rsqr    -- rsquare threshold
%            .realign -- reliagnment threshold
%            .all     -- all combined
% cfg     -- struct with configuration options

% =========================================================================
% RATIONALE 
%
% Spike Analyzer is a program to help you to detect (unnatural) spike
% (outliers) in an EPI time series. It is NOT a replacement for visual
% inspection, but it'll help you in selecting certain suspicious scan for
% closer visual inspection, e.g. with spm_movie.
%
% Spikes are detected in two ways: (1) by thresholding the 1st derivative
% of the global activity of an EPI series or (2) by thresholding the
% movement velocity. All scans that exceed these thresholds are marked as
% suspicious (candidates for removal by nuisance regressors). Images
% detected by (1) are shown as vertical bars in light blue in the plots,
% images detected by (2) are shown in light red.
%
% Visually you can see spike as "jumps" in the global activity, or as
% sudden jerks in the 1st derivative. The 1st derivative is better suited
% for defining a threshold because except for some jerks, this curves will
% be close to zero, even when the original global mean curve settles on a
% different level after some jumps.
%
% These "jumps" can be caused e.g. by scanner malfunction, but also by
% subject movement. In order evaluate the latter cause, Spike Analyzer uses
% the subject movement parameters to predict the global mean in a defined
% windows of scans (defaults: 20). It then moves through the entire time
% series with the same GLM. R-Square are color-coded and are used to
% evaluate if some jumps in the global mean could have been caused by
% movement (hot colorbar from dark red - orange - yellow - white). If the
% jumps in the global mean curve coincide with a hot color in the R-Square
% bar (named: Correlation between global activity and movement parameters)
% then there is a good chance that it was cause by movement.
% 
% Movement parameters are also shown to evaluate, if the subject has shown
% extraordinary large movement (rule of thumb: exclude subject if
% tranlational movement exceed 2x voxel size or 1 degree of rotational
% movement. These are independent of whether the movement is correlated
% with the global mean.
%
% The global mean is also shown with effects of movement removed (again
% through a GLM). Jumps in the global mean curve that are correlated
% with movement should not be visible in the last plot
%
% Finally, Spike Analyzer will also shown you how many images are marked as
% suspicious. You can then decide to adjust the thresholds (cfg.globthr and
% cfg.movethr) in order to exclude or include more scans.
%
% Spike Analyzer will output matrices of dummy regressors, in which a
% suspicious scan is marked with 1 and all other scans are marked with
% zero. The 4 matrices are for global mean based nuisance regressors
% (gmNuis), movement velocity based nuisance regressors (mNuis), the
% correlation-based neuisance regressors (rsqrNuis), or all combined
% (allNuis). To include these nuisance regressors in your design 
% matrix you need to have to enter them in SPM.Sess(i).C during design
% specification, e.g.
% SPM.Sess(1).C.C = nuis.all;
% SPM.Sess(1).C.name = repmat({'BS'},1,size(SPM.Sess(1).C.C,2));
% (since one is usually not interested in looking at the effect of a single
% nuisance regressor, they are all assigned the same name (BS = bad scan))
% =========================================================================
% EXPLANATION OF PLOTS (from top to bottom)
%
% 1. Normalized Global activity (blue) and 1st derivative of global
%    activity (gray) 
%    Configuration Options:
%    cfg.globthr  = threshold (in SD) for excluding scans that exceed this
%                   threshold from the 1st derivative of global mean
%                   (default: 0.5) 
%
% 2. small bar with R-square values of predicting the normalized global
%    mean from the movement parameters
%    cfg.rsqrthr  = threshold (range 0-1) on R-Square value between
%                   movement parameters as a predictor of global mean
%                   (default: 0.8)
%
% 3. translational movement parameters
%    cfg.ralignthr = threshold (1 or 2) for translational and rotational
%    parameters (rule of thumb). Conservative: translational threshold
%    1 * voxel size, rotational threshold: 1 degree. Liberal: 2* voxel
%    size, 2 degrees
%
% 4. rotational movement parameters
%
% 5. Fast movement detection (using the algorithm of ArtRepair Toolbox)
%    Configuration Options:
%    cfg.movethr = threshold for excluding scan based on fast movement
%                  (mm/TR) (default: 0.2)
%
% 6. Normalized Global Activity with effects of movement parameters removed
%    (by multiple linear regression)
% 
% =========================================================================
% SOME TECHNICAL INFORMATION
% Spike Analyzer used to compute the global mean exclusively using 
% spm_global. However, several test and simulations have revealed that 
% spm_global is prone to include outlying voxels, which tend to distort the
% global mean, thus leading to the potential exclusion of too many "bad
% scans". spm_global computes the global mean on an image-by-image basis,
% which almost certainly leads to an inclusion of different voxels for each
% image. This is not optimal for the computation of a global mean time
% series. Therefore, Spike Analyzer does not rely on spm_global anymore and
% offers 3 options for selection the voxels to compute the global mean:
% 1. Segment the first EPI image of the series using spm_segment and
%    computing a brain mask from the gray and white matter segments.
%    Optionally, the gray and white matter segments can be used exclusively
%    for voxel selection. This option is slow, but provides the most
%    accurate voxel selection
% 2. Computing the global mean on a constant number of voxels. This
%    requires several passes through the data: (a) compute and
%    image-by-image threshold by dividing the mean of each image by 8 (like
%    in the spm_global), (b) setting the overall threshold to the mean of
%    the image-spcific thresholds, (c) finding the voxels in each image
%    that exceed the overall threshold, and (d) selecting only those voxels
%    common to all images that exceed the overall threshold. Although this
%    seems a bit complicated, the results of this method are much better
%    than the plain spm_global computations and are similar to the
%    brainmask approach in 1.
% 3. Compute the global mean only in an anatomical ROI. If you are
%    interested in a specific brain region, then you can configure an
%    anatomical ROI image and only those voxels in the ROI will be selected
%    for the computation of global mean activity.
% 4. For the sake of comparison, cfg.voxselect can be set to 'spm_global',
%    which will then use the old spm_global approach. This method is
%    depreciated.
%
% The spm_global problem of selecting different voxels per image is
% especially aggravated for images that are not preprocessed yet. It is
% therefore recommended to use Spike Analyzer only on preprocessed images
% that will enter the statistical analysis (i.e. realign and warped
% images). The nuisance regressor created by Spike Analyzer will then fit
% to the images that will be eneter into the 1st level analysis, which is
% preferable. Happy Analyzing!
% =========================================================================
% (c) by Jan Glaescher (10/11/2012)

version = '1.6';

% === SETUP
spm('FnBanner',mfilename,version);
[junk,Fgraph] = spm('FnUIsetup','Spike Analyzer');
%spm_help('!ContextHelp',mfilename);

% Parse or prompt for user input
if nargin < 1
  str = sprintf('Select preprocessed EPI time series');
  P = spm_select(Inf,'image',str);
end

if nargin < 2
  str = sprintf('Select text file with movement parameters');
  rpfile = spm_select(1,'^rp.*\.txt',str);
end

if nargin < 3 
  cfg.voxselect = 'brainmask'; % or 'globmean' or 'roi'
  cfg.roiimg    = ''; % path to ROI image
  cfg.bmask     = 'both'; %use brainmask of gray and white matter
  cfg.globthr   = 0.75; % bad scan if activity exceeds +- 0.5 SD
  cfg.winsize   = 20; % window size for predicting global act. from movemt
  cfg.movethr   = 0.2; % bad scan if mm/TR > 0.2
  cfg.rsqrthr   = 0.8; % bad scan if R-Square(movement,globmean) > 0.8
  cfg.realignthr = 1;
  cfg.pflag     = 0; % do not print to psotscript file
end

blue    = [0.1211    0.4648    0.7031];
orange  = [0.9961    0.4961    0.0586];
red     = [0.8359    0.1523    0.1562];
green   = [0.1719    0.6250    0.1719];
violet  = [0.5781    0.4023    0.7383];
yellow  = [1.0000    0.7500         0];
magenta = [0.8945    0.2578    0.4570];
brown   = [0.4141    0.2031    0.1680];

% setup output files
[pa,nam,ext] = fileparts(P(1,:)); %location of images and rp_*.txt file
%psfile  = fullfile(pa,'spike_analyzer.ps');
psfile  = fullfile(pwd,'spike_analyzer.ps');

% load files
rp   = spm_load(rpfile);
nimg = size(P,1);
xl   = [0 nimg+1];

% map image files into memory
for f=1:length(P)
  Vin(f) 	= spm_vol(deblank(P(f,:)));
end

voxsize = min(abs(diag(Vin(1).mat(1:3,1:3))));

gx = zeros(1,nimg);
data = spm_read_vols(Vin);

switch cfg.voxselect
  case 'globmean' % I'll do my own thing there and don't rely on spm_global

    for f = 1:nimg
      d   = data(:,:,:,f);
      tmp(f) = mean(d(:))./8; % image average/8 (as in spm_global)
    end
    thr = mean(tmp); % use mean threshold for voxel selection
    for f = 1:nimg
      d   = data(:,:,:,f);
      i(:,f) = d(:)>thr; % find voxels exceed the mena threshold
    end
    ind = sum(i,2)==nimg; % use only above threshold voxels
    for f = 1:nimg
      d = data(:,:,:,f);
      gx(f) = mean(d(ind));
    end
          
  case 'brainmask'
    if ~exist(fullfile(pa,strcat('c1',nam,ext)),'file')
      pth = fileparts(P(1,:));
      segment_epi(P(1,:),pth);
    end
    cfg.roiimg = fullfile(pth,'epi_brainmask.nii');
    switch cfg.bmask
      case 'both'
        Vmask = spm_vol(fullfile(pa,'epi_brainmask.nii'));
        thr = 0;
      case 'gray'
        Vmask = spm_vol(fullfile(pa,strcat('c1',nam,ext)));
        thr = 0.2;
      case 'white'
        Vmask = spm_vol(fullfile(pa,strcat('c2',nam,ext)));
        thr = 0.2;
    end
    mask  = spm_read_vols(Vmask);
    ind  = mask(:)>thr;
    for f = 1:nimg
      d = data(:,:,:,f);
      gx(f) = mean(d(ind));
    end
    
  case 'roi'
    % reslice ROI to dimensions of EPI image
    Vroi = spm_vol(cfg.roiimg);
    d = zeros(Vin(1).dim);
    for p = 1:Vin(1).dim(3)
      B = spm_matrix([0 0 -p 0 0 0 1 1 1]);
      M = inv(B*inv(Vin(1).mat)*Vroi.mat);
      d(:,:,p) = spm_slice_vol(Vroi,M,Vin(1).dim(1:2),1);
    end
    ind = d(:)>0;
    for f = 1:nimg
      d = data(:,:,:,f);
      gx(f) = mean(d(ind));
    end
  
  case 'spm_global'
    for f=1:nimg
      gx(f) = spm_global(Vin(f));
    end
end

% === COMPUTATIONS ===


% first derivative
ngx = zscore(gx);
dngx = [0 diff(ngx)];

% find images exceeding gthr
gmX  = find(dngx>=cfg.globthr|dngx<=-cfg.globthr);

if ~isempty(gmX)
  for f=1:length(gmX)
    gmPx{f} = [gmX(f)-0.5 gmX(f)+0.5 gmX(f)+0.5 gmX(f)-0.5];
  end
end

% compute velocity
% code taken from ArtRepair5 (but implemented more efficiently)
rprev = [rp(end:-1:1,:); zeros(1,6)];
drp = diff(rprev).^2;
drp(:,4:6) = drp(:,4:6)*1.28;
mmtr = sqrt(sum(drp(end:-1:1,:),2));

% find images exceeding mthr
mX = find(mmtr>=cfg.movethr);

if ~isempty(mX)
  for f=1:length(mX)
    mPx{f} = [mX(f)-0.5 mX(f)+0.5 mX(f)+0.5 mX(f)-0.5];
  end
end

% compute movement parameter thresholds
diffrp = [zeros(1,6); diff(rp)];
mvX = [];
for f = 1:6
  mvX = [mvX; find(diffrp(:,f)>=cfg.realignthr)];
end
mvX = unique(mvX);

if ~isempty(mvX)
  for f=1:length(mvX)
    mvPx{f} = [mvX(f)-0.5 mvX(f)+0.5 mvX(f)+0.5 mvX(f)-0.5];
  end
end


% predict global mean from movement parameters
for f=1:nimg
  window = f-(cfg.winsize/2):f+(cfg.winsize/2);
  window(window<=0) = [];
  window(window>nimg) = [];
  y = ngx(window)';
  x = [rp(window,:) ones(length(window),1)];
  b = pinv(x)*y;
  res = y-(x*b);
  rss = res'*res;
  tss = (y-mean(y))'*(y-mean(y));
  rsq(f) = 1-(rss/tss);
end

rsqrX = find(rsq>=cfg.rsqrthr);
if ~isempty(rsqrX)
  for f=1:length(rsqrX)
    rsqrPx{f} = [rsqrX(f)-0.5 rsqrX(f)+0.5 rsqrX(f)+0.5 rsqrX(f)-0.5];
  end
end

% === DISPLAY RESULTS ===
figure(Fgraph); spm_clf; 
colormap hot;

tsz = 12; % title fontsize
falpha = 0.6;

% === NORMALIZED GLOBAL MEAN
ax(1) = subplot(7,1,1);
hold on

% suspicious scans
yl = [min(ngx)-0.1 max(ngx)+0.1];
% if ~isempty(rsqrX)
%   for f = 1:length(rsqrX)
%     rsqrPy1{f} = [yl(1) yl(1) yl(2) yl(2)];
%     rsqrPhdl1(f) = patch(rsqrPx{f},rsqrPy1{f},'g');
%     set(rsqrPhdl1(f),'facecolor',yellow,'facealpha',falpha,'edgecolor','none')
%   end    
% end
% if ~isempty(mvX)
%   for f = 1:length(mvX)
%     mvPy1{f} = [yl(1) yl(1) yl(2) yl(2)];
%     mvPhdl1(f) = patch(mvPx{f},mvPy1{f},'g');
%     set(mvPhdl1(f),'facecolor',green,'facealpha',falpha,'edgecolor','none')
%   end    
% end
% if ~isempty(gmX)
%   for f=1:length(gmPx)
%     gmPy1{f} = [yl(1) yl(1) yl(2) yl(2)];
%     gmPhdl1(f) = patch(gmPx{f},gmPy1{f},'b');
%     set(gmPhdl1(f),'facecolor',red,'facealpha',falpha,'edgecolor','none')
%   end
% end
% if ~isempty(mX)
%   for f=1:length(mPx)
%     mPy1{f} = [yl(1) yl(1) yl(2) yl(2)];
%     mPhdl1(f) = patch(mPx{f},mPy1{f},'r');
%     set(mPhdl1(f),'facecolor',violet,'facealpha',falpha,'edgecolor','none')
%   end
% end

plot(1:nimg,ngx,'k','linewidth',2);
set(gca,'xlim',xl,'ylim',yl,'ygrid','on','box','on','xtick',[])
title('Normalized Global Mean','fontsize',tsz)
ylabel('SD')

% subject directory
text((nimg+2)/2,yl(2)+1.5,sprintf('Subject Directory: %s',pa),'HorizontalAlignment','center',...
	'VerticalAlignment','bottom','FontSize',14,'interpreter','none')
text((nimg+2)/2,yl(2)+3,sprintf('Spike Analyzer %s',version),'HorizontalAlignment','center',...
  'VerticalAlignment','bottom','FontSize',20,'Interpreter','none')

% === FIRST DERIVATIVE
ax(2) = subplot(7,1,2);
hold on

% threshols as gray area
xrange = [1:nimg nimg:-1:1];
yrange = [ones(1,nimg)*-cfg.globthr ones(1,nimg)*cfg.globthr];
ptch = patch(xrange,yrange,'b');
set(ptch,'facecolor',[0.8 0.8 0.8],'edgecolor',[0.8 0.8 0.8]);

% suspicious scans
yl = [min(dngx)-0.1 max(dngx)+0.1];
if ~isempty(rsqrX)
  for f = 1:length(rsqrX)
    rsqrPy2{f} = [yl(1) yl(1) yl(2) yl(2)];
    rsqrPhdl2(f) = patch(rsqrPx{f},rsqrPy2{f},'g');
    set(rsqrPhdl2(f),'facecolor',yellow,'facealpha',falpha,'edgecolor','none')
  end    
end
if ~isempty(mvX)
  for f = 1:length(mvX)
    mvPy2{f} = [yl(1) yl(1) yl(2) yl(2)];
    mvPhdl2(f) = patch(mvPx{f},mvPy2{f},'g');
    set(mvPhdl2(f),'facecolor',green,'facealpha',falpha,'edgecolor','none')
  end    
end
if ~isempty(gmX)
  for f=1:length(gmPx)
    gmPy2{f} = [yl(1) yl(1) yl(2) yl(2)];
    gmPhdl2(f) = patch(gmPx{f},gmPy2{f},'b');
    set(gmPhdl2(f),'facecolor',red,'facealpha',falpha,'edgecolor','none')
  end
end
% if ~isempty(mX)
%   for f=1:length(mPx)
%     mPy2{f} = [yl(1) yl(1) yl(2) yl(2)];
%     mPhdl2(f) = patch(mPx{f},mPy2{f},'r');
%     set(mPhdl2(f),'facecolor',violet,'facealpha',falpha,'edgecolor','none')
%   end
% end
% derivative as black line
plot(1:nimg,dngx,'k','linewidth',2);
% set axis limits
set(gca,'xlim',xl,'ylim',yl,'ygrid','on','box','on','xtick',[])
title('1st Derivative of Global Mean','fontsize',tsz)
ylabel('SD')


% === CORRELATION MOVEMENT PARAMETERS WITH GLOBAL MEAN(R-SQUARE
ax(3) = subplot(7,1,3);
hold on
colormap hot
imagesc(rsq,[0 1])
title(sprintf('Correlation between global mean and movement parameters (window: %d scans)',cfg.winsize),'fontsize',tsz)
set(gca,'xlim',xl,'ytick',[],'box','on','xtick',[])
set(gca,'position',[0.13 0.675 0.775 0.01])

% === TRANSLATIONAL MOVEMENT PARAMETERS
ax(4) = subplot(7,1,4);hold on
pos = get(gca,'position');
pos(2) = pos(2) + 0.08;
set(gca,'position',pos)
if ~isempty(mvX)
  for f = 1:length(mvX)
    mvPy5{f} = [yl(1) yl(1) yl(2) yl(2)];
    mvPhdl5(f) = patch(mvPx{f},mvPy5{f},'g');
    set(mvPhdl5(f),'facecolor',green,'edgecolor','none')
  end    
end
plot(1:nimg,rp(:,1),'r-','linewidth',2)
plot(1:nimg,rp(:,2),'g-','linewidth',2)
plot(1:nimg,rp(:,3),'b-','linewidth',2)
set(gca,'xlim',xl,'box','on','xtick',[])
title('Translational realignment parameters','fontsize',tsz)
ylabel('Millimeter')

% === ROTATIONAL MOVEMENT PARAMETERS
ax(5) = subplot(7,1,5); hold on
pos = get(gca,'position');
pos(2) = pos(2) + 0.08;
set(gca,'position',pos)
if ~isempty(mvX)
  for f = 1:length(mvX)
    mvPy6{f} = [yl(1) yl(1) yl(2) yl(2)];
    mvPhdl6(f) = patch(mvPx{f},mvPy6{f},'g');
    set(mvPhdl6(f),'facecolor',green,'edgecolor','none')
  end    
end
plot(1:nimg,rp(:,4),'r-','linewidth',2)
plot(1:nimg,rp(:,5),'g-','linewidth',2)
plot(1:nimg,rp(:,6),'b-','linewidth',2)
set(gca,'xlim',xl,'box','on','xtick',[])
title('Rotational realignment parameters','fontsize',tsz)
ylabel('Degrees')

% === MOVEMENT VELOCITY
ax(6) = subplot(7,1,6); 
hold on
pos = get(gca,'position');
pos(2) = pos(2) + 0.08;
set(gca,'position',pos)

yl = [0 max(mmtr)+0.05];
% if ~isempty(rsqrX)
%   for f = 1:length(rsqrX)
%     rsqrPy3{f} = [yl(1) yl(1) yl(2) yl(2)];
%     rsqrPhdl3(f) = patch(rsqrPx{f},rsqrPy3{f},'g');
%     set(rsqrPhdl3(f),'facecolor',yellow,'facealpha',falpha,'edgecolor','none')
%   end    
% end
% if ~isempty(mvX)
%   for f = 1:length(mvX)
%     mvPy3{f} = [yl(1) yl(1) yl(2) yl(2)];
%     mvPhdl3(f) = patch(mvPx{f},mvPy3{f},'g');
%     set(mvPhdl3(f),'facecolor',green,'facealpha',falpha,'edgecolor','none')
%   end    
% end
% if ~isempty(gmX)
%   for f=1:length(gmPx)
%     gmPy3{f} = [yl(1) yl(1) yl(2) yl(2)];
%     gmPhdl3(f) = patch(gmPx{f},gmPy3{f},'b');
%     set(gmPhdl3(f),'facecolor',red,'facealpha',falpha,'edgecolor','none')
%   end
% end
if ~isempty(mX)
  for f=1:length(mPx)
    mPy3{f} = [yl(1) yl(1) yl(2) yl(2)];
    mPhdl3(f) = patch(mPx{f},mPy3{f},'r');
    set(mPhdl3(f),'facecolor',violet,'facealpha',falpha,'edgecolor','none')
  end
end
plot(1:nimg,mmtr,'color',red,'linewidth',2)
plot(1:nimg,ones(1,nimg)*cfg.movethr,'k-.');
set(gca,'xlim',xl,'ylim',yl,'box','on','xtick',[])
set(gca,'ygrid','on');
ylabel('mm/TR')
title('Movement Velocity','fontsize',tsz)

% === GLOBAL MEAN WITH MOMVEMENT REGRESSORS REMOVED
beta  = pinv(rp)*ngx(:);
fit   = rp*beta;
resid = ngx(:) - fit;

ax(7) = subplot(7,1,7);
pos = get(gca,'position');
pos(2) = pos(2) + 0.08;
set(gca,'position',pos)
hold on

yl = [min(resid)-0.1 max(resid)+0.1];
if ~isempty(rsqrX)
  for f = 1:length(rsqrX)
    rsqrPy4{f} = [yl(1) yl(1) yl(2) yl(2)];
    rsqrPhdl4(f) = patch(rsqrPx{f},rsqrPy4{f},'g');
    set(rsqrPhdl4(f),'facecolor',yellow,'facealpha',falpha,'edgecolor','none')
  end    
end
if ~isempty(mvX)
  for f = 1:length(mvX)
    mvPy4{f} = [yl(1) yl(1) yl(2) yl(2)];
    mvPhdl4(f) = patch(mvPx{f},mvPy4{f},'g');
    set(mvPhdl4(f),'facecolor',green,'facealpha',falpha,'edgecolor','none')
  end    
end
if ~isempty(gmX)
  for f=1:length(gmPx)
    gmPy4{f} = [yl(1) yl(1) yl(2) yl(2)];
    gmPhdl4(f) = patch(gmPx{f},gmPy4{f},'b');
    set(gmPhdl4(f),'facecolor',red,'facealpha',falpha,'edgecolor','none')
  end
end
if ~isempty(mX)
  for f=1:length(mPx)
    mPy4{f} = [yl(1) yl(1) yl(2) yl(2)];
    mPhdl4(f) = patch(mPx{f},mPy4{f},'r');
    set(mPhdl4(f),'facecolor',violet,'facealpha',falpha,'edgecolor','none')
  end
end
plot(1:nimg,resid,'color','black','linewidth',2);
set(gca,'xlim',xl,'ylim',yl,'box','on')
set(gca,'ygrid','on');
xlabel('scan')
title('Normalized global mean (without effects of movement parameter)','fontsize',tsz)
ylabel('SD')

newy = y(1) - 1.1*(yl(2)-yl(1));
newx = xl(1);
dy = 0.95;

str1 = sprintf('No. of candidate scans suggested for removal through separate modeling:');
str1a = sprintf('global mean/mvmt correlation (%1.1f): %d/%d (%1.2f%%)',...
  cfg.rsqrthr,length(rsqrX),nimg,(length(rsqrX)/nimg)*100);
str2 = sprintf('1st derivative of global mean (%1.1f SD): %d/%d (%1.2f%%)',...
  cfg.globthr,length(gmX),nimg,(length(gmX)/nimg)*100);
str2a = sprintf('movement parameters (%1.1f * voxsize): %d/%d (%1.2f%%)',...
  cfg.globthr,length(mvX),nimg,(length(mvX)/nimg)*100);
str3 = sprintf('movement velocity (%1.2f mm/TR): %d/%d (%1.2f%%)',...
  cfg.movethr,length(mX),nimg,length(mX)/nimg);
str4 = sprintf('TOTAL: %d/%d (%1.2f%%)',...
  length(unique([gmX(:); mX(:); rsqrX(:); mvX(:)])),nimg,(length(unique([gmX(:); mX(:); rsqrX(:); mvX(:)]))/nimg)*100);

text(newx,newy,str1,'horizontalalignment','left','verticalalignment','top',...
  'color','k','fontsize',14)
text(newx,newy-dy,str1a,'horizontalalignment','left','verticalalignment','top',...
  'color',yellow,'fontsize',14)
text(newx,newy-2*dy,str2,'horizontalalignment','left','verticalalignment','top',...
  'color',red,'fontsize',14)
text(newx,newy-3*dy,str2a,'horizontalalignment','left','verticalalignment','top',...
  'color',green,'fontsize',14)
text(newx,newy-4*dy,str3,'horizontalalignment','left','verticalalignment','top',...
  'color',violet,'fontsize',14)
text(newx,newy-5*dy,str4,'horizontalalignment','left','verticalalignment','top',...
  'color','k','fontsize',14)

linkaxes(ax(:)','x'); %for zooming on all subplot only in x direction

if cfg.pflag
	print(sprintf('-f%f',Fgraph),'-append','-dpsc2',psfile)
end

if ~isempty(rsqrX)
  rsqrNuis = zeros(nimg,length(rsqrX));
  for f=1:length(rsqrX)
    rsqrNuis(rsqrX(f),f) = 1;
  end
else
  rsqrNuis = [];
end

if ~isempty(gmX)
  gmNuis = zeros(nimg,length(gmX));
  for f=1:length(gmX)
    gmNuis(gmX(f),f) = 1;
  end
else
  gmNuis = [];
end

if ~isempty(mX)
  mNuis = zeros(nimg,length(mX));
  for f=1:length(mX)
    mNuis(mX(f),f) = 1;
  end
else
  mNuis = [];
end

if ~isempty(mvX)
  realignNuis = zeros(nimg,length(mvX));
  for f=1:length(mvX)
    realignNuis(mvX(f),f) = 1;
  end
else
  realignNuis = [];
end


allNuis = unique([gmNuis mNuis rsqrNuis realignNuis]','rows')';

nuis.gm = gmNuis;
nuis.m  = mNuis;
nuis.rsqr = rsqrNuis;
nuis.realign = realignNuis;
nuis.all = allNuis;

return

function segment_epi(img,pth)
% Use spm_segment and imcalc to compute a brain mask

matlabbatch{1}.spm.tools.oldseg.data = {deblank(img)};
matlabbatch{1}.spm.tools.oldseg.output.GM = [0 0 1];
matlabbatch{1}.spm.tools.oldseg.output.WM = [0 0 1];
matlabbatch{1}.spm.tools.oldseg.output.CSF = [0 0 0];
matlabbatch{1}.spm.tools.oldseg.output.biascor = 0;
matlabbatch{1}.spm.tools.oldseg.output.cleanup = 0;
matlabbatch{1}.spm.tools.oldseg.opts.tpm = {
  fullfile(spm('dir'),'toolbox','OldSeg','grey.nii')
  fullfile(spm('dir'),'toolbox','OldSeg','white.nii')
  fullfile(spm('dir'),'toolbox','OldSeg','csf.nii')};
matlabbatch{1}.spm.tools.oldseg.opts.ngaus = [2 2 2 4];
matlabbatch{1}.spm.tools.oldseg.opts.regtype = 'mni';
matlabbatch{1}.spm.tools.oldseg.opts.warpreg = 1;
matlabbatch{1}.spm.tools.oldseg.opts.warpco = 25;
matlabbatch{1}.spm.tools.oldseg.opts.biasreg = 0.0001;
matlabbatch{1}.spm.tools.oldseg.opts.biasfwhm = 60;
matlabbatch{1}.spm.tools.oldseg.opts.samp = 3;
matlabbatch{1}.spm.tools.oldseg.opts.msk = {''};
matlabbatch{2}.spm.util.imcalc.input(1) = cfg_dep('Old Segment: c1 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','c1', '()',{':'}));
matlabbatch{2}.spm.util.imcalc.input(2) = cfg_dep('Old Segment: c2 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','c2', '()',{':'}));
matlabbatch{2}.spm.util.imcalc.output = 'epi_brainmask';
matlabbatch{2}.spm.util.imcalc.outdir = {pth};
matlabbatch{2}.spm.util.imcalc.expression = '(i1+i2)>0.2';
matlabbatch{2}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{2}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{2}.spm.util.imcalc.options.mask = 0;
matlabbatch{2}.spm.util.imcalc.options.interp = 1;
matlabbatch{2}.spm.util.imcalc.options.dtype = 4;
matlabbatch{3}.cfg_basicio.file_dir.file_ops.file_move.files(1) = cfg_dep('Old Segment: c1 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','c1', '()',{':'}));
matlabbatch{3}.cfg_basicio.file_dir.file_ops.file_move.files(2) = cfg_dep('Old Segment: c2 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','c2', '()',{':'}));
matlabbatch{3}.cfg_basicio.file_dir.file_ops.file_move.action.delete = false;

spm_jobman('run',matlabbatch);

return

