function rez = applyDredge(rez)
% -------------------------------------------------------------------------------------------
% MODIFIED COPY. Part of the NPIX_MSIT DREDge + kriging pipeline; see matlab/README.md for what
% changed and ../CREDITS.md for provenance.
% Original: InterpolationAfterDREDge by Richard Hardstone and William Munoz (CC0 1.0),
%   https://github.com/williamunoz/InterpolationAfterDREDge
% -------------------------------------------------------------------------------------------
%% Adapted from Kilosort2.5 function datashift2.m
%% by Richard Hardstone (2022_11_30)

%%
outDir = rez.ops.dredge.outDir;
if ~exist(outDir,'dir')
    mkdir(outDir);
end

%%
if  getOr(rez.ops, 'nblocks', 1)==0
    rez.iorig = 1:rez.temp.Nbatch;
    return;
end

ops = rez.ops;

% The min and max of the y and x ranges of the channels
ymin = min(rez.yc);
ymax = max(rez.yc);
xmin = min(rez.xc);
xmax = max(rez.xc);

% Determine the average vertical spacing between channels.
% Usually all the vertical spacings are the same, i.e. on Neuropixels probes.
dmin = median(diff(unique(rez.yc)));
fprintf('pitch is %d um\n', dmin)
rez.ops.yup = ymin:dmin/2:ymax; % centers of the upsampled y positions

% Determine the template spacings along the x dimension
xrange = xmax - xmin;
npt = floor(xrange/16); % this would come out as 16um for Neuropixels probes, which aligns with the geometry.
rez.ops.xup = linspace(xmin, xmax, npt+1); % centers of the upsampled x positions

spkTh = 8; % same as the usual "template amplitude", but for the generic templates

% Extract all the spikes across the recording that are captured by the
% generic templates. Very few real spikes are missed in this way.
% ERC: the drift map runs standalone_detector over the whole recording, which is
% expensive and unnecessary when validating a short batch range.
skipDriftMap = getOr(rez.ops.dredge, 'skipDriftMap', false);
if ~skipDriftMap
[st3, rez] = standalone_detector(rez, spkTh);

%%
close all
figure(194);
set(gcf, 'Color', 'w')
% raster plot of all spikes at their original depths
st_shift = st3(:,2); %+ imin(batch_id)' * dd;
for j = spkTh:100
    % for each amplitude bin, plot all the spikes of that size in the
    % same shade of gray
    ix = st3(:, 3)==j; % the amplitudes are rounded to integers
    plot(st3(ix, 1)/ops.fs, st_shift(ix), '.', 'color', [1 1 1] * max(0, 1-j/40)) % the marker color here has been carefully tuned
    hold on
end
axis tight

figDir = fullfile(outDir,'DriftMapFigures');
if ~exist(figDir,'dir')
    mkdir(figDir);
end
xlabel('time (sec)')
ylabel('spike position (um)')
title('Drift map')
print('-dpng',fullfile(figDir, 'DriftMap_BeforeRealignment.png'))
saveas(194,fullfile(figDir, 'DriftMap_BeforeRealignment.fig'))
drawnow;
end


%%
hread = readmatrix(rez.ops.dredge.Fname)';
hread = resample(hread, rez.ops.dredge.dredgeSamples, rez.ops.dredge.BatchSamplesNT);

% ERC: optional zero-phase temporal low-pass on the displacement estimate.
% Measured on this recording, 98.2% of the estimate's variance is below 0.05 Hz, but
% ~20 um RMS sits above 0.1 Hz and ~8 um RMS above 2 Hz. Tissue does not move at 5 Hz,
% so that is estimator noise and the pipeline was applying it to the data at the batch
% rate. filtfilt so the drift is not time-shifted relative to the recording.
lpHz = getOr(rez.ops.dredge, 'motionLowpassHz', []);
if ~isempty(lpHz)
    fsBatch = rez.ops.fs / rez.ops.dredge.BatchSamplesNT;   % 234.375 Hz at NT = 128
    [blp, alp] = butter(3, lpHz/(fsBatch/2), 'low');
    pre = median(std(hread, 0, 1));
    hread = filtfilt(blp, alp, hread);
    fprintf('motion low-passed at %.2f Hz: median SD %.1f -> %.1f um\n', ...
            lpHz, pre, median(std(hread, 0, 1)));
end

% ERC: optional smoothing of the displacement ACROSS DEPTH windows.
% The nonrigid fit at win_step_um=200 produced differential motion between adjacent
% windows of median 52 um over 200 um (strain 0.26, >1.0 for 3.6% of samples), which
% no real tissue does. Smoothing across windows, or collapsing to rigid, constrains it.
% Prefer re-running DREDge with larger windows; this is the salvage path for an
% existing estimate.
dsm = getOr(rez.ops.dredge, 'depthSmoothWindows', 0);
if dsm > 0 && size(hread,2) > 1
    k = ones(1, min(dsm, size(hread,2))) / min(dsm, size(hread,2));
    hread = conv2(hread, k, 'same') ./ conv2(ones(size(hread)), k, 'same');
    fprintf('displacement smoothed across %d depth windows\n', dsm);
end
if getOr(rez.ops.dredge, 'forceRigid', false) && size(hread,2) > 1
    hread = repmat(mean(hread, 2), 1, size(hread,2));
    fprintf('displacement collapsed to rigid (mean across depth)\n');
end

% ERC: the displacement reference. Default keeps the original behaviour (each run
% centres on its own median), but a run over a SUBSET of batches would then use a
% different reference than a full pass. Passing dredgeParams.refDisplacement
% (1 x nWindows, from the full recording) makes a segment run bit-equivalent to
% the corresponding batches of a full run.
refDisp = getOr(rez.ops.dredge, 'refDisplacement', []);
if isempty(refDisp)
    refDisp = median(hread);
else
    refDisp = refDisp(:)';
    if numel(refDisp) ~= size(hread,2)
        error('refDisplacement has %d entries but motion has %d windows', numel(refDisp), size(hread,2));
    end
    fprintf('using supplied displacement reference (not this segment''s median)\n');
end
hread = refDisp - hread;
rez.ops.NT  = rez.ops.dredge.BatchSamplesNT;

% ERC: length() returns the LARGEST dimension, which is the time axis only by
% luck. Use size(...,1), and clamp to the number of batches the file can actually
% supply -- the resampled motion trace overruns the recording by a few tens of
% samples, and the old code wrote that overrun past EOF, growing the output file.
% ERC: take the length from the SOURCE. With streaming output, an existing output file may
% be stale or partial from an earlier run, and must not set how much gets processed.
d = dir(rez.ops.fbinary);
nFramesFile = floor(d.bytes / (2*rez.ops.NchanTOT));
Nbatches = min(size(hread,1), ceil(nFramesFile / rez.ops.dredge.BatchSamplesNT));
fprintf('motion trace supplies %d batches; file holds %d frames -> using %d batches\n', ...
        size(hread,1), nFramesFile, Nbatches);
rez.temp.Nbatch = Nbatches;
sig = rez.ops.sig;

% Load spatial window centers from CSV (one value per column of hread).
% This is the depth (um) of each DREDge non-rigid window, saved alongside
% dredge.csv from Python as: np.savetxt('dredge_window_centers.csv', motion.spatial_bins_um, delimiter=',')
ysamp = readmatrix(rez.ops.dredge.windowCentersFname);
ysamp = ysamp(:); % ensure column vector

% Validate: number of window centers must match columns of hread
nWindows = size(hread, 2);
if length(ysamp) ~= nWindows
    error('applyDredge: dredge.csv has %d columns but dredge_window_centers.csv has %d entries. They must match.', ...
          nWindows, length(ysamp));
end

% If only one window (rigid), pass [] so shift_batch_on_disk2_dredge skips
% interpolation and applies a single scalar shift uniformly.
if nWindows == 1
    ysamp = [];
end

% temporal shift-interpolation window = full batch -> shift ramps continuously
% batch-to-batch instead of holding constant (removes the 234 Hz staircase)
rez.ops.ntbuff = getOr(rez.ops.dredge, 'blendSamples', rez.ops.dredge.BatchSamplesNT);
ops = rez.ops;                                   % re-sync local copy so dprev is sized right
sig = rez.ops.sig;
dprev = gpuArray.zeros(ops.ntbuff, ops.Nchan, 'single');

copyfile('chanMap.mat', fullfile(outDir,'chanMap.mat'));

% ERC: Kxx depends only on the electrode geometry, which never changes, so factor it
% once here instead of once per batch (282,143 times over the full recording).
xp0 = cat(2, rez.xc(:), rez.yc(:));
Kxx = kernel2D(xp0, xp0, sig);
% nugget raised from .01 to .05 given poor conditioning of Kxx - ERC
rez.ops.dredge.iKxx = inv(Kxx + .05 * eye(size(Kxx,1)));

br = getOr(rez.ops.dredge, 'batchRange', [1 Nbatches]);
b0 = max(1, br(1)); b1 = min(Nbatches, br(2));

% ERC: stream source -> dest. The old code copyfile'd the entire recording (27.8 GB)
% and then rewrote nearly all of it in place. Reading from a pristine source and
% writing through a separate handle removes the copy, halves the write I/O, and makes
% read/write aliasing impossible. Handles stay open instead of fclose all/fopen per
% batch.
nRaw  = rez.ops.NchanTOT;
fidIn = fopen(rez.ops.fbinary, 'r');
if fidIn < 0, error('cannot open source %s', rez.ops.fbinary); end
if b0 == 1 || ~exist(rez.ops.dredge.outDat, 'file')
    fidOut = fopen(rez.ops.dredge.outDat, 'w');
else
    fidOut = fopen(rez.ops.dredge.outDat, 'r+');   % partial range into an existing file
end
if fidOut < 0, fclose(fidIn); error('cannot open output %s', rez.ops.dredge.outDat); end
closer = onCleanup(@() [fclose(fidIn), fclose(fidOut)]);

% ERC: per-channel DC offsets, estimated once from 20 short probes spread across the file,
% so kriging does not mix offsets between channels (see shift_batch_on_disk2_dredge).
if getOr(rez.ops.dredge, 'removeOffset', true)
    nProbe = 20; nS = 3000;
    probeStarts = round(linspace(0, max(0, nFramesFile - nS), nProbe));
    offAcc = zeros(nProbe, rez.ops.Nchan);
    for ip = 1:nProbe
        fseek(fidIn, 2*nRaw*probeStarts(ip), 'bof');
        blk = fread(fidIn, [nRaw nS], '*int16');
        offAcc(ip, :) = median(double(blk(1:rez.ops.Nchan, :)), 2)';
    end
    rez.ops.dredge.chanOffset = median(offAcc, 1);
    disp(sprintf('per-channel DC offset removed before kriging: median %.0f, SD across channels %.0f', ...
        median(rez.ops.dredge.chanOffset), std(rez.ops.dredge.chanOffset)));
end

clear shift_batch_on_disk2_dredge   % drop any cached kriging matrix from a prior run

fprintf('realigning batches %d..%d of %d (streaming, no pre-copy)\n', b0, b1, Nbatches);
tic;
for ibatch = b0:b1
    if mod(ibatch - b0, 5000) == 0 || ibatch == b1
        el = toc; frac = (ibatch-b0+1)/(b1-b0+1);
        fprintf('  %.2f%% (%d/%d, %.0f s elapsed, ~%.0f s remaining)\n', ...
                100*frac, ibatch-b0+1, b1-b0+1, el, el/max(frac,1e-9)-el);
    end
    dprev = shift_batch_on_disk2_dredge(rez, ibatch, hread(ibatch, :), ysamp, sig, dprev, fidIn, fidOut);
end
fprintf('realignment done in %.0f s\n', toc);

% ERC: copy any frames past the last batch verbatim, so the output is exactly the same
% length as the input rather than truncated or over-extended.
tailStart = Nbatches * rez.ops.NT;
if tailStart < nFramesFile
    fprintf('copying %d trailing frames verbatim\n', nFramesFile - tailStart);
    fseek(fidIn,  2*nRaw*tailStart, 'bof');
    fseek(fidOut, 2*nRaw*tailStart, 'bof');
    while true
        blk = fread(fidIn, [nRaw 100000], '*int16');
        if isempty(blk), break; end
        fwrite(fidOut, blk, 'int16');
    end
end
