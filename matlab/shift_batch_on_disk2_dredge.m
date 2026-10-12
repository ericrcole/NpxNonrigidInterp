function dprev = shift_batch_on_disk2_dredge(rez, ibatch, shifts, ysamp, sig, dprev, fidIn, fidOut)
% -------------------------------------------------------------------------------------------
% MODIFIED COPY. Part of the NPIX_MSIT DREDge + kriging pipeline; see matlab/README.md for what
% changed and ../CREDITS.md for provenance.
% Original: InterpolationAfterDREDge by Richard Hardstone and William Munoz (CC0 1.0),
%   https://github.com/williamunoz/InterpolationAfterDREDge
% -------------------------------------------------------------------------------------------

%% Adapted from Kilosort2.5 function shift_batch_on_disk2.m
%% by Richard Hardstone (2022_11_30)
%
% Reads one batch from the RAW source file and writes the realigned batch to a
% SEPARATE output file. Nothing is copied up front and the handles stay open.
%
% CHANNEL COUNTS (ERC): the Kilosort original works on ops.fproc, which is exactly
% ops.Nchan wide, so it used ops.Nchan for all file arithmetic. This version works on
% the RAW file, which is ops.NchanTOT wide (384 neural + 1 sync on Open Ephys NP1.0).
% File geometry must use NchanTOT; only the kriging uses Nchan. Getting this wrong
% sheared each batch diagonally across the probe and mixed the sync line into the
% neural channels.

persistent Mcache shiftsCache

ops    = rez.ops;
Nbatch = rez.temp.Nbatch;
NT     = ops.NT;
ntb    = ops.ntbuff;
nRaw   = ops.NchanTOT;          % channels per frame in the file (incl. sync/aux)
ichan  = 1:ops.Nchan;           % channels that get interpolated

batchstart = 0:NT:NT*Nbatch;
offset = 2 * nRaw * batchstart(ibatch);

% ---- upsample the shift to every channel ----------------------------------
if length(ysamp) > 1
    % pchip + edge clamp rather than makima 'extrap', which explodes cubically
    % at the edges of the fit - ERC
    ys = ysamp(:);  sh = shifts(:);
    [ys, iu] = unique(ys, 'stable');  sh = sh(iu);
    [ys, is] = sort(ys);              sh = sh(is);
    yq = min(max(rez.yc(:), ys(1)), ys(end));
    shifts = interp1(ys, sh, yq, 'pchip');
    shifts = reshape(shifts, size(rez.yc));
end

if isscalar(shifts)
    shifts = repmat(shifts, size(rez.yc));   % rigid estimate: same shift on every channel
end

% ---- taper the shift to zero above the brain surface ----------------------
% ERC: above maxInterpDepth the probe is in saline. Holding the clamped edge-window
% shift there moved those channels and generated spurious threshold crossings; a HARD
% cut instead left a step between corrected and uncorrected channels, which showed up
% as a spurious band right at the boundary. Tapering the shift smoothly to zero keeps
% the saline channels effectively unmoved with no discontinuity anywhere.
maxD = getOr(ops.dredge, 'maxInterpDepth', inf);
if isfinite(maxD)
    taper = getOr(ops.dredge, 'taperUm', 200);
    y = rez.yc(:);
    w = ones(size(y));
    ramp = (y > maxD - taper) & (y <= maxD);
    w(ramp) = 0.5 * (1 + cos(pi * (y(ramp) - (maxD - taper)) / taper));
    w(y > maxD) = 0;
    shifts = shifts .* reshape(w, size(shifts));
end

% ---- kriging prediction matrix (cached across batches) --------------------
% Kxx depends only on the electrode geometry, so its inverse is precomputed once in
% applyDredge. Kyx moves with the shift, but at NT=128 the shift changes by ~0.1 um
% between batches, far below the kernel scale, so M is reused until it moves by mTol.
mTol = getOr(ops.dredge, 'mTolUm', 0.5);
if isempty(Mcache) || numel(shiftsCache) ~= numel(shifts) || ...
        max(abs(shifts(:) - shiftsCache(:))) > mTol
    xp = cat(2, rez.xc(:), rez.yc(:));
    yp = xp;  yp(:, 2) = yp(:, 2) - shifts(:);
    Mcache = kernel2D(yp, xp, sig) * ops.dredge.iKxx;
    shiftsCache = shifts;
end

% ---- read from SOURCE ------------------------------------------------------
fseek(fidIn, offset, 'bof');
dat   = fread(fidIn, [nRaw NT+ntb], '*int16')';     % [nread x nRaw]
nread = size(dat, 1);
if nread == 0, return; end

% ERC: remove each channel's DC offset before kriging and restore the OUTPUT channel's own
% offset afterwards. Raw AP channels carry different DC offsets (msit004: SD 194 counts
% across channels). Kriging with time-varying weights mixes those offsets, producing a
% slow, channel-specific signal that follows the shift trajectory - measured on msit004
% as a 3x rise in <150 Hz power, predicted from the offsets to within 1% (r = 0.98).
% It sits below Kilosort's high-pass, but it is an artifact and costs nothing to remove.
off = getOr(ops.dredge, 'chanOffset', []);
if isempty(off)
    dati = gpuArray(single(dat(:, ichan))) * gpuArray(Mcache)';
else
    offg = gpuArray(single(off(:)'));
    dati = (gpuArray(single(dat(:, ichan))) - offg) * gpuArray(Mcache)' + offg;
end

% raised cosine (was linear) -> zero slope at both ends, so the effective shift is
% C1-continuous across batch joins, not merely continuous - ERC
w_edge = 0.5 * (1 - cos(pi * linspace(0, 1, ntb)'));
nb = min(ntb, size(dati,1));
dati(1:nb, :) = w_edge(1:nb) .* dati(1:nb, :) + (1 - w_edge(1:nb)) .* dprev(1:nb, :);

if size(dati,1) >= NT+ntb
    dprev = dati(NT+(1:ntb), :);
else
    dprev(:) = 0;
    if size(dati,1) > NT
        tmp = dati(NT+1:end, :);
        dprev(1:size(tmp,1), :) = tmp;
    end
end

% ---- write full frames to DEST, sync/aux channel passed through untouched --
% The sync line is the trial-alignment reference: shifting it would desynchronise the
% behavioural data. Leaving it verbatim also gives a free integrity check.
nwrite = min(NT, nread);
out = dat(1:nwrite, :);
out(:, ichan) = gather(int16(dati(1:nwrite, :)));

fseek(fidOut, offset, 'bof');
fwrite(fidOut, out', 'int16');
