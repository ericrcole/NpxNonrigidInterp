function run_interp(segdir, repoDir)
% RUN_INTERP  Apply a DREDge motion estimate to an AP segment by kriging interpolation.
%
%   run_interp('<segment folder>')
%   run_interp('<segment folder>', '<InterpolationAfterDREDge folder>')
%
% The segment folder is what make_segment.py wrote: continuous.dat, chanMap.mat,
% nonrigid_displacement.csv, dredge_window_centers.csv, ref_displacement.csv,
% inter_sample_shift.csv and interp_params.json. All parameters come from that json, so this
% file never needs editing. The ADC phase shift is applied before kriging by default.
% From a terminal:
%
%   matlab -batch "addpath('<this matlab folder>'); run_interp('<segment folder>')"
%
% Output: <segment folder>/<outName>/<outName>.dat, the same size as the input.
if nargin < 2 || isempty(repoDir)
    repoDir = getenv('INTERP_REPO');
end
if isempty(repoDir)
    error(['Pass the InterpolationAfterDREDge folder as the 2nd argument, or set the ' ...
           'INTERP_REPO environment variable.']);
end
if exist(fullfile(segdir, 'interp_params.json'), 'file') ~= 2
    error('No interp_params.json in %s - run make_segment.py first.', segdir);
end
restoredefaultpath;
addpath(genpath(fullfile(repoDir, 'Scripts', 'AP')));
addpath(genpath(fullfile(repoDir, 'Scripts', 'ExternalPackages')));
% This folder goes on LAST so it comes FIRST on the path: the modified applyDredge,
% shift_batch_on_disk2_dredge, RealignRecordingUsingDredgeAndKilosort2point5 and
% Kilosort25Configuration_Dredge here shadow the repo's originals, and the phase-shift
% helpers are found at all. Without this, restoredefaultpath above would drop whatever the
% caller added and the stock repo versions would run instead - silently, since they take the
% same arguments. The repo still supplies everything else (kernel2D, getOr, loadChanMap, ...).
addpath(fileparts(mfilename('fullpath')));

p = jsondecode(fileread(fullfile(segdir, 'interp_params.json')));
dredgeParams.Fname              = 'nonrigid_displacement.csv';
dredgeParams.windowCentersFname = 'dredge_window_centers.csv';
dredgeParams.fs                 = p.motion_fs;      % sample rate of the displacement trace
dredgeParams.BatchSamplesNT     = p.batch_samples;  % samples per interpolation batch
dredgeParams.NchanRaw           = p.nchan_raw;      % channels in the FILE, sync included
dredgeParams.skipDriftMap       = true;             % motion comes from DREDge, not Kilosort
dredgeParams.removeOffset       = true;             % no DC-offset mixing between channels
dredgeParams.mTolUm             = p.m_tol_um;       % kriging operator cache tolerance
dredgeParams.refDisplacement    = readmatrix(fullfile(segdir, 'ref_displacement.csv'));
dredgeParams.outName            = p.out_name;

% Neuropixels ADC phase shift, applied to the AP data BEFORE kriging. ON unless the segment
% says otherwise, because kriging mixes channels that were digitised at different times
% within a sample period and nothing downstream can undo that. make_segment.py writes
% inter_sample_shift.csv (one value per neural channel, in units of 1/fs) from the
% recording's own metadata; the correction itself matches spikeinterface.phase_shift.
dredgeParams.phaseShift = true;
if isfield(p, 'phase_shift')
    dredgeParams.phaseShift = logical(p.phase_shift);
end
if isfield(p, 'phase_shift_margin_ms')
    dredgeParams.psMarginMs = p.phase_shift_margin_ms;
end
if dredgeParams.phaseShift
    issFile = fullfile(segdir, 'inter_sample_shift.csv');
    if exist(issFile, 'file') ~= 2
        error(['No inter_sample_shift.csv in %s. Re-run make_segment.py to write it, or set ' ...
               '"phase_shift": false in interp_params.json to interpolate without the ADC ' ...
               'correction (not recommended).'], segdir);
    end
    dredgeParams.interSampleShift = readmatrix(issFile);
end
RealignRecordingUsingDredgeAndKilosort2point5(segdir, 'Dredge', dredgeParams);
disp(['INTERPOLATION COMPLETE -> ' fullfile(segdir, p.out_name)]);
end
