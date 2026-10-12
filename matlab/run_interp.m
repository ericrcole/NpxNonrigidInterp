function run_interp(segdir, repoDir)
% RUN_INTERP  Apply a DREDge motion estimate to an AP segment by kriging interpolation.
%
%   run_interp('<segment folder>')
%   run_interp('<segment folder>', '<InterpolationAfterDREDge folder>')
%
% The segment folder is what make_segment.py wrote: continuous.dat, chanMap.mat,
% nonrigid_displacement.csv, dredge_window_centers.csv, ref_displacement.csv and
% interp_params.json. All parameters come from that json, so this file never needs editing.
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
RealignRecordingUsingDredgeAndKilosort2point5(segdir, 'Dredge', dredgeParams);
disp(['INTERPOLATION COMPLETE -> ' fullfile(segdir, p.out_name)]);
end
