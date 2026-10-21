function dump_phase_shift(segdir, f0, n, outfile)
% DUMP_PHASE_SHIFT  Write this pipeline's phase-shifted samples for one window, for testing.
%
%   dump_phase_shift('<segment folder>', f0, n, '<out.bin>')
%
% Reads frames [f0, f0+n) through phase_shift_reader in the same 128-sample batch pattern the
% interpolation uses, so the block-cache boundaries are exercised, and writes the result as
% float32 [n x nChan]. scripts/check_phase_shift.py compares it with spikeinterface.phase_shift.
p = jsondecode(fileread(fullfile(segdir, 'interp_params.json')));
iss = readmatrix(fullfile(segdir, 'inter_sample_shift.csv'));
d = dir(fullfile(segdir, 'continuous.dat'));
nRaw = p.nchan_raw;
nFramesFile = floor(d.bytes / (2 * nRaw));
nChan = numel(iss);
fs = 30000;
margin = round(p.phase_shift_margin_ms * fs / 1000);
blockFrames = 60000;
fid = fopen(fullfile(segdir, 'continuous.dat'), 'r');
cl = onCleanup(@() fclose(fid));
phase_shift_reader();                                    % reset the cache
X = zeros(n, nChan, 'single');
NT = p.batch_samples;
for s = 0:NT:(n-1)
    len = min(NT, n - s);
    X(s + (1:len), :) = single(phase_shift_reader(fid, f0 + s, len, nRaw, nChan, ...
                                                  iss, margin, blockFrames, nFramesFile));
end
fo = fopen(outfile, 'w');
fwrite(fo, X', 'single');
fclose(fo);
fprintf('wrote %s: %d frames x %d channels\n', outfile, n, nChan);
end
