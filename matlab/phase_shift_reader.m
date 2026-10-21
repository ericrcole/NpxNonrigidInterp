function X = phase_shift_reader(fidIn, f0, n, nRaw, nChan, shiftSamples, margin, blockFrames, nFramesFile)
% PHASE_SHIFT_READER  Read frames [f0, f0+n) with the ADC phase shift already applied.
%
%   X = phase_shift_reader(fidIn, f0, n, nRaw, nChan, shiftSamples, margin, blockFrames, nFramesFile)
%       returns [nRead x nChan] double of the NEURAL channels, phase-shifted.
%       Call phase_shift_reader() with no arguments (or `clear phase_shift_reader`) to reset.
%
% Why a cache: the interpolation works in 128-sample batches, but a frequency-domain shift
% needs a margin around the data it transforms (SpikeInterface uses 40 ms, which it documents
% as giving a very small chunking error). Transforming 128 samples with a 1200-sample margin
% on each side would do ~20x redundant work per batch, so this reads one large block at a
% time, shifts it once, and serves the sequential batch requests out of it.
%
% The block itself is prepared exactly as SpikeInterface does per chunk
% (spikeinterface.core.get_chunk_with_margin with add_zeros=True, window_on_margin=True):
% read the block plus `margin` samples on each side, zero-pad where that runs off the file,
% taper both margins with a raised cosine, apply the shift, then keep the centre.
persistent BUF B0 BLEN KEY
if nargin == 0
    BUF = []; B0 = []; BLEN = []; KEY = [];
    return;
end
key = [double(fidIn), nRaw, nChan, margin, blockFrames, numel(shiftSamples), sum(shiftSamples)];
if isempty(KEY) || ~isequal(KEY, key)
    BUF = []; B0 = []; BLEN = []; KEY = key;
end

n = min(n, max(0, nFramesFile - f0));
if n <= 0
    X = zeros(0, nChan);
    return;
end

if isempty(BUF) || f0 < B0 || (f0 + n) > (B0 + BLEN)
    b0 = f0;
    blen = min(max(blockFrames, n), nFramesFile - b0);
    r0 = max(0, b0 - margin);                       % first frame actually in the file
    r1 = min(nFramesFile, b0 + blen + margin);      % one past the last
    fseek(fidIn, 2 * nRaw * r0, 'bof');
    raw = fread(fidIn, [nRaw, r1 - r0], '*int16');
    buf = zeros(blen + 2 * margin, nChan);          % zero padding outside the file
    i0 = margin - (b0 - r0);                        % where the real data starts in buf
    buf(i0 + (1:size(raw, 2)), :) = double(raw(1:nChan, :))';
    taper = (1 - cos((0:margin-1)' / margin * pi)) / 2;
    buf(1:margin, :)         = buf(1:margin, :)        .* taper;
    buf(end-margin+1:end, :) = buf(end-margin+1:end, :) .* flipud(taper);
    buf = apply_frequency_shift(buf, shiftSamples);
    BUF = buf(margin + (1:blen), :);
    B0 = b0; BLEN = blen;
end

X = BUF(f0 - B0 + (1:n), :);
end
