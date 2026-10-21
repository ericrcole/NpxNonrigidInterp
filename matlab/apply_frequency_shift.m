function y = apply_frequency_shift(x, shiftSamples)
% APPLY_FREQUENCY_SHIFT  Sub-sample time shift of each channel, in the frequency domain.
%
%   y = apply_frequency_shift(x, shiftSamples)
%       x            [nSamples x nChannels] double
%       shiftSamples [1 x nChannels] shift per channel, in units of 1/fs
%
% Port of spikeinterface.preprocessing.phase_shift.apply_frequency_shift (which follows
% IBL's ibllib.dsp.fourier and SpikeGLX), kept deliberately line-for-line equivalent:
%
%   frequency_grid = 2*pi*rfftfreq(N)              (numpy rfftfreq: k/N, k = 0..floor(N/2))
%   rotations      = exp(-1j * frequency_grid * shift_samples)
%   y              = irfft(rfft(x) * rotations, N)
%
% MATLAB has no rfft, so the full spectrum is built with the conjugate-symmetric frequency
% convention f(k) = k/N for k <= N/2 and (k-N)/N above it. For real input that is identical
% to numpy's rfft/irfft round trip: the negative half carries the conjugate rotation, and the
% Nyquist bin's imaginary part cancels in the inverse transform exactly as irfft discards it.
% Verified against SpikeInterface on real data by scripts/check_phase_shift.py.
%
% This corrects the Neuropixels ADC multiplexing delay: channels are digitised at staggered
% times within one sample period (inter_sample_shift, 0 to ~0.85 samples = 28 us on NP1.0),
% and any operation that MIXES channels - kriging interpolation here, as well as CAR or
% whitening - must come after this correction, because mixing time-skewed copies distorts
% the waveform in a way nothing downstream can undo.
N = size(x, 1);
if N == 0 || all(shiftSamples == 0)
    y = x;
    return;
end
k = (0:N-1)';
k(k > floor(N/2)) = k(k > floor(N/2)) - N;      % conjugate-symmetric frequency index
freq = 2 * pi * k / N;                          % 2*pi*rfftfreq, mirrored
rot = exp(-1i * (freq * shiftSamples(:)'));     % [N x nChannels]
y = real(ifft(fft(x, N, 1) .* rot, N, 1));
end
