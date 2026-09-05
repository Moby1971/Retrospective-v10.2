% How far a spectral peak stands above the floor of its own band
%
% Author : Gustav Strijkers
% Date   : 2026-09-04

function p = bandProminence(x, sf, bandBpm)

    % Purpose:
    %   Score a navigator by how clearly it shows a rhythm inside a given rate
    %   band, so that several candidate coil combinations can be compared on the
    %   thing that matters rather than on their variance.
    %
    % How it works:
    %   The score is the height of the strongest peak in the band divided by the
    %   median level of the band. A ratio rather than an absolute power, because
    %   the candidates being compared have wildly different amplitudes -- a coil
    %   close to the surface can carry ten times the signal of one further away
    %   and none of the motion -- and what is being asked is not how loud a
    %   candidate is but how far its rhythm stands out of its own noise.
    %
    %   The median is the floor rather than the mean, since the mean is pulled up
    %   by the very peak being measured.
    %
    %   The periodogram is lightly smoothed first, so that a single stray bin
    %   cannot win. Five bins is enough to join a peak that leakage has split in
    %   two and short enough not to flatten a real one.
    %
    %   A band too narrow to hold a peak, or one whose floor is zero, scores zero:
    %   such a candidate is never chosen over one that has something to show.
    %
    % Inputs:
    %   x       - navigator waveform
    %   sf      - sampling frequency in Hz
    %   bandBpm - [lo hi] rate band in beats or breaths per minute
    %
    % Output:
    %   p - peak height as a multiple of the band floor, 0 when nothing is
    %       measurable

    x = x(:) - mean(x(:));
    n = numel(x);

    if n < 16
        p = 0;
        return
    end

    pwr = abs(fft(x)) .^ 2 / n;
    freq = (0:n - 1) * sf * 60 / n;                  % bpm
    inBand = find(freq >= bandBpm(1) & freq <= bandBpm(2) & freq < sf * 30);

    if numel(inBand) < 5
        p = 0;
        return
    end

    pwr = movmean(pwr(inBand), 5);
    floorLevel = median(pwr);

    if floorLevel <= 0
        p = 0;
        return
    end

    p = max(pwr) / floorLevel;

end % bandProminence
