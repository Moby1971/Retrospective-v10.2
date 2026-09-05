% Which way up a respiratory navigator sits
%
% Author : Gustav Strijkers
% Date   : 2026-09-04

function [tf, confidence] = respirationPeaksUp(x, sf, rateHz)

    % Purpose:
    %   Decide whether the breathing excursions in a navigator point upwards, so
    %   that the peak finder that follows lands on the breaths and not on the
    %   rests between them.
    %
    % How it works:
    %   Breathing is not a wave. It is a brief excursion followed by a long rest,
    %   and that asymmetry is what tells the polarity: the correct way up is the
    %   one where the signal spends LESS of its time high than low. Get it wrong
    %   and the peak finder triggers on the rests, which are longer, more numerous
    %   and in the wrong place, and the respiratory gating goes with them.
    %
    %   Two things have to be arranged before that comparison means anything, and
    %   measured over 200 constructed records apiece they matter more than the
    %   choice of statistic:
    %
    %   The signal is bandpassed around the breathing rate first, not lowpassed.
    %   A lowpass at the respiratory rate passes everything below it as well, and
    %   what lies below it is slow wander -- heating, settling, gradual drift --
    %   which is symmetric, often larger than the breathing, and swamps any
    %   statistic of the amplitude distribution. Judged on a lowpassed trace the
    %   polarity came out right 56% of the time, which is to say not at all;
    %   through a band from four tenths to three times the rate, 100%. The upper
    %   edge keeps the first two harmonics, which is what carries the asymmetry --
    %   a single harmonic is a sine and has no way up.
    %
    %   The record is then judged in windows of about ten breaths and the windows
    %   vote. A sigh, a swallow, a twitch is one large one-sided excursion, and on
    %   a whole record taken at once a handful of them decides the answer: with
    %   artefacts present the single-window figure was 72%, and the vote 100%. A
    %   vote also survives a baseline that moves between windows.
    %
    %   Within a window the statistic is the fraction of samples above the
    %   midpoint of the first and ninety-ninth percentiles. Percentiles rather
    %   than the extremes so that one sample cannot set the midpoint, and a
    %   fraction rather than a mean-median or a skewness because it states the
    %   physiology directly: brief excursions, long rests.
    %
    %   The vote is returned with its margin. A trace with no asymmetry to read --
    %   a continuously breathing animal, or one whose breathing the navigator
    %   barely sees -- gives a margin near zero, and that is the case where the
    %   operator should reach for the flip switch.
    %
    % Inputs:
    %   x      - navigator waveform, real
    %   sf     - sampling frequency in Hz
    %   rateHz - respiratory rate in Hz, from the detected or entered rate
    %
    % Output:
    %   tf         - true when the breathing excursions already point upwards
    %   confidence - 0 to 1, the margin of the vote; 0 when nothing was readable

    tf = true;
    confidence = 0;

    x = x(:);

    if numel(x) < 64 || ~isfinite(rateHz) || rateHz <= 0
        return
    end

    % A band that rejects the wander below the breathing and the heart above it,
    % while keeping the harmonics the asymmetry lives in
    band = [0.4 * rateHz, 3 * rateHz] / (sf / 2);
    band(2) = min(band(2), 0.95);

    if band(1) <= 0 || band(1) >= band(2)
        return
    end

    try
        band1 = retro.util.zeroPhaseFilter(x, 4, band, 'bandpass');
    catch
        return          % an unusable band decides nothing rather than guessing
    end

    % About ten breaths to a window, and at least a few hundred samples
    width = max(round(10 * sf / rateHz), 200);
    nrWindows = max(floor(numel(band1) / width), 1);

    votes = zeros(1, nrWindows);

    for k = 1:nrWindows

        first = (k - 1) * width + 1;
        last = min(k * width, numel(band1));
        segment = band1(first:last);

        edges = prctile(segment, [1 99]);
        midpoint = mean(edges);

        if edges(2) <= edges(1)
            votes(k) = 0;               % flat: nothing to read here
        else
            % less than half the time above the midpoint means the excursions
            % are the brief ones, which is the right way up
            votes(k) = sign(0.5 - mean(segment > midpoint));
        end

    end

    if all(votes == 0)
        return
    end

    tally = mean(votes);
    tf = tally >= 0;
    confidence = abs(tally);

end % respirationPeaksUp
