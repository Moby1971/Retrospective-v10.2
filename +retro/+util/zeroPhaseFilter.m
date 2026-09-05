% Zero-phase Butterworth filter, designed in second-order sections
%
% Author : Gustav Strijkers
% Date   : 2026-09-03

function y = zeroPhaseFilter(x, order, Wn, type)

    % Purpose:
    %   Band-filter a navigator without shifting it in time, using a filter that
    %   stays numerically stable at the narrow bands the navigator needs.
    %
    % How it works:
    %   The filter is designed as zeros, poles and gain and converted to
    %   second-order sections, never to a single transfer function. That is the
    %   whole point of this function.
    %
    %   A navigator band is narrow relative to the sampling rate -- a respiratory
    %   band can be a hundredth of Nyquist -- and expanding such a design into
    %   [b, a] polynomials is ill-conditioned. The resulting coefficients no
    %   longer describe the intended filter: measured on this project's own
    %   settings, a 4th-order bandpass at TR 4 ms and 90 bpm with a 10 bpm
    %   bandwidth gives a pole at |p| = 1.0032, which is outside the unit circle
    %   and therefore unstable, and filtfilt returns values of order 1e22.
    %
    %   Those values are finite, so a check for NaN or Inf does not catch them.
    %   The filtered navigator is then numerical noise that looks like data, and
    %   every trigger taken from it is meaningless.
    %
    %   In second-order sections each pole pair is applied separately, so the
    %   conditioning problem does not arise and the same designs come out stable.
    %
    %   filtfilt runs the filter forwards and backwards, which is what makes it
    %   zero phase: any phase distortion is undone by the reverse pass. That
    %   matters more here than the magnitude response, since the whole purpose is
    %   to find when a peak occurred.
    %
    % Inputs:
    %   x     - navigator waveform, real
    %   order - Butterworth order, as passed to butter
    %   Wn    - normalised cutoff: a scalar for 'low', a two-element [lo hi] for
    %           'bandpass', each in (0, 1) where 1 is Nyquist
    %   type  - 'bandpass' or 'low'
    %
    % Output:
    %   y - the filtered waveform, same size and orientation as x
    %
    % Raises:
    %   Retrospective:filterBand - Wn is not a valid normalised band

    if any(Wn <= 0) || any(Wn >= 1) || (numel(Wn) == 2 && Wn(1) >= Wn(2))
        error('Retrospective:filterBand', ...
            'Filter band %s is not inside (0,1) with a rising edge order', mat2str(Wn, 4));
    end

    [z, p, k] = butter(order, Wn, type);
    sos = zp2sos(z, p, k);

    y = filtfilt(sos, 1, x);

end % zeroPhaseFilter
