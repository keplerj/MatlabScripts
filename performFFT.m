function [Fss, Magnitude, Phase] = performFFT(y, dqrate, window)
% PERFORMFFT performs n-point fft on input signal y of given frequency
%            rate dqrate applying a specified window.
%
%   (c) Georgios Printezis 2022

    NFFT = height(y);

    %% Apply Window    
    if window == "hanning"
        W = hann(NFFT);
    elseif window == "hanning2"
        W = hann(NFFT) .^ 2;
    elseif window == "flattop"
        W = flattopwin(NFFT);
    elseif window == "rectangular"
        W = ones(NFFT, 1);
    end

    y = W .* y / sum(W);

    %% Two-Sided
    Yts = fft(y, NFFT, 1);
    
    Fts = (0 : 1 / NFFT : 1 - 1 / NFFT)' * dqrate;

    %% Single-Sided
    Yss = Yts(1 : round(NFFT / 2), :);
    
    Yss(2 : end - 1, :) = 2 * Yss(2 : end - 1, :);
    
    Fss = Fts(1 : round(NFFT / 2));

    Magnitude = abs(Yss);        % Magnitude of the FFT
    
    Phase = unwrap(angle(Yss));  % Phase of the FFT
end
