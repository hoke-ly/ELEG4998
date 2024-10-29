% Radix-2 FFT algorithms simply dividing the DFT 
% of length N into two DFTs of length N/2 each, and iterating, 
% until N smaller or equal to 1. 
% Decimation-in-time (DIT) means the process 
% is applied to the time-domain samples (input data).

function X = myFFT(x)
    N = length(x);
    
    if N <= 1
        X = x;
        return;
    end
    
    % zero-padding to adjust the length to a power of 2
    nextPow2 = 2^nextpow2(N);
    x = [x, zeros(1, nextPow2 - N)]; 

    % divide and conquer
    even = myFFT(x(1:2:end));
    odd = myFFT(x(2:2:end));

    X = zeros(1, N);
    for k = 0:(N/2 - 1)
        % Twiddle factor for phase shift on odd DFTs
        t = exp(-2*pi*1i*k/N) * odd(k+1); 
        % Positive Frequencies
        X(k+1) = even(k+1) + t;
        % Negative Frequencies
        X(k+1+N/2) = even(k+1) - t;
    end
end
