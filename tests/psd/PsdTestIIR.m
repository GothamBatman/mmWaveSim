%% IIR based low-resolution filter design

% Add paths
addpath('../../txrx');

% Parameters
plotFilt = 0;       % plot the filter response
nsamp = 2^15;       % num samples to test
nbitsQTest = [0 4 5]; % number of quantizer bits to test.
pbFreq = 0.5;       % passband freq
sbFreq = 0.6;       % stopband freq
sbAttn = 30;        % stopband attenuation in dB
pbRip = 2;          % passband ripple in dB
nfft = 1024;        % num points in the FFT
dither = 0;         % use dithering to remove limit cycles

% Limit cycle detection filter
lclen = 16;         % number of samples in detecting limit cycles
lcprob = 1e-5;      % detection probability


%% Main simulation loop

% Initialize vectors
ntest = length(nbitsQTest);
Pyy = zeros(nfft/2+1,2);
mseQ = zeros(ntest,1);
legStr = cell(ntest,1);

for it = 1:ntest
    
    % Print progress
    nb = nbitsQTest(it);
    fprintf(1,'Test %d of %d: nbits=%d\n', it, ntest, nb);
    
    % Construct the filter object
    txfilt = TxFiltIIR();    
    txfilt.set('nbdac',nb);
    
    % Design the filter and optimize the scale parameters
    txfilt.designFilt();
    
    % Plot the frequency response
    if plotFilt
        npts = 256;
        [hfilt,w] = freqz(txfilt.bfilt,txfilt.afilt,npts);
        f = w/2/pi;
        plot(f,10*log10(abs(hfilt).^2));
        grid on;
        axis([0 0.5 -60 0]);
    end
    
    % Generate a random input
    x = randn(nsamp/2,1);
    
    % Filter the data
    y = txfilt.filt(x);
    
    % Compute PSD
    [Pyyi,w] = pwelch(y,[],[],nfft);
    
    % Normalize so that the passband frequencies are 0 dB.
    f = w/2/pi;
    I = (f < pbFreq/2);
    Plow = mean(Pyyi(I));
    Pyyi = Pyyi/Plow;
    
    % Save the PSD
    Pyy(:,it) = Pyyi;
    
    % Create the legend string
    if (nb==0)
        legStr{it} = 'infinite';
    else
        legStr{it} = sprintf('nb=%d',nb);
    end
end

% Plot the PSD
plot(f,10*log10(Pyy));
axis([0 0.5 -70 10]);
set(gca,'FontSize',16);
legend(legStr,'Location','SouthWest');
grid on;








