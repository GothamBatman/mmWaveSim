%% PDPthreshOpt:  Optimizes threshold selection for PDP estimation

%% Parameters
addpath('../common');

nt = 100;       % number of trials
EsN0 = -35;     % Sample SNR in dB
dlynsMax = 50;          % maximum delay
freqErrHzMax = 50e3;    % maximum frequency error
savedat = 1;            % save results

% Cluster parameters represented as a vector with one
% component per cluster
angc = 0;     % center angle in radians
angspd = 0;   % angular spread in radians
dlycns = 0;   % excess delay of first path in cluster nsec
dlyspdns = 0; % delay spread within cluster in nsec
powcdB = 0;   % gain power in dB
fadec = 0;    % fadec(i)=1 => cluster i is non-fading

% Carrier and sample rate
fsampMHz = 750;     % sample frequency in MHz
fcGHz = 28;         % carrier freq in GHz

% Phase noise parameters
addPN = 1;          % indicates if phase noise is to be added
pndBcHz = -80;      % dBc/Hz at freq pnFreqkHz
pnFreqkHz = 100;

%% Construct test class and set parameters
% Contruct test class
sim = PDPSim();


% Set the multipath parameters
chan = sim.chan;
chan.set('fsampMHz', fsampMHz, 'dlycns', dlycns, 'fadec', fadec);
chan.set('EsN0', EsN0, 'fcGHz', fcGHz);

% Set the noise and phase noise parameters
chan.set('addPN', addPN, 'pndBcHz', pndBcHz, 'pnFreqkHz', pnFreqkHz);


%% Main simulation loop

% Initialize vectors
nfft = sim.nfft;
newdat = 1;
if (newdat)
    % Generate random frequency and delay offsets
    dly0 = rand(nt,1)*dlynsMax;
    freqErr = (2*rand(nt,1)-1)*freqErrHzMax;

    % Initialize vectors
    pdpTot = zeros(nfft,nt);
    dlyEst = zeros(nt,1);     
    powEst = zeros(nt,1);
    
    for it = 1:nt
        
        % Set parameters
        chan.set('dlycns', dly0(it), 'freqErrHz', freqErr(it) );
        
        % Run the test
        pdpi = sim.run();
        
        % Extract the PDP
        pdpTot(:,it) = pdpi;
        
        % Get the delay estimate
        [pdpMax,im] = max(pdpi);
        dlyEst(it) = im/fsampMHz*1e3;
        powEst(it) = pdpMax;
        
        % Print progress
        fprintf(1,'it=%d dly error=%12.4e\n', it, dly0(it)-dlyEst(it));
        
    end
    
end


%% Optimize the threshold
% Find noise level in each run
pdpSort = sort(pdpTot);
mufrac = 0.5;                    % fraction of samples to measure noise level
n1 = round(mufrac*nfft);
noisePow = mean(pdpSort(1:n1,:));
sigPowMax = max(pdpSort);
muMax = min(sigPowMax./noisePow);

% Set threshold level
mutest = linspace(0.25*muMax,0.9*muMax,100);     % Threshold levels to test
nmu = length(mutest);


dlySamp = (0:nfft-1)'*1e3/fsampMHz;
dlyErrTot = zeros(nt,nmu);

for imu = 1:nmu
    
    % Find the threshold level
    mu = mutest(imu);
    t = mu*noisePow;
    
    % Compute the mean and std deviation of
    I = (pdpTot > repmat(t,nfft,1));
    d1 = mean( I.*pdpTot.*repmat(dlySamp,1,nt) );
    d0 = mean( I.*pdpTot );
    dlyEsti = d1./d0;
    dlyErrTot(:,imu) = dlyEsti'-dly0;
    
end

% Find threshold with optimal delay
dlyErrAvg = std(dlyErrTot)';
[dlyErrMin, im] = min(dlyErrAvg);
muOpt = mutest(im);

% Compute power
t = muOpt*noisePow;
I = (pdpTot > repmat(t,nfft,1));
powEst = sum(I.*pdpTot)';
powErr = std( 10*log10(powEst) );

p = 10*log10(powEst);
p = p - mean(p);
plot(sort(p),(1:nt)/nt,'-','LineWidth', 2);
grid on;


