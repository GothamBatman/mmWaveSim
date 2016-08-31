% Add paths
addpath('../../txrx');

% Simulation parameters
nsym = 2^11;
nbdac = 4;
nbadc = 4;
nfiltTx = 3;    % TX filter order
snr = 15;       % SNR
nsyminit = 3;   % symbol index of the first symbol in the preamble
nsympre = 16;    % number of symbols used for initial sync
nsync = 2048;   % number of samples used in initial synchronization

% OFDM parameters.  Get these from the VZ spec.
phyp = VZParams();
nfft = phyp.nfft;
ncp = phyp.ncp1;
nsc = phyp.nscTot;
nsampsym = ncp+nfft;

% Create OFDM transmitter
ofdmtx = OFDMTx();
ofdmtx.set('nfft',nfft,'nsc',nsc,'ncp',ncp);

% Create modulation symbols
x = exp(1i*pi/2*(randi(4,nsc*nsym,1)+0.5));

% OFDM modulation
xtd = ofdmtx.mod(x);

% Compute the preamble time-domain and freq-domain
[xpref,xpret] = ofdmtx.getPre(x,nsyminit,nsympre);

% TX filter
txfilt = TxFiltIIR();
nov = 1;                    % oversampling ratio
pbFreq = nsc/nfft/nov;      % passband freq
sbFreq = 1.1*pbFreq;        % stopband freq
txfilt.set('nbdac',nbdac,'pbFreq',pbFreq,'sbFreq',sbFreq,'nov',nov, ...
    'nfilt',nfiltTx);
txfilt.designFilt();

% Run the TX filter
y0 = txfilt.filt(xtd);

if 0
npts = 1024;
[Pyy,w] = pwelch(y0,[],[],npts);
f = w/2/pi;
plot(f,10*log10(fftshift(Pyy)));
grid on;
end

% Add noise
snr1 = snr - 10*log10(nfft/nsc);
ns = size(y0,1);
yvar = mean(abs(y0).^2);
wvar = 10^(-0.1*snr1)/2*yvar;
y = y0 + sqrt(wvar)*(randn(ns,1) + 1i*randn(ns,1));

% Design the RX filter
rxfilt = RxFilt();
rxfilt.set('nbadc',nbadc,'fc',pbFreq,'nov',nov);
rxfilt.designFilt();

% Run the RX filter
r = rxfilt.filt(y);

% Create the OFDM RX object
ofdmrx = OFDMRx();
ofdmrx.set('nfft',nfft,'nsc',nsc,'ncp',ncp);

% Set the synchronization parameters
ipre0 = (nsyminit-1)*nsampsym;      % first sample to start in preamble search
npresearch = 2^12;             % max delay range
ofdmrx.set('ipre0',ipre0,'npresearch',npresearch,'xpret',xpret,...
    'xpref',xpref);

% Get the initial delay and channel estimate
ofdmrx.sync(r);

% Demodulate to create unequalized symbols
rs = ofdmrx.demod(r);

% Get equalized symbols
xeq = ofdmrx.equalize(rs);

% Compute MSE
xtr = x((nsyminit-1)*nsc+1:(nsym-1)*nsc);
mse = mean(abs(xtr-xeq(:)).^2);
snrEq = -10*log10(mse);

fprintf(1,'Max SNR=%f Act SNR=%f\n', snr, snrEq);


