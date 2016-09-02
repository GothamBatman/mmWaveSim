classdef OFDMSim < hgsetget
    % Test to evaluate front end by TX and RX OFDM modulation symbols and
    % then measuring the SNR
    properties
        % Simulation parameters
        nsym = 2^11;
        nbdac = 4;
        nbadc = 4;
        nfiltTx = 3;    % TX filter order
        nfiltRx = 3;    % RX filter order
        snr = 15;       % SNR
        nsyminit = 3;   % symbol index of the first symbol in the preamble
        nsympre = 16;   % number of symbols used for initial sync
        nsync = 2048;   % number of samples used in initial synchronization
        
        % OFDM parameters
        nfft;       % number of FFT points
        ncp;        %
        nsc;
        
        % Output of the test
        y0;         % TX waveform
        snrEq;      % Post-equalized symbol SNR
    end
    
    methods
        % Constructor
        function obj = OFDMSim()
            % OFDM parameters.  Get these from the VZ spec.
            phyp = VZParams();
            obj.nfft = phyp.nfft;
            obj.ncp = phyp.ncp1;
            obj.nsc = phyp.nscTot;
        end
        
        % Run test
        function run(obj)
            % Add paths
            addpath('../../txrx');

            % Compute dimensions
            nsampsym = obj.ncp+obj.nfft;                                   

            % Create OFDM transmitter
            ofdmtx = OFDMTx();
            ofdmtx.set('nfft',obj.nfft,'nsc',obj.nsc,'ncp',obj.ncp);
            
            % Create modulation symbols
            x = exp(1i*pi/2*(randi(4,obj.nsc*obj.nsym,1)+0.5));
            
            % OFDM modulation
            xtd = ofdmtx.mod(x);
            
            % Compute the preamble time-domain and freq-domain
            [xpref,xpret] = ofdmtx.getPre(x,obj.nsyminit,obj.nsympre);
            
            % TX filter
            txfilt = TxFiltIIR();
            nov = 1;                    % oversampling ratio
            pbFreq = obj.nsc/obj.nfft/nov;      % passband freq
            sbFreq = 1.1*pbFreq;        % stopband freq
            txfilt.set('nbdac',obj.nbdac,'pbFreq',pbFreq,'sbFreq',sbFreq,'nov',nov, ...
                'nfilt',obj.nfiltTx);
            txfilt.designFilt();
            
            % Run the TX filter
            obj.y0 = txfilt.filt(xtd);
            
            if 0
                npts = 1024;
                [Pyy,w] = pwelch(obj.y0,[],[],npts);
                f = w/2/pi;
                plot(f,10*log10(fftshift(Pyy)));
                grid on;
            end
            
            % Add noise
            snr1 = obj.snr - 10*log10(obj.nfft/obj.nsc);
            ns = size(obj.y0,1);
            yvar = mean(abs(obj.y0).^2);
            wvar = 10^(-0.1*snr1)/2*yvar;
            y = obj.y0 + sqrt(wvar)*(randn(ns,1) + 1i*randn(ns,1));
            
            % Design the RX filter
            rxfilt = RxFilt();
            rxfilt.set('nbadc',obj.nbadc,'fc',pbFreq,'nov',nov);
            rxfilt.designFilt();
            
            % Run the RX filter
            r = rxfilt.filt(y);
            
            % Create the OFDM RX object
            ofdmrx = OFDMRx();
            ofdmrx.set('nfft',obj.nfft,'nsc',obj.nsc,'ncp',obj.ncp);
            
            % Set the synchronization parameters
            ipre0 = (obj.nsyminit-1)*nsampsym;      % first sample to start in preamble search
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
            xtr = x((obj.nsyminit-1)*obj.nsc+1:(obj.nsym-1)*obj.nsc);
            mse = mean(abs(xtr-xeq(:)).^2);
            obj.snrEq = -10*log10(mse);
            
            %fprintf(1,'Max SNR=%f Act SNR=%f\n', obj.snr, obj.snrEq);
        end
    end
end


