%% PDPSim:  Simulates one run of the PDP computation
classdef PDPSim < hgsetget
    properties
        
        %% Simulation parameters
        % Parameters
        nfft = 1024;        % num FFT points
        nsym = 2^10;        % number of symbols
        EsN0 = -35;         % Sample SNR in dB
        sampppm = 0;        % Sample clock error in ppm
        loppm = 0;          % LO error in ppm
        nfftT = 32;         % number of FFT symbols over time.
        
        % Simulation objects
        chan;       % Channel object
        rx;         % Receiver
        
    end
    
    methods
        % Constructor
        function obj = PDPSim()
            
            % Construct the simulation objects
            obj.chan = MPChan();
            obj.rx = RxCorr();
        end
        
        % Run test
        function pdp = run(obj)
            
            %% Transmitter
            
            % Generate random QPSK symbols
            b = randi(4,obj.nfft,1)-0.5;
            xf0 = exp(2*pi*1i*b/4);
            
            xt0 = ifft(xf0);             % Take IFFT
            xt = repmat(xt0,obj.nsym,1); % Repeat
            xt = xt(:);
            
            %% Multipath channel                       
            
            % Generate the subpaths
            obj.chan.genSubPath();
            
            % Filter the data through the multipath channel
            yt = obj.chan.filt(xt);
            
            %% Receiver            
            obj.rx.set('xf0', xf0);
            pdp = obj.rx.computePDP(yt);
        end
    end
end


