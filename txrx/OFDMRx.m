classdef OFDMRx < hgsetget
    properties
        % Default parameters, taken from VZ spec
        nfft = 2048;    % number of FFT
        nsc = 1200;     % number of used sub-carriers
        ncp = 144;      % number of CP
        
        % Synchronization parameters
        maxDly = 2^13;
        ipre0;          % First sample to start in preamble search
        npresearch;     % max delay
        xpret;          % time-domain to search
        xpref;          % frequency domain preamble
        idly;  % Starting sample
        rho;   % sync correlation
    end
    
    methods
        % Constructor
        function obj = OFDMRx()
        end
        
        % Synchronization
        % Searches for preamble xpre in a received vector r
        % i0 is the first sample
        function sync(obj,r)
            
            % Extract sample to search
            npre = length(obj.xpret);
            rs = r(obj.ipre0:obj.ipre0+npre+obj.npresearch-1);
            obj.rho = xcorr(rs,obj.xpret);
            obj.rho = abs(obj.rho).^2;
            [~,im] = max(obj.rho);           
            obj.idly = im -length(rs)+obj.ipre0-1;            
        end
        
        % OFDM demodulation.  Recovers unequalized symbols
        function rs = demod(obj,r)
            
            % Strip initial estimated delay
            r = r(obj.idly:end,:);
            nr = size(r,1);
            
            % Reshape to symbols
            nsampsym = obj.ncp+obj.nfft;
            nsym = floor(nr/nsampsym);
            r = reshape(r(1:nsym*nsampsym),nsampsym,nsym);
            
            % Remove CP
            r = r(obj.ncp/2:obj.ncp/2+obj.nfft-1,:);
            
            % Take FFT and derotate
            rf1 = fft(r);
            delf = -(0:obj.nfft-1)'*obj.ncp/2/obj.nfft;
            rf1 = rf1.*repmat(exp(2*pi*1i*delf),1,nsym);
            
            % Extract the used subcarriers
            rs = zeros(obj.nsc,nsym);
            rs(obj.nsc/2+1:obj.nsc,:) = rf1(1:obj.nsc/2,:);
            rs(1:obj.nsc/2,:) = rf1(obj.nfft-obj.nsc/2+1:obj.nfft,:);
                          
        end
        
        % Create equalized modulation symbols
        function xeq = equalize(obj,rs)
            
            % Extract preambles symbols
            npre = size(obj.xpref,2);
            ns = size(rs,2);
            
            % Compute channel estimate on the symbols
            h = rs(:,1:npre)./obj.xpref;
            
            % Average over time.  Right now, there is no averaging over
            % frequency
            h = mean(h,2);
            
            % Equalize the symbols
            xeq = rs./repmat(h,1,ns);
                        
        end
    end
end