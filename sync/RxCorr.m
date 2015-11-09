classdef RxCorr < hgsetget
    % RxCorr:  Receiver correlator class
    properties
        xf0;            % Transmitted freq-domain symbols
        ncoh=32;        % number of symbol of coherent integration
        nfreq=9;        % number of frequency offsets to test
        freqErrMaxkHz=50;   % frequency offset max
        fsampMHz=750;   % sampling frequency in MHz
    end
    
    methods
        % Constructor
        function obj = RxCorr()
        end
        
        % Estimate the PDP
        % pdpT is the 2D matrix of the PDPs for each frequency offset
        % pdp is the column of pdpT with the highest peak
        function [pdp,pdpT] = computePDP(obj,yt)
            
            % Num freq samples in FFT
            nfft = length(obj.xf0);    
            
            % Compute frequency offsets to test
            delMax = 2*pi*obj.freqErrMaxkHz*1e-3*nfft/obj.fsampMHz;
            del = linspace(-delMax,delMax,obj.nfreq)';
            
            % Number of periods over which to non-coherently integrate
            nnoncoh = floor(length(yt)/nfft/obj.ncoh);
            
            % Other dimensions
            nsym = nnoncoh*obj.ncoh;   % Total number of symbols
            ny = nsym*nfft;            % Total number of samples
            yt = reshape(yt(1:ny),nfft,nsym);  % Reshape y
            
            % Compute phase rotations
            % Apply D(ifreq,i) on symbol i for frequency offset ifreq
            D = exp(1i*del*(0:obj.ncoh-1));
            
            % FFT correlation
            zf = fft(yt).*repmat(conj(obj.xf0),1,nsym);
            pdpT = zeros(nfft, obj.nfreq);
            for iper = 1:nnoncoh
                
                % Get symbols to perform the coherent integration over
                I = ((iper-1)*obj.ncoh+1:iper*obj.ncoh);
                
                % Loop over frequency offsets
                for ifreq = 1:obj.nfreq
                    % Derotate the symbols within the coherent period
                    zsum = sum(zf(:,I).*repmat(D(ifreq,:),nfft,1),2);
                    
                    % Add the power non-coherently
                    pdpT(:,ifreq) = pdpT(:,ifreq) + abs(ifft(zsum)).^2;
                end
            end
            
            % Find PDP with highest frequency offset
            zmax1 = max(pdpT);
            [~,im] = max(zmax1);
            pdp = pdpT(:,im);                       
        end
        
    end
    
end

