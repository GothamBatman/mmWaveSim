classdef TxFiltIIR < hgsetget
    % Transmit IIR-based filter with D/A
    
    properties
        % architecture
        % 'feedbackk':  quantizer is in feedback loop
        % 'feedfwd':  quantizer is prior to an analog filter
        arch = 'feedfwd';
        
        nfilt = 5;      % number of filter taps
        afilt, bfilt;   % filter coefficients
        nbdac = 0;      % number of bits in the DAC
                
        % Filter specifications used in design
        nov = 2;            % oversampling ratio
        pbFreq = 0.5;       % passband freq
        sbFreq = 0.6;       % stopband freq
        sbAttn = 30;        % stopband attenuation in dB
        pbRip = 2;          % passband ripple in dB
        
        % DAC calibration parameters
        nscal = 1e4;  % number of samples used for calibration
        mseQ;         % quantizer MSE
        aq;           % optimizer quantizer scale value
        autoscale = true;   % autoscale input prior to filtering
        
        % Limit cycle detection filter
        lclen = 16;         % number of samples in detecting limit cycles
        lcprob = 1e-5;      % detection probability
        lct;                % limit cycle threshold
    end
    
    % Design the filter
    methods
        % Filter design.  This should be called before running
        function designFilt(obj)
            
            % Use Chebyshev to find
            [obj.bfilt,obj.afilt] = cheby1(obj.nfilt,obj.pbRip,obj.pbFreq);
            
            % Optimize quantizer scale value
            if obj.nbdac > 0
                % Find the optimal scale level for the quantizer assuming
                % the input has variance 1
                [aqopt,obj.mseQ] = Quant.optScale(obj.nbdac,obj.nscal);
                obj.mseQ = 10*log10(obj.mseQ);
                
                if strcmp(obj.arch,'feedfwd')
                    obj.aq = aqopt;
                elseif strcmp(obj.arch, 'feedbk')
                    
                    % Run the filter at infinite resolution to estimate the
                    % output standard deviation
                    x = randn(obj.nscal/obj.nov,1);
                    x = upsample(x,obj.nov);
                    y = filter(obj.bfilt,obj.afilt,x);
                    ystd = sqrt(mean(abs(y).^2));
                    obj.aq = aqopt/ystd;
                    
                    % Compute limit cycle detection level
                    ymag1 = filter(ones(obj.lclen,1),1,abs(y));
                    ymag = sort(ymag1,'descend');
                    obj.lct = ymag(max(round(obj.nscal*obj.lcprob),1));
                else
                    error('Unknown architecture %s', obj.arch);
                end
            end
        end
        
        % Run the filter in a batch manner
        function y = filt(obj,x)                        
            
            % Get dimensions (should we make this two dim)
            [ntin,nx] = size(x);
            ntout = ntin * obj.nov;
            na = length(obj.afilt);
            nb = length(obj.bfilt);
            
            % Rescale input to unit norm
            % Note that for complex inputs the rms value of each dimension
            % is set to unit variance.
            if (obj.autoscale)
                xstd = sqrt(mean(abs(x).^2));
                if isreal(x)
                    x = x / xstd;
                else
                    x = sqrt(2)*x / xstd;
                end
            end
            
            % Upsample x
            if (obj.nov > 1)
                x = upsample(x,obj.nov);            
            end
                            
            % Run infinite precision filter if desired
            if (obj.nbdac==0)
                y = filter(obj.bfilt,obj.afilt,x);
                return
            end
            
            % For feedforward filter, quantize and run the filter            
            if strcmp(obj.arch,'feedfwd')
                x = Quant.qsat(x,obj.nbdac,obj.aq);
                y = filter(obj.bfilt,obj.afilt,x);
                return
            end

            
            % Run the interpolation filter
            y = zeros(ntout,nx);
            for t=max(na,nb):ntout
                vp = x(t:-1:t-nb+1,:);
                yp = y(t-1:-1:t-na+1,:);
                yi = -obj.afilt(2:na)*yp + obj.bfilt*vp;
                y(t,:) = Quant.qsat(yi,obj.nbdac,obj.aq);
                
                % To remove limit cycles, test if magnitude of past samples
                % exceeds a threshold.  In that case, set output to zero.
                if (t > obj.lclen)
                    if (sum(abs(y(t-obj.lclen+1:t)))>obj.lct)
                        y(t-obj.lclen+1:t) = 0;
                    end
                end
                
            end
            
        end
    end
    
    
    
end


