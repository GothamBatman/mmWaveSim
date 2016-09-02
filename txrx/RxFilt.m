classdef RxFilt < hgsetget
    properties
        nfilt = 64;        % filter order
        bfilt;             % num coeffs
        afilt;             % denom coeffs
        nov = 2;           % oversampling ratio
        fc = 0.5;          % cut off frequency
        nbadc;             % number of bits in the ADC
        
        % ADC calibration parameters
        nscal = 1e4;  % number of samples used for calibration
        mseQ;         % quantizer MSE
        aq;           % optimizer quantizer scale value   
        autoscale = true;  % scale input (as if there is an oracle AGC)
    end
    
    methods
        
        % Constructor
        function obj = RxFilt()
        end
        
        % Design the filter and set the quantizer values
        function designFilt(obj)
            % fir1 creates an FIR filter
            % scale normalizes the coefficients so that the magnitude
            % response of the filter at the center of the passband is 1 (0dB)            
            obj.bfilt = fir1(obj.nfilt, obj.fc, 'low', 'scale');
            obj.afilt = 1;
            
            % Optimize quantizer scale value
            % Find the optimal scale level for the quantizer assuming
            % the input has variance 1
            if obj.nbadc > 0
                [obj.aq,obj.mseQ] = Quant.optScale(obj.nbadc,obj.nscal);
                obj.mseQ = 10*log10(obj.mseQ);                          
            end
        end
        
        % Filter of input performed blockwise
        function y = filt(obj, x)
            
            % Autoscale.  This would be done by the AGC.
            % Note that in the complex case, the input is scaled
            % unit variance per I/Q.
            if obj.autoscale
                xvar = mean(abs(x).^2);
                if isreal(x)
                    x = x / sqrt(xvar);
                else
                    x = x / sqrt(xvar/2);
                end
                
            end
            
            % Quantize the input.
            xq = Quant.qsat(x,obj.nbadc,obj.aq);
            
            % Filter the quantized input
            y = filter(obj.bfilt, obj.afilt, xq);
            if (obj.nov > 1)
                y = downsample(y,obj.nov);           
            end
        end
    end
end

