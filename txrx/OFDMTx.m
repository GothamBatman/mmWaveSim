classdef OFDMTx < hgsetget
    properties
        % Default parameters, taken from VZ spec
        nfft = 2048;    % number of FFT 
        nsc = 1200;     % number of used sub-carriers
        ncp = 144;      % number of CP
    end
    
    methods
        % Constructor
        function obj = OFDMTx()            
        end
        
        % Modulate symbols
        function xtd = mod(obj,xsym)
            
            % Compute number of symbols and zero pad symbols
            nx = size(xsym,1);             
            
            % Zero pad modulation symbol to an even number of symbols
            nsym = ceil(nx/obj.nsc);
            xsym = [xsym; zeros(nsym*obj.nsc-nx,1)];
            xsym = reshape(xsym,obj.nsc,nsym);
                        
            % Map to used subcarriers 
            xf1 = zeros(obj.nfft,nsym);
            xf1(1:obj.nsc/2,:) = xsym(obj.nsc/2+1:obj.nsc,:);
            xf1(obj.nfft-obj.nsc/2+1:obj.nfft,:) = ...
                xsym(1:obj.nsc/2,:);
            
            % Take IFFT
            xtd = ifft(xf1);
            
            % Add CP
            xtd = [xtd; xtd(1:obj.ncp,:)];                       
            xtd = xtd(:);
        end
        
         % Compute the time-domain and frequency-domain preamble
        function [xpref,xpret] = getPre(obj,x,n0,npre)
            i0 = (n0-1)*obj.nsc+1;
            i1 = (n0+npre-1)*obj.nsc;
            xpref = reshape(x(i0:i1),obj.nsc,npre);
            xpret = obj.mod(xpref(:));
        end
           
    end
end