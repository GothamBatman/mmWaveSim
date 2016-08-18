classdef MPChan < hgsetget
    %MPChan Multipath fading channel class
    
    properties
        nsub = 20;  % number of subpaths per cluster
        fmaxHz = 0;     % max Doppler spread in Hz
        freqErrHz = 0;  % constant frequency error in Hz
        
        % Cluster parameters represented as a vector with one
        % component per cluster
        angc = 0;     % center angle in radians
        angspd = 0;   % angular spread in radians
        dlycns = 0;   % excess delay of first path in cluster nsec
        dlyspdns = 0; % delay spread within cluster in nsec
        powcdB = 0;   % gain power in dB
        fadec = 0;    % fadec(i)=1 => cluster i is non-fading
        
        % Angles and delay of each path.
        angp;       % AoA in radians
        fdHz;       % Doppler shift in Hz = fmaxHz*cos(angp)
        gain;       % complex gain of each path
        dlypns;     % delay in us
        
        nov = 4;     % oversample ratio
        EsN0 = 10;   % SNR per sample
        
        % Carrier and sample rate
        fsampMHz = 750;     % sample frequency in MHz
        fcGHz = 76;         % carrier freq in GHz
        
        % Phase noise parameters
        addPN = 1;          % indicates if phase noise is to be added
        pndBcHz = -85;      % dBc/Hz at freq pnFreqkHz
        pnFreqkHz = 100;   
        
        
    end
    
    methods
        % Constrcutor
        function obj = MPChan()
        end
        
        % Generate the sub path coefficients
        function genSubPath(obj)
            
            % Convert scalar parameters to vectors
            nc = length(obj.angc);      % num clusters
            if length(obj.angspd) == 1
                obj.angspd = obj.angspd*ones(nc,1);
            end
            if length(obj.dlycns) == 1
                obj.dlycns = obj.dlycns*ones(nc,1);
            end
            if length(obj.dlyspdns) == 1
                obj.dlyspdns = obj.dlyspdns*ones(nc,1);
            end
            if length(obj.powcdB) == 1
                obj.powcdB = obj.powcdB*ones(nc,1);
            end
            
            % Compute total number of paths
            ncfade = sum(obj.fadec);    % num fading clusters
            np = obj.nsub*ncfade + nc-ncfade;       % total num of subpaths
            
            % Select random parameters for the sub-paths within each
            % cluster
            obj.angp = zeros(np,1);
            obj.dlypns = zeros(np,1);
            powp = zeros(np,1);     % power of each path in linear scale
            ip = 0;
            for ic = 1:nc
                % Compute number of subpaths in the cluster
                % 1 for non-fading, nsub for fading
                if obj.fadec(ic)
                    nsubi = obj.nsub;
                else
                    nsubi = 1;
                end
                I = (ip+1:ip+nsubi)';
                
                % Set the sub-path angles, delay spread
                obj.angp(I) = obj.angc(ic) ...
                    + obj.angspd(ic)*(rand(nsubi,1)-0.5);
                obj.dlypns(I) = obj.dlycns(ic)  ...
                    - obj.dlyspdns(ic)*log(rand(nsubi,1));
                powp(I) = 10^(0.1*obj.powcdB(ic))/nsubi;
                ip = ip + nsubi;
            end
            
            % Generate complex gains and Doppler shifts
            obj.gain = sqrt(powp).*exp(2*pi*1i*rand(np,1));
            obj.fdHz = obj.fmaxHz*cos(obj.angp) + obj.freqErrHz;
        end
        
        % Generate time-varying frequency response
        function H = genTimeFreqResp(obj,tms,fMHz)
            nt = length(tms);
            nf = length(fMHz);
            np = length(obj.dlypns);
            H = zeros(nt,nf);
            T = repmat(tms,1,nf)*1e-3;   % Time matrix in seconds
            F = repmat(fMHz',nt,1); % Freq matrix in MHz
            for ip=1:np
                H = H + exp(2*pi*1i*(T*obj.fdHz(ip)+obj.dlypns(ip)*F*1e-3))*obj.gain(ip);
            end
        end
        
        % Filter the signal x.  
        % The signal x should be a vector of complex baseband samples at
        % the sample rate
        function [y,y0] = filt(obj,x)
            
            % Get dimensions and initialize y
            tsms = 1e-3/obj.fsampMHz;    % sample period in ms
            nt = size(x,1);
            np = length(obj.fdHz);
            nt1 = nt*obj.nov;
            ts1 = tsms/obj.nov; % sample period of upcoverted signal in ms
            
            % Upsample
            x1 = resample(x,obj.nov,1);
            y = zeros(nt1,1);
            
            % Filter the signal by summing each path
            t = (0:nt1-1)'*ts1*1e-3;    % Time in seconds
            for ip=1:np
                idly = round(obj.dlypns(ip)*1e-6/ts1);
                xs = [zeros(idly,1); x1(1:nt1-idly)];
                y = y + exp(2*pi*1i*t*obj.fdHz(ip))*obj.gain(ip).*xs;
            end
            
            % Downsample
            y0 = resample(y,1,obj.nov);
                        
            % Add phase noise
            if (obj.addPN)
                y0 = obj.addPhaseNoise(y0);
            end
            
            % Add noise
            yvar = mean(abs(y0).^2);
            wvar = 10.^(-0.1*obj.EsN0)*yvar;
            ny = length(y0);
            w = sqrt(wvar/2)*(randn(ny,1) + 1i*randn(ny,1));
            y = y0 + w;

            
        end
        
        % Adds phase noise
        function [y1,theta] = addPhaseNoise(obj,y)
            
            % Compute phase variance per sample
            % variance in delay over a period T = c*T, 
            %   c = 10^(0.1*pndBcHz)*(f0/fc)^2
            % Hence the phase rotation in radians
            %   = sqrt(c)*fc*2*pi = 2*pi*10^(0.05*pndBcHz)*f0*sqrt(Tsamp)
            phaseStd = 2*pi*obj.pnFreqkHz*1e3 ...
                *sqrt(10^(0.1*obj.pndBcHz)*1e-6/obj.fsampMHz);
            
            % Add phase noise
            ny = length(y);
            theta = cumsum(randn(ny,1)*phaseStd);
            y1 = y.*exp(1i*theta);             
        end
        
    end
    
    
end

