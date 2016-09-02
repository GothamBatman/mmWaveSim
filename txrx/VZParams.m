classdef VZParams < hgsetget
    properties
        nfft = 2048;        % number of FFT points
        delfkHz = 75;       % sub-carrier spacing
        fsMHz;              % sample rate in MHz       
        ncp0 = 160;         % num CP samples for symbol 0
        ncp1 = 144;         % num CP samples for symbols 1,...,6
        Tsym0us;            % symbol periods
        Tsym1us;
        Tslotus;            % slot time in us
        nsymslot0 = 1;      % number of symbols of type 0 in slot
        nsymslot1 = 6;      % number of symbols of type 1 in slot
        nscTot = 1200;      % total number of occupied subcarriers
        WTotMHz;            % total occupied bandwidth
    end
    
    methods
        
        % Constructor.  Derives parameters
        function obj = VZParams()
            obj.fsMHz = obj.delfkHz * obj.nfft / 1000;
            obj.Tsym0us = (obj.ncp0 + obj.nfft)/ obj.fsMHz;
            obj.Tsym1us = (obj.ncp1 + obj.nfft)/ obj.fsMHz;
            obj.Tslotus = obj.Tsym0us*obj.nsymslot0  + obj.Tsym1us*obj.nsymslot1;        
            obj.WTotMHz = obj.delfkHz * obj.nscTot / 1000;
        end
    end
end