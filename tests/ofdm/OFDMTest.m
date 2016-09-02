% Add path
addpath('../../txrx');

% Simulation parameters
nbTest = [0 5 4 3]';        % number of bits in ADC and DAC
snrTest = (0:5:30)';        % SNR (dB)
nt = length(nbTest);        
nsnr = length(snrTest);
savedat = 1;    % Save data


% Run the simulation
snrEq = zeros(nsnr,nt);
for isnr = 1:nsnr
    snr = snrTest(isnr);
    
    for it = 1:nt
        nb = nbTest(it);
        
        % Create a simulation object
        sim = OFDMSim();
        sim.set('nbadc',nb,'nbdac',nb,'snr',snr);
        
        % Run the simulation
        sim.run();        
        
        % Save the results
        fprintf(1,'nb=%d snr=%12.4e snrEq=%12.4e\n', nb, snr, sim.snrEq);
        snrEq(isnr,it) = sim.snrEq;
    end
end

% Save results
if savedat
    save OFDMTestRes snrTest snrEq nbTest;
end

