% Load data
load OFDMTestRes;

% Create legend string
nt = length(nbTest);
legStr = cell(nt+1,1);
legStr{1} = 'Ideal';
for it =1:nt
    nb = nbTest(it);
    if nb == 0
        legStr{it+1} = 'nb=infinite';
    else
        legStr{it+1} = sprintf('nb=%d', nb);
    end
end

% Plot results
plot(snrTest,snrTest,'--', snrTest, snrEq, '-o', 'LineWidth', 2);
grid on;
set(gca,'FontSize',16);
xlabel('Input SNR');
ylabel('SNR after eq');
legend(legStr,'Location','SouthEast');

print -dpng SnrQuant