function [DefProb] = defprob(CDS_settle, CDSDates, CDSSpreads, ZeroData, SimDates)
% CDS info:
% Day count: ACT/360
% Frequency: Q
% Recovery rate: 0.4
% Reco Override: 0.4
% Red info update time: 29/9
% Unit: bps
% CDS-type sub or senoior?
% Import CDS market information for each counterparty
% CDS = readtable(swapFile,'Sheet','CDS Spreads');
% disp(CDS);


% Calibrate default probabilities for each counterparty

DefProb = zeros(length(SimDates), size(CDSSpreads,2));
% i goes for number of counterparties
for i = 1:size(DefProb,2)
    probData = cdsbootstrap(ZeroData, [CDSDates CDSSpreads(:,i)],...
        CDS_settle, 'probDates', SimDates);
    DefProb(:,i) = probData(:,2);
end
end

