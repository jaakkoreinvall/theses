function [discEE] = discexpexp(exposures,nTrials,DF,SimDates)
% We compute the discounted expected exposures using the discount factors from each simulated interest rate scenario. 
% Get discounted exposures per counterparty, for each scenario

discExp = zeros(size(exposures));
for i = 1:nTrials
    discExp(:,:,i) = bsxfun(@times,[1 diag(DF(:,:,i),-1)']',exposures(:,:,i)); %diag(DF(:,:,i),-1)  includes discount factors to the very beginning.
end

% Discounted expected exposure
discProfiles = exposureprofiles(SimDates,discExp,'ProfileSpec','EE');

% Aggregate the discounted EE for each counterparty into a matrix
discEE = [discProfiles.EE];

end

