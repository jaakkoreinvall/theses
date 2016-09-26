function [ColExp, ColNexp] = collateralisation(exposures, nexposures, m, H, w, varargin)
%COLLATERALISATION Summary of this function goes here
%   Detailed explanation goes here
% exposures: NumDates-by-NumCounterparties-by-NumScenarios "cube"
% m: is allowed to be margin period of risk and meaning period between
% collateral calculation and collateral posting.
% H: threshold value
    nTrials = size(exposures, 3);
    SwaptionLength = size(exposures, 1);
    Threshold = H*ones((SwaptionLength-m),1,nTrials);
    if w == 1
        Col = max(exposures(1:(end-m),1,:)-Threshold,0); % collateral calculation
        ColTotal = cat(1,zeros(1:m,1,nTrials),Col); % leave m first as is, i.e., uncollateralized
        ColExp = max(exposures-ColTotal,0); % collateralized exposures
        ColNexp = max(nexposures+ColTotal,0);
    elseif w == 2
        Col = max(exposures(1:(end-m),1,:)-Threshold,0); % collateral calculation
        ColTotal = cat(1,zeros(1:m,1,nTrials),Col); % leave m first as is, i.e., uncollateralized
        ColCP = max(nexposures(1:(end-m),1,:)-Threshold,0); % collateral calculation
        ColCPTotal = cat(1,zeros(1:m,1,nTrials),ColCP); % leave m first as is, i.e., uncollateralized
        ColExp = max(exposures-ColTotal,0); % collateralized exposures taking into account CP posting collateral
        ColExp = max(ColExp+ColCPTotal,0); % collateralized exposures taking into account P posting collateral
        ColNexp = max(nexposures-ColCPTotal,0); % collateralized exposures taking into account P posting collateral
        ColNexp = max(ColNexp+ColTotal,0); % collateralized exposures taking into account CP posting collateral
    else
        disp('Value is neither 1 nor 2')
    end
        
%0)  NEXT STEPS: 
%* Create collateralisation part
%* Create minimalistic results section, for example one graph for CVA and DVA, where different collateralisation arrangements have been studied
%* This gives new insight into what how much work is still left to be done.

%1) RESULTS SECTION NEXT STEPS:
%* check Nik for style
%* also a little bit of other writing, but not too much have to start
%* get 2 versions done per day

end

