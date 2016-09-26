function [dva_value] = dva(R_P, discNEE, DefProb, DefProb2)
%DVA:
% copy from above
% subtractable
SP_P = [1-DefProb];
dva_value = (1-R_P) * sum(discNEE(2:end,:) .* diff(DefProb2) .* SP_P(1:(end-1)));
end

