function [cva_value] = cva(R_C, discEE, DefProb, DefProb2)
%CVA COMPUTATION
SP_P = [1-DefProb2];
cva_value = (1-R_C) * sum(discEE(2:end,:) .* diff(DefProb) .* SP_P(1:(end-1)));
%for i = 1:numel(CVA)
%    fprintf('CVA for counterparty %d = $%.2f\n',i,CVA(i));
%end
end

