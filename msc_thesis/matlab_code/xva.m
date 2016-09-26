% CASES
% 1. No collateral
% 2. 1-way collateral - counterparty
% 3. 1-way collateral - self
% 4. 2-way collateral - threshold = 0
% 5. 2-way collateral - threshold = 1000

% BASIC BCVA DATA
R_C=0.4;
R_P=0.4;

%CVA CALCULATION
CVA = cva(R_C, discEE, DefProbP, DefProbCP) 
CVAc = cva(R_C, discEE_C, DefProbP, DefProbCP)
CVAp = cva(R_C, discEE_P, DefProbP, DefProbCP)
CVAcpt0 = cva(R_C, discEE_CPT0, DefProbP, DefProbCP)
CVAcpt1000 = cva(R_C, discEE_CPT1000, DefProbP, DefProbCP)

%DVA CALCULATION
DVA = dva(R_P, discNEE, DefProbP, DefProbCP)
DVAc = dva(R_P, discNEE_C, DefProbP, DefProbCP)
DVAp = dva(R_P, discNEE_P, DefProbP, DefProbCP)
DVAcpt0 = dva(R_P, discEE_CPT0, DefProbP, DefProbCP)
DVAcpt1000 = dva(R_P, discEE_CPT1000, DefProbP, DefProbCP)

% XVA plotting
% plotXVAs(CVA, DVA, CVAc, DVAc, CVAp, DVAp, CVAcpt0, ...
%    DVAcpt0, CVAcpt1000, DVAcpt1000)



