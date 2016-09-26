% CDSD BASIC DATA FOR COUNTERPARTY
CDS_settle = datenum('29-Sep-2014');
CDSDates = daysadd(CDS_settle,360*([0.5 1 2 3 4 5 7 10]),1);
ZeroData = [RateSpec.EndDates RateSpec.Rates];

% COUNTERPARTY SPECIFIC DATA
CDSSpreadsCP = [33.040 36.890 51.820 67.010 83.720 94.000 113.020 130.020]'; % CDSSpreads should include other counterparties also if there are several. Commerzbank.
CDSSpreadsP = [8.010 14.067 26.213 38.260 50.715 63.170 70.572 81.620]'; % Nordea.

% Calibrate default probabilities for each counterparty. Change only
% CDSSpreads accordingly.
DefProbCP = defprob(CDS_settle, CDSDates, CDSSpreadsCP, ZeroData, SimDates);
DefProbP = defprob(CDS_settle, CDSDates, CDSSpreadsP, ZeroData, SimDates);

% We plot of the cumulative probability of default for each counterparty.
plotdefprob(DefProbCP, SimDates)
plotdefprob(DefProbP, SimDates)