%COLLATERALIZED AND UNCOLLATERALIZED EXPOSURE DATES SEPARATELY

%ColExpDates = daysadd(sd,360*DeltaTime*(0:nPeriods),1); % not used
%m = 20 % margin period of risk in days or months - not used
%UncolExpDates = daysadd(ColExpDates(2:end), -m, 1) % not used
%SimDates = union(UncolExpDates,ColExpDates) % para ok - not used


% SOMETHING ABOUT DATES

DeltaTime = 10/360;
SimDates=dategrid(sd, im, DeltaTime);
SimTimes = diff(yearfrac(SimDates(1),SimDates));
nPeriods = length(SimDates)-1

%Tenor = ((1/4):(1/4):10)'; %laskettavat maturiteetit jokaiselle simuloidulle pvm:lle
%Tenor = [1/12 Tenor']'

Step= round(360*DeltaTime)
IntegerTenors = ceil(Step:Step:Step*nPeriods)'
Tenor = IntegerTenors/360 % ACT/360 mukaan


thedates = daysadd(sd,(0:360*DeltaTime:daysbetween),1);

nPeriods = 5;
DeltaTime = 1;
nTrials = 1000;

Tenor = (1:10)';

SimDates = daysadd(Settle,360*DeltaTime*(0:nPeriods),1)


% reverse diag
A = randi(10,5,4)
s = size(A,1)
A(s:s-1:end-1)

s = size(HW1FSimPaths,1)
IntRates = HW1FSimPaths(exRow:s-1:s*(endCol-1)+1,:)