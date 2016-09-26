clear

%Constructing zero curve

%Settlement date. The date when zero curve was retrieved from market data
sd=datenum('30-Sep-2014');

%Curve dates daysadd: Date = daysadd(D1, NUM, BASIS i.e. day count basis).
%Zero rates are here percentages. Eg. 0.1198 % and 0.2075 %.
%Q1: Which day count basis to use? Use 2 = actual/360 at least for now.
%Step years taken into account below.
%3M,6M,1Y,2Y,3Y,4Y,5Y,7Y,8Y,9Y,10Y,15Y,20Y
%cd=daysadd(sd,([91 181 365 2*365+1 3*365+1 4*365+1 5*365+1 7*365+2 8*365+2 9*365+2 10*365+3 15*365+4 20*365+5]),2);
%zr=[0.1198 0.2075 0.4411 0.9088	1.3527 1.7481 2.0902 2.6094 2.8097 2.9902 3.1541 3.7152 4.0508]'/100;
%cd=daysadd(sd, 360*([1:20]),1);
%zr = [0.1482	0.617	1.1203	1.5246	1.8308	2.0638	2.2458	2.3922	2.5131	2.6152	2.7028	2.7791	2.846	2.9053	2.9582	3.0055	3.0482	3.0867	3.1216	3.1534]'/100;
cd=daysadd(sd, 360*([1 2 3 4 5 7 10 12 15 20]),1);
zr=[0.2000	0.6374	1.1058	1.4708	1.7341	2.1078	2.4576	2.6227	2.7912	2.949]'/100;

plot(cd,zr)
xlabel('Maturity (year)')
ylabel('Interest rate (%)')
datetick
% interest rate curve represented with data
% IRDataCurve parameters: 1-type of the curve, 2-settle, 3-dates
% corresponding to rates data, 4-zero rates themselves. 
irc=IRDataCurve('zero',sd,cd,zr);
% RateSpec is the annualized zero rate term structure. It's confirmed that
% it includes annualized discount factors with assumed semi-annual
% compounding. Confirmed at least somewhat that us treasury bonds compound
% semi-annually.
% RateSpec.Disc = (1 + Z/F)^(-T), where F is the compounding frequency, Z is the zero rate, and T is the time in 
% periodic units (in this case it is 2*half-year); for example, T = F is 1 year.

RateSpec=intenvset('Rates',zr,'EndDates',cd,'StartDate',sd,'Basis',2,'EndMonthRule',1,'Compounding',2);

% Defining swaption parameters. The price of a European swaption is computed
% with an exercise date of 5 years and underlying swap of 5 years.

% Instrument exercise date. This is 5 years from the retrieval data
% date.
ied=datenum('30-Sep-2019');

% Instrument maturity. This is the maturity of the instrument if the
% swaption holder enters into swap. This one is 5 years after the expiration date.
im=datenum('30-Sep-2024');

% Instrument strike price [bps]. As for the strike, swaption vols are quoted ATM. 
% So the strike for each of these is the forward par swap rate at the option expiry 
% date for the tenor specified on the contract. This means that every
% element is its own swaption and it has its own strike price (i.e. in this
% case forward swap rate). Those swap rates can be calculated from zero
% rates curve.

%Q2: Is bid the same as buy? Yes bid is buy, ask is sell.
%Q3: "is" means instrument strike. Where to derive it? 
%is=0.045;
%Q4: are swap rates different for different tenors?

%Q5: Which one to use for swaption black vols: ask/mid/bid? Mid used ATM
%{
SwaptionBlackVol = [57.00	51.47	41.73	35.59	32.13	29.75	27.67	26.40	25.22	24.11
57.08	48.61	40.29	35.25	32.34	30.20	28.40	27.35	26.12	24.98
53.76	44.17	37.97	34.37	32.05	30.25	28.75	27.85	26.79	25.82
41.29	35.92	33.08	31.07	29.37	28.31	27.41	26.75	26.05	25.15
35.25	32.57	30.48	29.00	27.72	26.94	26.32	25.82	25.26	24.64
31.82	30.07	28.81	27.67	26.65	26.03	25.39	25.01	24.50	23.86
29.78	28.23	27.48	26.47	25.70	25.07	24.48	24.14	23.77	23.36
28.21	26.86	26.15	25.27	24.59	24.07	23.62	23.36	23.05	22.63
26.72	25.56	24.88	24.11	23.54	23.13	22.81	22.62	22.38	21.93
25.57	24.51	23.86	23.22	22.70	22.35	22.05	21.92	21.69	21.25
24.44	23.46	22.89	22.36	21.90	21.59	21.33	21.23	21.01	20.59
23.28	22.45	21.95	21.53	21.12	20.89	20.64	20.57	20.34	19.95]/100;
%}

SwaptionBlackVol = [52.87	44.12	38.82	35.52	33.38	31.40	30.01	27.85	26.79	26.77
41.25	36.58	34.11	32.23	30.74	29.59	28.56	26.76	26.06	26.14
35.84	33.90	31.72	30.19	28.92	28.04	27.24	25.83	25.26	25.42
32.88	31.06	29.77	28.71	27.66	26.91	26.29	25.02	24.51	24.72
30.45	28.87	28.01	27.22	26.44	25.80	25.27	24.15	23.78	24.06
28.93	27.62	26.78	26.04	25.35	24.82	24.37	23.37	23.06	23.37
27.54	26.45	25.60	24.90	24.33	23.88	23.51	22.63	22.39	22.71
26.34	25.35	24.59	23.98	23.47	23.07	22.73	21.92	21.69	21.99
25.17	24.27	23.62	23.10	22.65	22.30	21.99	21.24	21.01	21.29
23.97	23.24	22.70	22.25	21.87	21.57	21.27	20.58	20.35	20.61]/100;


%Below dates chosen a bit ambigiously
ExerciseDates = [1:10];
Tenors = [1:10];

%Plotting the surface of Black volatilities
figure
surf(Tenors,ExerciseDates,(100*SwaptionBlackVol)) 
title(['USD ATM European swaption Black volatilities (bps)'])
xlabel('Tenor (Years)')
ylabel('Expiry (Years)')
ax = gca
ax.XTick = [1:10]
ax.YTick = [1:10]

% B = repmat(A,M,N), creates a large matrix B consisting of an M-by-N tiling of copies of A
% EurExDatesFull are exercise dates. 12 (originally 14) rows corresponding to tenors.
% Exercise dates change according to columns.
EurExDatesFull = repmat(daysadd(sd, 360*ExerciseDates, 1)', length(Tenors),1);

% EurMatFull are maturity dates
% repmat(Tenors,1,length(ExerciseDates)) describes matrix in which there
% are 10 (=no. exercise dates) times the 14 tenors in one row.
% daysadd(EurExDatesFull, repmat(Tenors,1,length(ExerciseDates)),2): adds
% days so that tenor lengths are taken into account. We get final dates for
% swaps. Code works in a way that tenor lenghts are added to exercise dates
% one-by-one.
% reshape(X,M,N) or reshape(X,[M,N]) returns the M-by-N matrix whose elements are taken columnwise from X. 
% An error results if X does not have M*N elements. This function makes the
% column back to its original matrix form which is the same as for
% EurExDatesFull, 14 X 10.

EurMatFull = reshape(daysadd(EurExDatesFull, repmat(360*Tenors,1,length(ExerciseDates)),1),size(EurExDatesFull));



% Swaption maturities in readable form
swaptionmats=tomatrix(EurMatFull,10,10);
% Select Calibration Instruments
% Find the swaptions that expire on or before the maturity date of the
% sample swaption.  In this case, all swaptions having an underlying tenor that 
% matures before the maturity of the swaption to be priced (23-Sep-2024) are used 
% in the calibration. So can use for example 6 X 3 swaption or 3 X 6
% swaption.

% Find indices that correspond to the rule.
% Confirmed that this works ok.
relidx = find(EurMatFull <= im);

% Compute Swaption Prices Using Black's Model

% SBP for SwaptionBlackPrices
% SS for SwaptionStrike 
SBP = zeros(size(SwaptionBlackVol));
SS = zeros(size(SwaptionBlackVol));
CFd1 = cell(size(SwaptionBlackVol));
CFd2 = cell(size(SwaptionBlackVol));
CF1 = cell(size(SwaptionBlackVol));
CF2 = cell(size(SwaptionBlackVol));
% Good to know: NINST is number of instruments

for iSwaption=1:length(ExerciseDates) % over swaption exercisedates (= 10)
    for iTenor=1:length(Tenors) % over swap tenors/maturities (= 14)
        % When ~ appears in the left side of an assignment statement with multiple outputs, it means that
        % that particular output should be discarded
        % swapbyzero prices a vanilla swap by a set of zero curve(s).(I think that floating rate comes from forward
        % interests, which can be obtained from zero curve, at least it has been checked that swapbyzero computes 
        % correctly). Or more like pricing strike swap rate values, which
        % go to swaptionbyblk in Strike.
        % RateSpec is here used to discount both paying and receiving legs.
        % ATTENTION: [NaN 0] as legrate specify that swap rate is solved so
        % that swap is priced to zero. Also, spread is zero here.
        % LegReset=[1 1] means payments once per year.
        % I am confident this calculates correct strike swap rate as in
        % Glasserman
        [~,SS(iTenor,iSwaption),~,CF1{iTenor,iSwaption},CFd1{iTenor,iSwaption},CFd2{iTenor,iSwaption},CFd2{iTenor,iSwaption}] = swapbyzero(RateSpec,[NaN 0], sd, EurMatFull(iTenor,iSwaption),...
            'StartDate',EurExDatesFull(iTenor,iSwaption),'LegReset',[4 4]);
        
        % swaptionbyblk: prices swaptions using the Black option pricing model.
        % Strike is third input in the function, hence, it comes from the
        % above function.
        % A 'call swaption' entitles the buyer to pay the fixed rate. A 'put swaption' entitles the buyer 
        % to receive the fixed rate. Fixed rate as currency.
        % A Call swaption or Payer swaption allows the option buyer to enter into an interest rate swap 
        % in which the buyer of the option pays the fixed rate and receives the floating rate.
        % A Put swaption or Receiver swaption allows the option buyer to enter into an interest rate swap
        % in which the buyer of the option receives the fixed rate and pays the floating rate.
        
        SBP(iTenor,iSwaption) = swaptionbyblk(RateSpec, 'call', SS(iTenor,iSwaption),sd, ...
            EurExDatesFull(iTenor,iSwaption), EurMatFull(iTenor,iSwaption), SwaptionBlackVol(iTenor,iSwaption));
        % [SBP(iTenor,iSwaption), strk(iTenor,iSwaption), sr(iTenor,iSwaption), Ae(iTenor,iSwaption), dd1(iTenor,iSwaption), dde(iTenor,iSwaption)] = sbb(RateSpec, 'call', SS(iTenor,iSwaption),sd, ...
        %    EurExDatesFull(iTenor,iSwaption), EurMatFull(iTenor,iSwaption), SwaptionBlackVol(iTenor,iSwaption));   
    end
end



is=SS(5,5)

% NEXT STEPS: 
% check that main steps (not every little detail) is ok after these changes in swaption valuation
% Create collateralisation part

nYears = 10;
DeltaTime = 10/360;
nPeriods = nYears/DeltaTime;
SimDates = daysadd(sd, 360*DeltaTime*(0:nPeriods), 1);
SimTimes = diff(yearfrac(SimDates(1),SimDates,1));
datestr(SimDates)

nTrials = 100; % number of Monte Carlo simulations
Step = 360*DeltaTime
Tenor = ceil(Step:Step:Step*nPeriods)'/360 % 30/360 SIA
Tenor = (DeltaTime:DeltaTime:DeltaTime*nPeriods)

% Use this form, because it works somehow. It is required that this multiplied by 360 is an integer number.
% Yo. tenor sisaltaa kaikki t -> t+k missa k kay lapi 1 kk to 120 kk


% For 1 year periods and an evenly spaced tenor, the exercise row will be
% the sixth row and the swaption maturity will be the 5th column
%exRow = length(SimDates)/2; % 61: monthly
%endCol = length(nPeriods/2); % 60 means 5th year, because 1/12 month intervals

exRow = find(SimDates==datenum(ied)); % ied="instrument exercise date"
endCol = exRow-1;


%Simulate Interest-Rate Paths Using the Hull-White One-Factor Model
%TimeSpec = hwtimespec(sd,EurMatFull(relidx), 2); %
TimeSpec = hwtimespec(sd,daysadd(sd,360*(1:11),1), 2);
% Defines the observation dates of the Hull-White tree and the compounding rule 
% for date to time mapping and price-yield formulas.


%2nd and 4th arguments in hwvolspec can be anything no change in the values
%of the output parameters.
%test:HWTree=hwtree(hwvolspec(sd,'11-Aug-2015',0.01,'11-Aug-2015',0.1), RateSpec, TimeSpec)
%test:treeviewer(HWTree)
HW1Fobjfun = @(x) SBP(relidx) - ...
    swaptionbyhw(hwtree(hwvolspec(sd,'23-Sep-2021',x(2),'23-Sep-2021',x(1)), RateSpec, TimeSpec), 'call', SS(relidx),...
    EurExDatesFull(relidx), 0, sd, EurMatFull(relidx)); %does not seem to include basis
options = optimset('disp','iter','MaxFunEvals',1000,'TolFun',1e-5);
% Find the parameters that minimize the difference between the observed and
% predicted prices
x0 = [.1 .01];
lb = [0 0];
ub = [1 1];
HW1Fparams = lsqnonlin(HW1Fobjfun,x0,lb,ub,options);
HW_alpha = HW1Fparams(1)
HW_sigma = HW1Fparams(2)

% CONSTRUCTING THE HULLWHITE1F MODEL USING THE HULLWHITE1F CONSTRUCTOR:

HW1F = HullWhite1Fv2(RateSpec,HW_alpha,HW_sigma) % HullWhite1F function works as follows: 
tic
HW1FSimPaths = HW1F.simTermStructs(nPeriods,'NTRIALS',nTrials,...
    'DeltaTime',DeltaTime,'Tenor',Tenor,'antithetic',true); % simulates future zero curve paths using a specified HullWhite1F object. 
toc
trialIdx = 1;
figure
surf(Tenor,SimDates,HW1FSimPaths(:,:,trialIdx)) % Plots the colored parametric surface defined by three matrix arguments. Because of trialIdx, only 1 trial shown here. SIZE=PERIODS X TENORS X TRIALS
datetick y keepticks keeplimits
title(['Evolution of the Zero Curve for Trial:' num2str(trialIdx) ' of Hull White Model'])
xlabel('Tenor (Years)')

% PRICE THE EUROPEAN SWAPTIONS:

N=100;
sigma=1/4;
CashFlowFrequency = 3/12;
CashFlowInterval = round(CashFlowFrequency/DeltaTime);
DF = exp(bsxfun(@times,-HW1FSimPaths,repmat(Tenor,[nPeriods+1 1]))); % Continuous discounting. To beginning of the swap.
Sum = sigma*sum(bsxfun(@times,1,DF(exRow,CashFlowInterval:CashFlowInterval:endCol,:)));
SwapRate = (1 - DF(exRow,endCol,:))./Sum; % Top equation at 564
PayoffValue = N*max(SwapRate-is,0).*Sum; % Discount to T_n
PayoffValueCopy = repelem(PayoffValue,1,exRow); % for swaption exposure dates
RealizedDFs = [];
for t=1:(exRow-1)
    RealizedDF = prod(exp(-HW1FSimPaths(t,endCol-(t-1),:)*sum(SimTimes(t:exRow-1))),1); % This can go wrong
    RealizedDFs = [RealizedDFs RealizedDF];
end
RealizedDFs = [RealizedDFs reshape([repelem(1,nTrials)],size(RealizedDF))];
Swaptions = RealizedDFs.*PayoffValueCopy; % Swaption values for each secenario & for each exposure date before and on exercise date.
SwaptionValues = mean(Swaptions,3);


% PRICE THE SWAPS: 
%Jos aikoo hinnoitella swaption suoraan diskontatuista kassavirroista (ei
%kaavaalla), niin tulee hieman eri arvo kuin kaavasta, koska max(x,0) tulee
%menemaan erilailla, kun arvot menevat nollan eri puolin. Pelkka swapin
%arvo on kuitenkin hyvin lahelle sama.


% INITIALATIZIONS
SwaptionLength = length(SimDates);
SwapLength = SwaptionLength - exRow;
FinalIndex = SwaptionLength; % Final index for the SimDates
FixingIndex = (exRow:CashFlowInterval:FinalIndex)'; % Fixing date indeces
%FixIn2Mths=((exRow+1):CashFlowInterval:(FinalIndex-2))' % Fixing date in 2 months time indeces
%FixIn1Mth=((exRow+2):CashFlowInterval:(FinalIndex-1))'% Fixing date in 1 month indeces
PayFixed=is; % Fixed payments
RecFloat=HW1FSimPaths(exRow:CashFlowInterval:(end-CashFlowInterval),CashFlowInterval,:); % Floating receivables
Swaps=[]; % Collect swap values here




for t=(exRow+1):(SwaptionLength-1) % from 1st real swap date to almost end
    if any(FixingIndex == t) % meaning 
        RecFloat=HW1FSimPaths(t:CashFlowInterval:(end-CashFlowInterval),CashFlowInterval,:);
        D=DF(t,CashFlowInterval:CashFlowInterval:(SwaptionLength-t),:); % ok
    else
        NextFixing = FixingIndex(find(FixingIndex>t,1));
        D=DF(t,(NextFixing-t):CashFlowInterval:(SwaptionLength-t),:);    
    %elseif any(FixIn2Mths == t)
    %    D=DF(t,2:3:(endCol-(t-exRow)),:);
    %else any(FixIn1Mth == t)
    %    D=DF(t,1:3:(endCol-(t-exRow)),:);
    end
    Float=sum(bsxfun(@times,permute(RecFloat,[2,1,3]),D),2);
    Fixed=sum(D,2)*is;
    Swaps=[Swaps (100*(1/4)*(Float-Fixed))];
end


% CORRECTION ZEROING SWAP VALUES FOR NON-EXERCISED SWAPTIONS
NonExercised = find(PayoffValue==0);
LenIndeces=length(NonExercised);
Swaps(1,:,NonExercised)=zeros(1,(SwapLength-1),LenIndeces); % SwapLength-1: because more zeros coming
SwapValues=mean(Swaps,3);

% EXPOSURES
% Remember to change exposures equation if multiple trades.
% Exposure of an unnetted contract is equal to the market value of the contract if it has positive value, otherwise it is zero.

DerivativeValues = [Swaptions Swaps zeros(1,1,nTrials)];

% Ehka kuva viela tahan peraan.

% We compute these exposures for the entire portfolio as well as each counterparty at each simulation date using the creditexposures function.
% NumDates-by-NumCounterparties-by-NumScenarios "cube" of credit exposures representing potential losses that could be incurred over all 
% dates, counterparties, and scenarios, if a counterparty defaulted (ignoring any post-default recovery).
[exposures, expcpty] = creditexposures(permute(DerivativeValues,[2 1 3]),'Commerzbank');
nexposures = permute(max(-DerivativeValues,0),[2 1 3]) % exposures from CP POV


% COLLATERAL EFFECT

% We now compute the BCVA in the collateralised case assuming a
% zero-threshold, two-way CSA. The situation is much more balanced, since 
% the impact of the CSA is to create similar amounts of residual exposure for both the institution and
% the counterparty. The impact of this is that the BCVA is approximately
% half the CVA, in line with the institution's CDS spread being half that
% of the counterparty.

% The results below rely on the assumption of a margin period of risk for
% the institution as well as the counterparty. This requires an institution
% to consider a benefit, not only from their own default, but also from
% their ability to stop posting collateral 10 days prior to their default.
% Finally, note that by signing a two-way CSA, the institution would
% realise a loss... Their counterparty would realise the same amount of
% gain.

% Collateralization parameters

% add flexibility about when collateral is posted?
% EXTRA:
% 1. check collateralized exposure definition
% 2. move stuff to the beginning in the collateralization
% 3. assumed that posting collateral is negative for a party

m = (10/360)/DeltaTime;
H = 0;
[E_C, NE_C] = collateralisation(exposures, nexposures, m, H, 1);
[NE_P, E_P] = collateralisation(nexposures, exposures, m, H, 1);
[E_CPT0, NE_CPT0] = collateralisation(exposures, nexposures, m, 0, 2);
[E_CPT1000, NE_CPT1000] = collateralisation(exposures, nexposures, m, 2, 2);


%new function
%collateralisation(exposures, nexposures, m, H, ?)


% Compute exposure profiles for each counterparty and entire portfolio
numDates = numel(SimDates);
cpProfiles = exposureprofiles(SimDates,exposures);
eprof=exposureprofiles(SimDates, E_C);


% Visualize portfolio exposure profiles
% Exposures are the same as expected exposures, of course, because only one
% contract is in question here.


figure;
plot(SimDates,cpProfiles.PFE,...
    SimDates,cpProfiles.MPFE * ones(numDates,1),...
    SimDates,cpProfiles.EE,...
    SimDates,cpProfiles.EPE * ones(numDates,1),...
    SimDates,cpProfiles.EffEE,...
    SimDates,cpProfiles.EffEPE * ones(numDates,1));
legend({'PFE (95%)','Max PFE','Exp Exposure (EE)','Time-Avg EE (EPE)',...
    'Max past EE (EffEE)','Time-Avg EffEE (EffEPE)'})
datetick('x','mmmyy')
title('Portfolio Exposure Profiles');
ylabel('Exposure ($)')
xlabel('Simulation Dates')



%{
%TAHAN JAIT (05/06/15)
PayFixed=sum(bsxfun(@times,1,DF(exRow,2:endCol,:)))*is
RecFloat=
SwapValue=

SwapValue=100*(1/4)*(sum(bsxfun(@times,permute(HW1FSimPaths(exRow:3:(end-3),2,:),[2,1,3]),DF(exRow,2:endCol,:)))-)
SwapPayoff=max(SwapValue,0)
mean(SwapPayoff)



%tarvittavat valinnat mitä cf:it ottaa huomioon (koska ero kuukausittaisen ja joka kolmannen kk välillä)

temp=(SwapRate-is).*(sigma*sum(bsxfun(@times,1,DF(exRow,2:endCol,:))))
%}

%{

%TESTI NO:2 ------------

temps=(SwapRate-is).*(sigma*sum(bsxfun(@times,1,DF(exRow,3:3:endCol,:))))
temp=mean(temps)

Float=sum(bsxfun(@times,permute(HW1FSimPaths(exRow:3:(end-3),3,:),[2,1,3]),DF(exRow,3:3:endCol,:)))
Fixed=sum(bsxfun(@times,1,DF(exRow,3:3:endCol,:)))*is
SwapValues=(1/4)*(Float-Fixed)
SwapValue=mean(SwapValues)

diff=mean(SwapValues-temps)

%TESTI----------

%The expiry time of a swaption is typically the first reset date T0 of the
%underlying interest rate swap.
%}





%{
datestr(CFd1{5,5})
HW1FSimPaths
R=is
L=
pay fixed, receive floating

%*(HW1FSimPaths(exRow+1,1,:)-is)

*R=is % Remember to change this for SS(??,??)
*Sum from i=n+1 to M+1
*B(T_n,T_i)=
*B(T_n,T_M+1)=
%}

% HW1F_SwaptionPrice = mean(RealizedDF.*PayoffValue); % Avg. over 1000 trials.

% EXPOSURE CALCULATION:
% 2 years after start
% That is 23-Sep-2016.
% Then it is discounting swap cash flows all the way from 23-Sep-2024 to 23-Sep-2020 back to 23-Sep-2019. 
% And discounting that to 23-Sep-2016.
% No change in SwapRate because exercise date is still the same and so are cashflow dates.
% No change in PayoffValue because is is strike that does not change and discounts are to T_n (exercise_date)
% Only RealizedDF changes.

%{
DF_2
PayoffValue_2
RealizedDF_2 = prod(exp(bsxfun(@times,-HW1FSimPaths(2:exRow,1,:),SimTimes(1:exRow-1))),1);
Exposures = RealizedDF_2.*PayoffValue_2
%}

% DISCOUNTED EXPECTED EXPOSURE
% Uuden discEE luominen tapahtuu s.e. vaihtaa discEE ja exposures nimet

discEE = discexpexp(exposures,nTrials,DF,SimDates)
discEE_C = discexpexp(E_C,nTrials,DF,SimDates)
discEE_P = discexpexp(E_P,nTrials,DF,SimDates)
discEE_CPT0 = discexpexp(E_CPT0,nTrials,DF,SimDates)
discEE_CPT1000 = discexpexp(E_CPT1000,nTrials,DF,SimDates)

discNEE = discexpexp(nexposures,nTrials,DF,SimDates)
discNEE_C = discexpexp(NE_C,nTrials,DF,SimDates)
discNEE_P = discexpexp(NE_P,nTrials,DF,SimDates)
discNEE_CPT0 = discexpexp(NE_CPT0,nTrials,DF,SimDates)
discNEE_CPT1000 = discexpexp(NE_CPT1000,nTrials,DF,SimDates)

% plotting

%plotdiscEE(discEE,SimDates)
%plotdiscEE(discEE_C,SimDates)


% DEFAULT PROBABILITIES COMMERZBANK:

cds


%BCVA laskeminen (Gregory, 2012, p. 267) mukaan

%Recovery values

%Dates
%j=1,...,n
%t_0=now

%Exposures
%EE_{t_{j}} is discounted expected exposure.
%EE>=0 always.
%NEE_{t_{j}} discounted expected negative exposure.
%NEE=-EE<=0

%Default probabilities
%DP_P(t_j-1,t_j)="probability that party defaults in the interval"
%DP_C(t_j-1,t_j)="probability that counterparty defaults in the interval"


%{
*EE(t_j)
*NEE(t_j)
*PD_I(0,t_j-1)
*PD_C(0,t_j-1)
*PD_I(t_j-1,t_j)
*PD_C(t_j-1,t_j)

BCVA=(1-R_C)*sum_(j=1)^n(EE(t_j)[1-PD_I(0,t_j-1)]PD_C(t_j-1,t_j)
+(1-R_P)*sum_(j=1)^n(NEE(t_j)[1-PD_C(0,t_j-1)]PD_I(t_j-1,t_j)
%}


xva
plotXVAs(CVA, DVA, CVAc, DVAc, CVAp, DVAp, CVAcpt0, ...
    DVAcpt0, CVAcpt1000, DVAcpt1000)





