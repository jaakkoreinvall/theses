%Constructing zero curve

%Settlement date. The date when zero curve was retrieved from market data
sd=datenum('23-Sep-2014');

%Curve dates daysadd: Date = daysadd(D1, NUM, BASIS i.e. day count basis).
%Zero rates are here percentages. Eg. 0.1198 % and 0.2075 %.
%Q1: Which day count basis to use? Use 2 = actual/360 at least for now. Not
%chosen yet in daysadd
%So assume now 365 days per year. 3m and 6m must be integer, so fixed that way. 
%365*[91 183 1 2 3 4 5 7 8 9 10 15 20]
%Step years taken into account below.
cd=daysadd(sd,([91 183 365 2*365+1 3*365+1 4*365+1 5*365+1 7*365+2 8*365+2 9*365+2 10*365+3 15*365+4 20*365+5]),2);
zr=[0.1198	0.2075	0.4411	0.9088	1.3527	1.7481	2.0902	2.6094	2.8097	2.9902	3.1541	3.7152	4.0508]'/100;
plot(cd,zr)
datetick
title(['Zero curve for ' datestr(sd)])
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
ied=datenum('23-Sep-2019');

% Instrument maturity. This is the maturity of the instrument if the
% swaption holder enters into swap. This one is 5 years after the expiration date.
im=datenum('23-Sep-2024');

% Instrument strike price [bps]. As for the strike, swaption vols are quoted ATM. 
% So the strike for each of these is the forward par swap rate at the option expiry 
% date for the tenor specified on the contract. This means that every
% element is its own swaption and it has its own strike price (i.e. in this
% case forward swap rate). Those swap rates can be calculated from zero
% rates curve.

%Q2: Is bid the same as buy? Yes bid is buy, ask is sell.
%Q3: "is" means instrument strike. Where to derive it? 
is=0.045;
%Q4: are swap rates different for different tenors?

%Q5: Which one to use for swaption black vols: ask/mid/bid?
SwaptionBlackVol = [48.71	50.91	41.81	34.32	30.69	28.24	26.35	25.27	24.09	22.81
57.00	51.47	41.73	35.59	32.13	29.75	27.67	26.40	25.22	24.11
57.08	48.61	40.29	35.25	32.34	30.20	28.40	27.35	26.12	24.98
55.96	46.84	39.60	35.43	32.29	30.27	28.67	27.57	26.47	25.47
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

%Below dates chosen a bit ambigiously
%ExerciseDates = [1:10]
ExerciseDates = [365 2*365+1 3*365+1 4*365+1 5*365+1 6*365+2 7*365+2 8*365+2 9*365+2 10*365+3]
%Tenors =[(1/12) (1/4) (1/2) (3/4) 1:10]
%Think about automating date process. KARKAUSVUODET!
Tenors = [30 91 181 273 365 2*365 3*365+1 4*365+1 5*365+1 6*365+1 7*365+2 8*365+2 9*365+2 10*365+2]

% B = repmat(A,M,N), creates a large matrix B consisting of an M-by-N tiling of copies of A
% EurExDatesFull are exercise dates. 14 rows corresponding to tenors.
% Exercise dates change according to columns.
EurExDatesFull = repmat(daysadd(sd, ExerciseDates, 2)', length(Tenors),1);

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

EurMatFull = reshape(daysadd(EurExDatesFull, repmat(Tenors,1,length(ExerciseDates)),2),size(EurExDatesFull));

%CORRECTIONS:
EurMatFull(3,1)=EurMatFull(3,1)+1;
EurMatFull(3,5)=EurMatFull(3,5)+1;
EurMatFull(3,9)=EurMatFull(3,9)+1;
EurMatFull(4,1)=EurMatFull(4,1)+1;
EurMatFull(4,5)=EurMatFull(4,5)+1;
EurMatFull(4,9)=EurMatFull(4,9)+1;
EurMatFull(5,1)=EurMatFull(5,1)+1;
EurMatFull(5,5)=EurMatFull(5,5)+1;
EurMatFull(5,9)=EurMatFull(5,9)+1;
EurMatFull(6,1)=EurMatFull(6,1)+1;
EurMatFull(6,4)=EurMatFull(6,4)+1;
EurMatFull(6,5)=EurMatFull(6,5)+1;
EurMatFull(6,8)=EurMatFull(6,8)+1;
EurMatFull(6,9)=EurMatFull(6,9)+1;
EurMatFull(7,2)=EurMatFull(7,2)-1;
EurMatFull(7,6)=EurMatFull(7,6)-1;
EurMatFull(7,10)=EurMatFull(7,10)-1;
EurMatFull(9,1)=EurMatFull(9,1)+1;
EurMatFull(9,5)=EurMatFull(9,5)+1;
EurMatFull(9,9)=EurMatFull(9,9)+1;
EurMatFull(10,1)=EurMatFull(10,1)+1;
EurMatFull(10,4)=EurMatFull(10,4)+1;
EurMatFull(10,5)=EurMatFull(10,5)+1;
EurMatFull(10,8)=EurMatFull(10,8)+1;
EurMatFull(10,9)=EurMatFull(10,9)+1;
EurMatFull(11,2)=EurMatFull(11,2)-1;
EurMatFull(11,6)=EurMatFull(11,6)-1;
EurMatFull(11,10)=EurMatFull(11,10)-1;
EurMatFull(13,1)=EurMatFull(13,1)+1;
EurMatFull(13,5)=EurMatFull(13,5)+1;
EurMatFull(13,9)=EurMatFull(13,9)+1;
EurMatFull(14,1)=EurMatFull(14,1)+1;
EurMatFull(14,4)=EurMatFull(14,4)+1;
EurMatFull(14,5)=EurMatFull(14,5)+1;
EurMatFull(14,8)=EurMatFull(14,8)+1;
EurMatFull(14,9)=EurMatFull(14,9)+1;

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
%OTHER UNUSED:
%strk = zeros(size(SwaptionBlackVol));
%sr = zeros(size(SwaptionBlackVol));
%Ae= zeros(size(SwaptionBlackVol));
%dd1= zeros(size(SwaptionBlackVol));
%dde=zeros(size(SwaptionBlackVol));

% Good to know: NINST is number of instruments

% To understand strike here. Read BSC thesis. []


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
        [~,SS(iTenor,iSwaption)] = swapbyzero(RateSpec,[NaN 0], sd, EurMatFull(iTenor,iSwaption),...
            'StartDate',EurExDatesFull(iTenor,iSwaption),'LegReset',[1 1]);
        
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

nPeriods = 5;
DeltaTime = 1;
nTrials = 1000;

Tenor = (1:10)';

SimDates = daysadd(sd,360*DeltaTime*(0:nPeriods),1)
SimTimes = diff(yearfrac(SimDates(1),SimDates))

% For 1 year periods and an evenly spaced tenor, the exercise row will be
% the sixth row and the swaption maturity will be the 5th column
exRow = 6;
endCol = 5;

%Simulate Interest-Rate Paths Using the Hull-White One-Factor Model
TimeSpec = hwtimespec(sd,daysadd(sd,360*(1:11),1), 2); %

%2nd and 4th arguments in hwvolspec can be anything no change in the values
%of the output parameters.
HW1Fobjfun = @(x) SBP(relidx) - ...
    swaptionbyhw(hwtree(hwvolspec(sd,'11-Aug-2015',x(2),'11-Aug-2015',x(1)), RateSpec, TimeSpec), 'call', SS(relidx),...
    EurExDatesFull(relidx), 0, sd, EurMatFull(relidx));
options = optimset('disp','iter','MaxFunEvals',1000,'TolFun',1e-5);
% Find the parameters that minimize the difference between the observed and
% predicted prices
x0 = [.1 .01];
lb = [0 0];
ub = [1 1];
HW1Fparams = lsqnonlin(HW1Fobjfun,x0,lb,ub,options);
HW_alpha = HW1Fparams(1)
HW_sigma = HW1Fparams(2)

