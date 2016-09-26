Settle = '01-Jan-2000';
Maturity = '01-Jan-2003';
Basis = 0;
Principal = 100;
LegRate = [0.06 0]; % [CouponRate Spread]
LegType = [1 0]; % [Fixed Float]
LegReset = [1 1]; % Payments once per year
[Price, SwapRate, AI, RecCF, RecCFDates, PayCF,PayCFDates]  = swapbyzero(ZeroRateSpec, LegRate, Settle, Maturity,...
LegReset, Basis, Principal, LegType)
load deriv.mat;

%2.example
Rate = 0.06;
Compounding  = 2;
ValuationDate = 'Jan-1-2010';
EndDates =   'Jan-1-2020';
Basis = 2;
RateSpec = intenvset('ValuationDate', ValuationDate,'StartDates', ValuationDate, ...
'EndDates', EndDates, 'Rates', Rate, 'Compounding', Compounding, 'Basis', Basis);

Settle = 'Jan-1-2011';
ExerciseDates = 'Jan-1-2012';
Maturity = 'Jan-1-2013';
Reset = 2;
Principal = 100;
Strike = 0.06;
Volatility = 0.2;
OptSpec = 'call';
Price= swaptionbyblk(RateSpec, OptSpec, Strike, Settle, ExerciseDates, Maturity, ...
Volatility, 'Reset', Reset, 'Principal', Principal, 'Basis', Basis)

y=2*1;
xl=zeros(y,1);                 %preallocation of the storage
for n=0:(y-1)
xl(n+1) = (1+Rate/2)^(-(3+n)); %save every iteration result
end
A=(1/Reset)*sum(xl)
d1=(log(Rate/Strike)+Volatility^2*(1/2))/(Volatility*sqrt(1))
d2=d1-Volatility*sqrt(1)
sprice=Principal*A*(Rate*normcdf(d1)-Strike*normcdf(d2))


%3.example
%First set:
StartDates = '01-Oct-2011';
EndDates = ['01-Oct-2012'; '01-Oct-2013'];
Rates = [0.05;0.2];
Basis=0;
RateSpec = intenvset('Rates', Rates, 'StartDates',StartDates,...
'EndDates', EndDates, 'Compounding', 1, 'Basis', Basis)
Settle = 'Oct-1-2011';
ExerciseDates = 'Oct-1-2012';
Maturity = 'Oct-1-2013';
Principal = 100;
LegRate = [NaN 0];
LegReset = [1 1]; % Payments once per year
[~, Strike] = swapbyzero(RateSpec, LegRate, Settle, Maturity,...
'StartDate', ExerciseDates, 'LegReset', LegReset)
Volatility = 0.2;
OptSpec = 'call';
Price= swaptionbyblk(RateSpec, OptSpec, Strike, Settle, ExerciseDates, Maturity, ...
Volatility, 'Principal', Principal, 'Basis', Basis)
%3.example
%Second set:
A=(1+0.2)^(-2)
d1=(log(Strike/Strike)+Volatility^2*(1/2))/(Volatility*sqrt(1))
d2=d1-Volatility*sqrt(1)
sprice=Principal*A*(Strike*normcdf(d1)-Strike*normcdf(d2))


%4.example
Principal = 100;
Volatility=0.2570
Reset=2
Ske=0.0425
Rate1=zr(7)+0.25*(zr(8)-zr(7));
Rate2=zr(7)+0.5*(zr(8)-zr(7));
Rate3=zr(7)+0.75*(zr(8)-zr(7));
Rate4=zr(8);
Rate5=zr(8)+0.5*(zr(9)-zr(8));
Rate6=zr(9);
Rate7=zr(9)+0.5*(zr(10)-zr(9));
Rate8=zr(10);
Rate9=zr(10)+0.5*(zr(11)-zr(10));
Rate10=zr(11);

A=(1/Reset)*((1+Rate1/2)^(-11)+(1+Rate2/2)^(-12)+(1+Rate3/2)^(-13)+(1+Rate4/2)^(-14)+(1+Rate5/2)^(-15)+...
(1+Rate6/2)^(-16)+(1+Rate7/2)^(-17)+(1+Rate8/2)^(-18)+(1+Rate9/2)^(-19)+(1+Rate10/2)^(-20));
d1=(log(Ske/Ske)+Volatility^2*(5/2))/(Volatility*sqrt(5))
d2=d1-Volatility*sqrt(5)
sprice=Principal*A*(Ske*normcdf(d1)-Ske*normcdf(d2))


%5.example
Principal = 100;
Volatility=SwaptionBlackVol(9,5);
Ske=SS(9,5);
R0=zr(7) %2019
R1=zr(7)+0.5*(zr(8)-zr(7)); %2020
R2=zr(8) %2021
R3=zr(9) %2022
R4=zr(10) %2023
R5=zr(11) %2024
%Ratee1=zr(3)+0.5*(zr(4)-zr(3))
Ab=(1+R1/2)^(-12)+(1+R2/2)^(-14)+(1+R3/2)^(-16)+(1+R4/2)^(-18)+(1+R5/2)^(-20)
swara=((1+R0/2)^(-10)-(1+R5/2)^(-20))/Ab
d1b=(log(swara/Ske)+Volatility^2*(5/2))/(Volatility*sqrt(5))
d2b=d1b-Volatility*sqrt(5)
sprice=Principal*Ab*(swara*normcdf(d1b)-Ske*normcdf(d2b))

%SwapRate = (Disc(1) - Disc(end))./ A;
%SBP(5,1)
%Ae(5,1)
%dd1(5,1)
%dde(5,1)
%sr(5,1)



%Code for dates, to have them in strings in "matrix"
BB=EurMatFull
for i=1:length(ExerciseDates)
    for j=1:length(Tenors)
        AA(j,i)=cellstr(datestr(BB(j,i)))
    end
end
    