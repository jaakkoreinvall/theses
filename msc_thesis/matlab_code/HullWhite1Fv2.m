classdef HullWhite1Fv2
%HULLWHITE1F Create a Hull White 1 Factor Model
%
% Syntax:
%
%   OBJ = HullWhite1F(ZeroCurve,alpha,sigma)
%
% Description:
%
%   Create a 1 factor Hull White model by specifying the zero curve, alpha
%   and sigma parameters for the following equation
%
%   dr = [theta(t) - a*r]dt + sigma*dz
%
% Properties:
%
%   ZeroCurve - IRDataCurve or RateSpec. This is the zero curve that is
%               used to evolve the path of future interest rates.
%
%   Alpha - Mean reversion; specified either as a scalar or as a function
%               handle which takes time as as input and returns a scalar
%               mean reversion value.
%
%   Sigma - Volatility; specified either as a scalar or as a function
%               handle which takes time as as input and returns a scalar
%               mean volatility.
%
% Methods:
%
%   [ZeroRates, ForwardRates] = simTermStructs(nPeriods,'name1','val1')
%
% Example:
%
%   Settle = datenum('15-Dec-2007');
%   CurveTimes = [1:5 7 10 20]';
%   ZeroRates = [.01 .018 .024 .029 .033 .034 .035 .034]';
%   CurveDates = daysadd(Settle,360*CurveTimes,1);
%
%   irdc = IRDataCurve('Zero',Settle,CurveDates,ZeroRates);
%   
%   alpha = .1;
%   sigma = .01;
%
%   HW1F = HullWhite1F(irdc,alpha,sigma);
% 
%   SimPaths = HW1F.simTermStructs(10,'nTrials',100);
%
% References:
%
%   [1] Brigo, D and F. Mercurio. Interest Rate Models - Theory and
%       Practice. Springer Finance, 2006.
%
%   [2] Hull, J. Options, Futures, and Other Derivatives. Prentice-Hall, 2011.  
%
% See also LIBORMARKETMODEL, LINEARGAUSSIAN2F

% Copyright 1999-2012 The MathWorks, Inc.

    properties
        ZeroCurve
        Alpha
        Sigma
    end
    properties (Access = private)
        SDE
        PM
    end
    
    methods (Access = public)
        function obj = HullWhite1Fv2(inCurve,inAlpha,inSigma)
        %HULLWHITE1F Create a Hull White 1 Factor Model
            
            narginchk(3,3);
            
            % If RateSpec, convert to be an IRDataCurve
            if isafin(inCurve,'RateSpec')
                obj.ZeroCurve = IRDataCurve('Zero',inCurve.ValuationDate,inCurve.EndDates,...
                    inCurve.Rates,'Basis',inCurve.Basis,'Compounding',inCurve.Compounding);
            elseif isa(inCurve,'IRDataCurve')
                obj.ZeroCurve = inCurve;
            else
                error(message('fininst:HullWhite1Fv2:invalidCurve'));
            end
            
            % If scalar, convert to be a function handle
            if isa(inAlpha,'function_handle')
                obj.Alpha = @(t,V) inAlpha(t);
            elseif isscalar(inAlpha)
                if inAlpha == 0
                    % if alpha is 0 set it to be eps to resolve numerical
                    % issues related to it appearing in denominator
                    obj.Alpha = @(t,V) eps;
                else
                    obj.Alpha = @(t,V) inAlpha;
                end
            else
                error(message('fininst:HullWhite1Fv2:invalidAlpha'));
            end
            
            % If scalar, convert to be a function handle
            if isa(inSigma,'function_handle')
                obj.Sigma = @(t,V) inSigma(t);
            elseif isscalar(inSigma)
                obj.Sigma = @(t,V) inSigma;
            else
                error(message('fininst:HullWhite1Fv2:invalidSigma'));
            end
            
            obj.PM = @(t) obj.ZeroCurve.getDiscountFactors(datemnth(obj.ZeroCurve.Settle,t*12))';
            
            obj.SDE = hwv(@(t,X) obj.Alpha(t),0,@(t,X) obj.Sigma(t),'StartState',0);
            
        end
        function [ZeroRates, ForwardRates] = simTermStructs(obj,nPeriods,varargin)
        %SIMTERMSTRUCTS Simulate Term Structures
        %
        % Syntax:
        %
        %   Paths = simTermStructs(obj,nPeriods)
        %   Paths = simTermStructs(obj,nPeriods)
        %
        % Description:
        %
        %   Simulate future zero curve paths using the specified 1 Factor
        %   Hull White model
        %
        % Required Input Arguments:
        %
        %   nPeriods - Number of simulation periods
        %
        % Option Input Arguments:
        %
        %   deltaTime - scalar time step betwen periods. Default is 1.
        %
        %   nTrials - scalar number of trials. Default is 1.
        %
        %   antithetic - Boolean scalar flag indicating whether antithetic
        %                sampling is used to generate the Gaussian random
        %                variates that drive the zero-drift, unit-variance
        %                rate Brownian vector dW(t). See
        %                hwv/simBySolution for more information.
        %
        %   Z - Direct specification of the dependent random noise process
        %       used to generate the zero-drift, unit-variance rate Brownian
        %       vector dW(t) that drives the simulation. See
        %       hwv/simBySolution for more information.
        %
        %   Tenor - numeric vector of maturities to be computed at each
        %           time step. Default is the tenor of the object's zero
        %           curve.
        % 
        % Output Arguments:
        %
        %   ZeroRates - nPeriods+1 X nTenors X nTrials matrix of simulated
        %               zero rate term structures.
        %
        %   ForwardRates - nPeriods+1 X nTenors X nTrials matrix of simulated
        %               forward rate term structures.
        %
        % Example:
        %
        %   CurveTimes = [1:5 7 10 20]';
        %   ZeroRates = [.01 .018 .024 .029 .033 .034 .035 .034]';
        %   
        %   alpha = .1;
        %   sigma = .01;
        %
        %   HW1F = HullWhite1F([CurveTimes ZeroRates],alpha,sigma);
        % 
        %   SimPaths = HW1F.simTermStructs(10,'nTrials',100);
            
            narginchk(2,12);
            
            p = inputParser;
            
            p.addParamValue('ntrials',1);
            p.addParamValue('deltatime',1);
            p.addParamValue('tenor',[]);
            p.addParamValue('antithetic',false);
            p.addParamValue('Z',[]);
            
            try
                p.parse(varargin{:});
            catch ME
                newMsg = message('fininst:HullWhite1Fv2:optionalInputError');
                newME = MException(newMsg.Identifier, getString(newMsg));
                newME = addCause(newME, ME);
                throw(newME)
            end
            
            nTrials = p.Results.ntrials;
            deltaTime = p.Results.deltatime;
            Tenor = p.Results.tenor;
            Antithetic = p.Results.antithetic;
            Z = p.Results.Z;
            
            if isempty(Tenor)
                Tenor = yearfrac(obj.ZeroCurve.Settle,obj.ZeroCurve.Dates,obj.ZeroCurve.Basis);
                InitZeroRates = obj.ZeroCurve.Data;
                InitForwardRates = obj.ZeroCurve.getForwardRates(obj.ZeroCurve.Dates);
            else
                TenorDates = daysadd(obj.ZeroCurve.Settle,round(360*Tenor),1);
                InitZeroRates = obj.ZeroCurve.getZeroRates(TenorDates);
                InitForwardRates = obj.ZeroCurve.getForwardRates(TenorDates);
            end
            
            Tenor = reshape(Tenor,1,length(Tenor));
            
            % Generate factors and short rates
            [Paths,SimTimes] = obj.SDE.simBySolution(nPeriods,'NTRIALS',nTrials,...
                'DeltaTime',deltaTime,'antithetic',Antithetic,'Z',Z);
            
            nTenors = length(Tenor);
            
            % Allocate interest rate paths
            ZeroRates = zeros(nPeriods+1,nTenors,nTrials);
            ForwardRates = ZeroRates;
            
            ZeroRates(1,:,:) = repmat(InitZeroRates',[1 1 nTrials]);
            ForwardRates(1,:,:) = repmat(InitForwardRates',[1 1 nTrials]);
            
            % Formula for V
            V = @(t,T) obj.Sigma(t)^2/obj.Alpha(t)^2*(T - t + ...
                2/obj.Alpha(t)*exp(-obj.Alpha(t)*(T-t)) - ...
                1/(2*obj.Alpha(t))*exp(-2*obj.Alpha(t)*(T-t)) - 3/2/obj.Alpha(t));
            
            A = @(t,T) obj.PM(T)./obj.PM(t) .*exp(1/2*(V(t,T) - V(0,T) + V(0,t)));
            
            B = @(t,T) (1 - exp(-obj.Alpha(t)*(T-t)))./obj.Alpha(t);
            
            for iPeriod=2:nPeriods+1
                t = SimTimes(iPeriod);
                DiscRates = bsxfun(@times,A(t,t+Tenor),...
                    exp(-bsxfun(@times,B(t,t+Tenor),Paths(iPeriod,:,:))));
                ZeroRates(iPeriod,:,:) = bsxfun(@rdivide,-log(DiscRates),Tenor);
                ForwardRates(iPeriod,:,:) = bsxfun(@rdivide,-log(cat(2,DiscRates(1,1,:),...
                    DiscRates(1,2:end,:)./DiscRates(1,1:end-1,:))),[Tenor(1) diff(Tenor)]);
            end
        end
    end
end